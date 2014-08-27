// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */

#include <Eigen/Dense> 
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/extensions/hscpsf/hscPsf.h"
#include "lsst/meas/extensions/hscpsf/GRMinimizer.h"

namespace afwImage = lsst::afw::image;


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


HscPcaPsfNoSM::HscPcaPsfNoSM(CONST_PTR(HscCandidateSet) cs, int nside, int npca, int lanczos_order)
    : HscPcaPsfBase(cs, nside, npca, lanczos_order)
{
    _pca_coeffs = std::vector<double>(_ncand * npca, 0.0);

    _init_pcas();
    _normalize_pcas();
    _optimize_pca_coeffs_and_update();
    _dump("initialized (no SM)", 1);
}


void HscPcaPsfNoSM::optimize()
{
    _optimize_pcas_and_update();
    _dump("optimize_pcas_and_update() called", 1);

    _optimize_xy();
    _optimize_pca_coeffs_and_update();
    _dump("optimize_xy() called", 1);
}


double HscPcaPsfNoSM::eval_shifted_chi2(int icand, double x, double y) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid icand in HscPcaPsfNoSM::eval_shifted_chi2()");

    std::vector<double> basis_images(_npca * _nx * _ny);
    
    for (int ipca = 0; ipca < _npca; ipca++)
        this->_interpolate(&basis_images[ipca*_nx*_ny], &_pca[ipca*_nn], icand, x, y);

    return fit_basis_images(NULL, _npca, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &basis_images[0]);
}


void HscPcaPsfNoSM::_init_pcas()
{
    Eigen::MatrixXd m(_ncand, _nn);

    for (int icand = 0; icand < _ncand; icand++) {
	std::vector<double> img(_nn);
	this->_extirpolate(&img[0], &_im[icand*_nx*_ny], icand);

	double w = 0.0;
	for (int i = 0; i < _nx*_ny; i++)
	    w += _iv[icand*_nx*_ny + i];

	double sqrtw = std::sqrt(w);
	for (int i = 0; i < _nn; i++)
	    m(icand,i) = sqrtw * img[i];
    }

    //
    // FIXME: SVD-decomposing the rectangular matrix M should be
    // more numerically stable than taking eigenvectors of M^T M
    //
    Eigen::MatrixXd a = m.transpose() * m;
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> ev(a);
    
    if (ev.info() != Eigen::Success)
	throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "eigenvector calculation in HscPcaPsfNoSM::_init_pcas() failed?!");

    //
    // A small sanity check on the eigenvector calculation
    //
    Eigen::VectorXd evals = ev.eigenvalues();    
    for (int i = _nn - _npca + 1; i < _nn; i++) {
	if (evals(i-1) >= evals(i))
	    throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "eigenvectors not sorted in HscPcaPsfNoSM::_init_pcas()?!");
    }
    for (int i = 0; i < _nn - _npca; i++) {
	if (evals(i) >= evals(_nn - _npca))
	    throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "eigenvectors not sorted in HscPcaPsfNoSM::_init_pcas()?!");
    }

    Eigen::MatrixXd evecs = ev.eigenvectors();

    for (int ipca = 0; ipca < _npca; ipca++)
	for (int i = 0; i < _nn; i++)
	    this->_pca[ipca*_nn + i] = evecs(i, _nn-ipca-1);
}


// FIXME rethink normalization scheme here...
void HscPcaPsfNoSM::_normalize_pcas()
{
    // Force PCA0 positive
    double t = 0.0;
    for (int i = 0; i < _nn; i++)
	t += _pca[i];
    if (t < 0.0)
	_scale_pca(&_pca[0], -1.0);

    for (int i = 0; i < _npca; i++) {
	for (int j = 0; j < i; j++) {
	    double t = _dot_pca(&_pca[i*_nn], &_pca[j*_nn]);
	    _shift_pca(&_pca[i*_nn], &_pca[j*_nn], -t);
	}
	double t = _dot_pca(&_pca[i*_nn], &_pca[i*_nn]);
	_scale_pca(&_pca[i*_nn], 1.0/sqrt(t));
    }
}


void HscPcaPsfNoSM::_optimize_pcas_and_update()
{
    for (int ipca = 0; ipca < _npca; ipca++) {
	// Invariant: at top of loop, this->ampl and this->psf_arr must be up-to-date

	std::vector<double> grad(_nn, 0.0);
	std::vector<double> t(_nx*_ny);

	for (int icand = 0; icand < _ncand; icand++) {
	    double w = this->_pca_coeffs[icand*_npca + ipca];
	    const double *iv = &this->_iv[icand*_nx*_ny];
	    const double *im = &this->_im[icand*_nx*_ny];
	    const double *psf = &this->_current_psfimg[icand*_nx*_ny];

	    for (int i = 0; i < _nx*_ny; i++)
		t[i] = w * iv[i] * (im[i] - psf[i]);

	    this->_extirpolate(&grad[0], &t[0], icand);
	}

	double a = 0.0;
	double b = 0.0;
	
	for (int icand = 0; icand < _ncand; icand++) {
	    double w = this->_pca_coeffs[icand*_npca + ipca];	    
	    const double *iv = &this->_iv[icand*_nx*_ny];
	    const double *im = &this->_im[icand*_nx*_ny];
	    const double *psf = &this->_current_psfimg[icand*_nx*_ny];

	    this->_interpolate(&t[0], &grad[0], icand);

	    for (int i = 0; i < _nx*_ny; i++) {
		a += w*w * iv[i] * t[i] * t[i];
		b += w * iv[i] * t[i] * (im[i] - psf[i]);
	    }
	}

	_shift_pca(&this->_pca[ipca*_nn], &grad[0], b/a);
	
	// Preserve loop invariant
	this->_normalize_pcas();
	this->_optimize_pca_coeffs_and_update();
    }
}


void HscPcaPsfNoSM::_optimize_pca_coeffs_and_update()
{
    std::vector<double> basis_images(_npca * _nx * _ny);
    std::vector<double> psf_img(_nx * _ny);

    for (int icand = 0; icand < _ncand; icand++) {
	for (int ipca = 0; ipca < _npca; ipca++)
	    this->_interpolate(&basis_images[ipca*_nx*_ny], &_pca[ipca*_nn], icand);

	fit_basis_images(&_pca_coeffs[icand*_npca], _npca, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &basis_images[0]);

        // FIXME could be removed with optional arg to fit_basis_images()
        memset(&psf_img[0], 0, _nx * _ny * sizeof(double));
        for (int ipca = 0; ipca < _npca; ipca++) {
            double w = _pca_coeffs[icand*_npca + ipca];
            for (int i = 0; i < _nx*_ny; i++)
                psf_img[i] += w * basis_images[ipca*_nx*_ny + i];
        }

	this->_updateCurrentPsfImg(icand, &psf_img[0]);        
    }
}


// -------------------------------------------------------------------------------------------------


namespace {
    struct XOptimizer : public GRMinimizer
    {
	const HscPcaPsfNoSM *_p;
	int _icand;

	XOptimizer(const HscPcaPsfNoSM *p, int icand) 
            : GRMinimizer(p->getX(icand)-1.0, p->getX(icand)+1.0, 1.0e-2, 1.0e-2), _p(p), _icand(icand)
        { }

        virtual double eval(double x) const 
        { 
            return _p->eval_shifted_chi2(_icand, x, _p->getY(_icand));
        }
    };

    struct YOptimizer : public GRMinimizer
    {
	const HscPcaPsfNoSM *_p;
	int _icand;

	YOptimizer(const HscPcaPsfNoSM *p, int icand) 
            : GRMinimizer(p->getY(icand)-1.0, p->getY(icand)+1.0, 1.0e-2, 1.0e-2), _p(p), _icand(icand)
        { }

        virtual double eval(double y) const 
        { 
            return _p->eval_shifted_chi2(_icand, _p->getX(_icand), y);
        }
    };
};


void HscPcaPsfNoSM::_optimize_xy()
{
    for (int icand = 0; icand < _ncand; icand++) {
        _current_xy[2*icand] = XOptimizer(this,icand).minimize(_current_xy[2*icand]);
        _current_xy[2*icand+1] = YOptimizer(this,icand).minimize(_current_xy[2*icand+1]);
    }
}


}}}}   // namespace lsst::meas::extensions::hscpsf
