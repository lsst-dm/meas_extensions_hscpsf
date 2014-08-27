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
 

#include "lsst/meas/extensions/hscpsf/hscPsf.h"
#include "lsst/meas/extensions/hscpsf/GRMinimizer.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


static inline void copyVector(std::vector<double> &dst, const std::vector<double> &src)
{
    assert(dst.size() == src.size());
    memcpy(&dst[0], &src[0], dst.size() * sizeof(double));
}


HscPcaPsf::HscPcaPsf(CONST_PTR(HscPcaPsfBase) psf0, CONST_PTR(HscSpatialModelBase) spatialModel)
    : HscPcaPsfBase(psf0->getCandSet(), psf0->getNside(), psf0->getNpca(), psf0->getLanczosOrder())
{
    _sm_ncoeffs = spatialModel->getNcoeffs();
    _sm_coeffs = std::vector<double>((_npca-1) * _sm_ncoeffs, 0.0);
    _candidate_ampl = std::vector<double>(_ncand, 0.0);
    _spatial_model = spatialModel;

    for (int icand = 0; icand < _ncand; icand++) {
        this->_current_xy[2*icand] = psf0->getX(icand);
        this->_current_xy[2*icand+1] = psf0->getY(icand);
    }

    copyVector(this->_pca, psf0->getPcas());

    this->_optimize_ampl_and_update();
    this->_dump("initialized", 1);
}


void HscPcaPsf::optimize()
{
    this->_optimize_spatial_model_and_update();
    this->_dump("optimize_spatial_model() called", 1);

    this->_optimize_pcas_and_update();
    this->_dump("optimize_pcas() called", 1);

    this->_optimize_xy();
    this->_optimize_ampl_and_update();
    this->_dump("optimize_xy() called", 1);
}


//
// FIXME nearly the same as HscPcaPsf::_optimize_pcas()
// How to combine both routines into a common routine in HscPcaPsfBase?
//
void HscPcaPsf::_optimize_pcas_and_update()
{
    // FIXME for now, we evaluate the spatial model at the _initial_ position!
    std::vector<double> sm_eval(_ncand * (_npca-1));
    _spatial_model->eval(&sm_eval[0], _ncand, &_initial_xy[0], _npca-1, &_sm_coeffs[0]);

    for (int ipca = 0; ipca < _npca; ipca++) {
	std::vector<double> grad(_nn, 0.0);
	std::vector<double> t(_nx*_ny);

	for (int icand = 0; icand < _ncand; icand++) {
            double w = _candidate_ampl[icand] * ((ipca > 0) ? sm_eval[icand*(_npca-1)+(ipca-1)] : 1.0);

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
            double w = _candidate_ampl[icand] * ((ipca > 0) ? sm_eval[icand*(_npca-1)+(ipca-1)] : 1.0);

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

	// Preserve invariants
        _spatial_model->normalizePcaImages(_npca, _nn, _ncand, &_pca[0], &_candidate_ampl[0], &_sm_coeffs[0]);	
	this->_optimize_ampl_and_update();
    }
}


//
// Warning: must be up-to-date before calling!
//
void HscPcaPsf::_optimize_spatial_model_and_update()
{
    for (int ism = 0; ism < _npca-1; ism++) {
        std::vector<double> dpsf(_nx*_ny);
        std::vector<double> a(_ncand, 0.0);
        std::vector<double> b(_ncand, 0.0);

	for (int icand = 0; icand < _ncand; icand++) {
	    const double *iv = &_iv[icand*_nx*_ny];
	    const double *im = &_im[icand*_nx*_ny];
	    const double *psf = &_current_psfimg[icand*_nx*_ny];
            double ampl = _candidate_ampl[icand];

	    this->_interpolate(&dpsf[0], &_pca[(ism+1)*_nn], icand);

	    for (int i = 0; i < _nx*_ny; i++) {
		a[icand] += ampl * ampl * iv[i] * dpsf[i] * dpsf[i];
		b[icand] += ampl * iv[i] * dpsf[i] * (im[i] - psf[i]);
	    }
	}

        std::vector<double> dcoeffs(_sm_ncoeffs);
	_spatial_model->optimize(&dcoeffs[0], _ncand, &_initial_xy[0], &a[0], &b[0]);

        for (int i = 0; i < _sm_ncoeffs; i++)
            _sm_coeffs[ism*_sm_ncoeffs + i] += dcoeffs[i];

	this->_optimize_ampl_and_update();
    }
}


void HscPcaPsf::_optimize_ampl_and_update()
{
    // FIXME for now, we evaluate the spatial model at the _initial_ position!
    std::vector<double> sm_eval(_ncand * (_npca-1));
    _spatial_model->eval(&sm_eval[0], _ncand, &_initial_xy[0], _npca-1, &_sm_coeffs[0]);

    std::vector<double> psf_img(_nx*_ny);
    std::vector<double> psf_preimg(_nn);

    for (int icand = 0; icand < _ncand; icand++) {
	memcpy(&psf_preimg[0], &_pca[0], _nn * sizeof(double));
	for (int ism = 0; ism < _npca-1; ism++)
            this->_shift_pca(&psf_preimg[0], &_pca[(ism+1)*_nn], sm_eval[icand*(_npca-1)+ism]);

	this->_interpolate(&psf_img[0], &psf_preimg[0], icand);

        double ampl = 0.0;
	fit_basis_images(&ampl, 1, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &psf_img[0]);

        _candidate_ampl[icand] = ampl;

        for (int i = 0; i < _nx*_ny; i++)
            psf_img[i] *= ampl;

        this->_updateCurrentPsfImg(icand, &psf_img[0]);
    }
}


// virtual
// FIXME remove redundancy between this routine and others
void HscPcaPsf::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    if (out == NULL)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "out == NULL in HscPcaPsf::eval()");
    if (nx_out <= 0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "nx_out <= 0 in HscPcaPsf::eval()");
    if (ny_out <= 0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "ny_out <= 0 in HscPcaPsf::eval()");

    std::vector<double> xy(2);
    xy[0] = x;
    xy[1] = y;

    std::vector<double> sm_eval(_npca-1);
    _spatial_model->eval(&sm_eval[0], 1, &xy[0], _npca-1, &_sm_coeffs[0]);
    
    std::vector<double> psf_preimg(_nn);
    memcpy(&psf_preimg[0], &_pca[0], _nn * sizeof(double));
    for (int ism = 0; ism < _npca-1; ism++)
        this->_shift_pca(&psf_preimg[0], &_pca[(ism+1)*_nn], sm_eval[ism]);

    int n = 2*_nside + 1;
    
    lanczos_offset_2d(_lanczos_order,
		      x0 - x + _nside,       // x0
		      y0 - y + _nside,       // y0
		      nx_out, ny_out, out, ny_out,
		      n, n, &psf_preimg[0], n,
		      true);                // zero_pad
}


// -------------------------------------------------------------------------------------------------
//
// FIXME mostly cut-and-paste from HscPcaPsfNoSM


// FIXME cut-and-paste with _optimize_ampl_and_update()
double HscPcaPsf::eval_shifted_chi2(int icand, double x, double y) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid icand in HscPcaPsf::eval_shifted_chi2()");
    
    // FIXME for now, we evaluate the spatial model at the _initial_ position!
    std::vector<double> sm_eval(_npca-1);
    _spatial_model->eval(&sm_eval[0], 1, &_initial_xy[2*icand], _npca-1, &_sm_coeffs[0]);

    std::vector<double> psf_img(_nx*_ny);
    std::vector<double> psf_preimg(_nn);

    memcpy(&psf_preimg[0], &_pca[0], _nn * sizeof(double));
    for (int ism = 0; ism < _npca-1; ism++)
	this->_shift_pca(&psf_preimg[0], &_pca[(ism+1)*_nn], sm_eval[ism]);

    this->_interpolate(&psf_img[0], &psf_preimg[0], icand, x, y);
    return fit_basis_images(NULL, 1, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &psf_img[0]);
}


namespace {
    struct XOptimizer : public GRMinimizer
    {
	const HscPcaPsf *_p;
	int _icand;

	XOptimizer(const HscPcaPsf *p, int icand) 
            : GRMinimizer(p->getX(icand)-1.0, p->getX(icand)+1.0, 1.0e-2, 1.0e-2), _p(p), _icand(icand)
        { }

        virtual double eval(double x) const 
        { 
            return _p->eval_shifted_chi2(_icand, x, _p->getY(_icand));
        }
    };

    struct YOptimizer : public GRMinimizer
    {
	const HscPcaPsf *_p;
	int _icand;

	YOptimizer(const HscPcaPsf *p, int icand) 
            : GRMinimizer(p->getY(icand)-1.0, p->getY(icand)+1.0, 1.0e-2, 1.0e-2), _p(p), _icand(icand)
        { }

        virtual double eval(double y) const 
        { 
            return _p->eval_shifted_chi2(_icand, _p->getX(_icand), y);
        }
    };
};


void HscPcaPsf::_optimize_xy()
{
    for (int icand = 0; icand < _ncand; icand++) {
        _current_xy[2*icand] = XOptimizer(this,icand).minimize(_current_xy[2*icand]);
        _current_xy[2*icand+1] = YOptimizer(this,icand).minimize(_current_xy[2*icand+1]);
    }
}



}}}}   // namespace lsst::meas::extensions::hscpsf
