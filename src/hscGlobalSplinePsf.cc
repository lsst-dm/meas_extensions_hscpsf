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
 

#include "lsst/pex/exceptions.h"
#include "lsst/meas/extensions/hscpsf/hscPsf.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif



HscGlobalSplinePsf::HscGlobalSplinePsf(CONST_PTR(HscCandidateSet) cs, int nr, double dr, double sigma)
    : HscSplinePsfBase(cs, nr, dr)
{
    _global_profile = std::vector<double>(_nr, 0.0);
    _candidate_ampl = std::vector<double>(_ncand, 0.0);

    if (sigma <= 0.0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "sigma < 0 in HscGlobalSplinePsf constructor");

    _global_profile[0] = 1.0;
    for (int i = 1; i < _nr; i++) 
	_global_profile[i] = exp(-0.5 * (i/sigma) * (i/sigma));

    this->_optimize_ampl();
    this->_updateGlobalSplineParams();
    this->_dump("initialized", 1);
}


void HscGlobalSplinePsf::optimize()
{
    this->_optimize_global_profile();
    this->_optimize_ampl();
    this->_updateGlobalSplineParams(); 
    this->_dump("optimize_global_profile", 1);

    this->_optimize_xy();
    this->_optimize_gamma();
    this->_optimize_kappa();
    this->_optimize_ampl();
    this->_updateGlobalSplineParams();
    this->_dump("update_xygk", 1);
}


void HscGlobalSplinePsf::_optimize_global_profile()
{
    //
    // It's convenient to implement this by fitting a "super image" with npix = (ncand * nx * ny)
    // using a single call to fit_basis_images().  This uses more memory than it needs to, but
    // for a typical HSC CCD, the (temporary) memory footprint is only ~7 MB, so I haven't bothered
    // to optimize memory usage.
    //
    std::vector<double> basis_images(_nr * _ncand * _nx * _ny);

    for (int ir = 0; ir < _nr; ir++) {
	for (int icand = 0; icand < _ncand; icand++) {
            // FIXME is this OK?
            assert(_candidate_ampl[icand] > 0.0);

            std::vector<double> profile(_nr, 0.0);
            profile[ir] = _candidate_ampl[icand];

            double *bim = &basis_images[(ir*_ncand+icand) * _nx * _ny];
            this->_fillSplineImage(bim, &profile[0], icand, getX(icand), getY(icand), getGamma1(icand), getGamma2(icand), getKappa(icand));
	}
    }

    fit_basis_images(&_global_profile[0], _nr, _ncand*_nx*_ny, _iv, _im, &basis_images[0]);
}


void HscGlobalSplinePsf::_optimize_ampl()
{
    std::vector<double> psf_image(_nx * _ny, 0.0);

    for (int icand = 0; icand < _ncand; icand++) {
        this->_fillSplineImage(&psf_image[0], &_global_profile[0], icand, getX(icand), getY(icand), getGamma1(icand), getGamma2(icand), getKappa(icand));
        fit_basis_images(&_candidate_ampl[icand], 1, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &psf_image[0]);
        assert(_candidate_ampl[icand] > 0.0);   // FIXME
    }
}


void HscGlobalSplinePsf::_updateGlobalSplineParams()
{
    // FIXME rethink this!
    double normalization = _global_profile[0];
    assert(normalization > 0.0);

    // normalize profiles
    for (int ir = 0; ir < _nr; ir++)
        _global_profile[ir] /= normalization;
    for (int icand = 0; icand < _ncand; icand++)
        _candidate_ampl[icand] *= normalization;

    for (int icand = 0; icand < _ncand; icand++)
        for (int ir = 0; ir < _nr; ir++)
            _profile[icand*_nr+ir] = _candidate_ampl[icand] * _global_profile[ir];

    this->_updateSplineParams();
}


void HscGlobalSplinePsf::_dump(const char *msg, int level) const
{
    if (level < 1)
        return;

    std::cerr << "HscGlobalSplinePsf: " << msg << ": ncand=" << _ncand << ", chi2=" << get_total_reduced_chi2() << std::endl;

    if (level < 2)
        return;

    std::cerr << "   global profile =";
    for (int i = 0; i < _nr; i++)
        std::cerr << " " << _global_profile[i];
    std::cerr << std::endl;

    for (int icand = 0; icand < _ncand; icand++) {
        std::cerr << "    candidate " << icand << ": " << _current_xy[2*icand] << " " << _current_xy[2*icand+1]
                  << " " << _gamma[2*icand] << " " << _gamma[2*icand+1] << " " << _kappa[icand]
                  << " " << _candidate_ampl[icand] 
                  << " " << (_residual_chi2[icand]/_ndof[icand]) << std::endl;
    }
}



}}}}   // namespace lsst::meas::extensions::hscpsf
