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


HscClipPsf::HscClipPsf(CONST_PTR(HscCandidateSet) cs, int nr, double dr)
    : HscSplinePsfBase(cs, nr, dr)
{
    this->_optimize_profile();
    this->_updateSplineParams();
    this->_dump("initialized", 1);
}


void HscClipPsf::optimize()
{
    this->_optimize_xy();
    this->_optimize_gamma();
    this->_optimize_profile();
    this->_updateSplineParams();
    this->_dump("optimize() called", 1);
}


void HscClipPsf::_optimize_profile()
{
    std::vector<double> basis_images(_nr * _nx * _ny);

    for (int icand = 0; icand < _ncand; icand++) {
	for (int i = 0; i < _nr; i++) {
            std::vector<double> t(_nr, 0.0);
            t[i] = 1.0;
	    this->_fillSplineImage(&basis_images[i*_nx*_ny], &t[0], icand, getX(icand), getY(icand), getGamma1(icand), getGamma2(icand), 0.0);
        }

        fit_basis_images(&_profile[icand*_nr], _nr, _nx*_ny, &_iv[icand*_nx*_ny], &_im[icand*_nx*_ny], &basis_images[0]);
    }
}


void HscClipPsf::_dump(const char *msg, int level) const
{
    if (level < 1)
        return;

    std::cerr << "HscClipPsf: " << msg << ": ncand=" << _ncand << ", chi2=" << get_total_reduced_chi2() << std::endl;

    if (level < 2)
        return;

    for (int icand = 0; icand < _ncand; icand++) {
        std::cerr << "    " << icand << ": " << _current_xy[2*icand] << " " << _current_xy[2*icand+1]
                  << " " << _gamma[2*icand] << " " << _gamma[2*icand+1]
                  << " " << (_residual_chi2[icand]/_ndof[icand]) << std::endl; 
    }
}


}}}}   // namespace lsst::meas::extensions::hscpsf
