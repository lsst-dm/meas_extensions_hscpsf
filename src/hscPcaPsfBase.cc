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

namespace afwImage = lsst::afw::image;


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


HscPcaPsfBase::HscPcaPsfBase(CONST_PTR(HscCandidateSet) cs, int nside, int npca, int lanczos_order)
    : HscPsfBase(cs,nside), _npca(npca), _lanczos_order(lanczos_order), _nn((2*nside+1)*(2*nside+1))
{
    if ((npca <= 1) || (npca > _ncand) || (npca > (2*nside+1)*(2*nside+1)) || (nside < 2) || (lanczos_order < 3))
	throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters to HscPcaPsfBase constructor");

    _pca = std::vector<double>(npca * (2*nside+1) * (2*nside+1), 0.0);
}


void HscPcaPsfBase::_interpolate(double *out, const double *in, int icand, double x, double y) const
{
    assert(out != NULL);
    assert(in != NULL);
    assert(icand >= 0 && icand < _ncand);

    int n = 2*_nside + 1;

    lanczos_offset_2d(this->_lanczos_order, 
		      this->_xy0[2*icand] - x + _nside,        // x0
		      this->_xy0[2*icand+1] - y + _nside,      // y0
		      this->_nx, this->_ny, out, this->_ny,
		      n, n, in, n, 
		      true);
}


void HscPcaPsfBase::_interpolate(double *out, const double *in, int icand) const
{
    assert(icand >= 0 && icand < _ncand);
    this->_interpolate(out, in, icand, this->_current_xy[2*icand], this->_current_xy[2*icand+1]);
}


void HscPcaPsfBase::_extirpolate(double *out, const double *in, int icand, double x, double y) const
{
    assert(in != NULL);
    assert(out != NULL);
    assert(icand >= 0 && icand < _ncand);

    int n = 2*_nside + 1;

    lanczos_offset_2d(this->_lanczos_order, 
		      -this->_xy0[2*icand] + x - _nside,
		      -this->_xy0[2*icand+1] + y - _nside,
		      n, n, out, n,
		      this->_nx, this->_ny, in, this->_ny,
		      true, true);     // last "true" is "accumulate"
}


void HscPcaPsfBase::_extirpolate(double *out, const double *in, int icand) const
{
    assert(icand >= 0 && icand < _ncand);
    this->_extirpolate(out, in, icand, this->_current_xy[2*icand], this->_current_xy[2*icand+1]);
}


void HscPcaPsfBase::_dump(const char *msg, int level) const
{
    if (level < 1)
	return;

    std::cerr << "HscPcaPsfBase: " << msg << ": ncand=" << _ncand << ", chi2=" << get_total_reduced_chi2() << std::endl;

    if (level < 2)
        return;

    for (int icand = 0; icand < _ncand; icand++) {
        std::cerr << "    candidate " << icand << ": " << _current_xy[2*icand] << " " << _current_xy[2*icand+1]
                  << " " << (_residual_chi2[icand]/_ndof[icand]) << std::endl;
    }
}


void HscPcaPsfBase::_scale_pca(double *out, double t) const
{
    assert(out != NULL);

    for (int i = 0; i < _nn; i++)
	out[i] *= t;
}

void HscPcaPsfBase::_shift_pca(double *out, const double *in, double t) const
{
    assert(out != NULL && in != NULL);
    for (int i = 0; i < _nn; i++)
	out[i] += t * in[i];
}

double HscPcaPsfBase::_dot_pca(const double *in1, const double *in2) const
{
    assert(in1 != NULL && in2 != NULL);
    
    double ret = 0.0;
    for (int i = 0; i < _nn; i++)
	ret += in1[i] * in2[i];

    return ret;
}



}}}}   // namespace lsst::meas::extensions::hscpsf
