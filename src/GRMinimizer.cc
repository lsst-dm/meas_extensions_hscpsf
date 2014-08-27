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
 

#include "lsst/meas/extensions/hscpsf/GRMinimizer.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


GRMinimizer::GRMinimizer(double xmin, double xmax, double xres, double xstep)
    : _xmin(xmin), _xmax(xmax), _xres(xres), _xstep(xstep)
{
    if ((xmin >= xmax) || (xres <= 0.0) || (xstep <= 0.0))
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "invalid call to GRMinimizer constructor");
}

double GRMinimizer::minimize(double x, bool allow_endpoints)
{
    if ((x < _xmin) || (x > _xmax))
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "GRMinimizer::minimze() argument is out of range");

    int res = _prebracket(x);

    if ((res != 0) && !allow_endpoints)
        throw GRMinimizer::OutOfRangeException();

    if (res > 0)
        return (_x1+_x2)/2.;
    if (res < 0)
        return (_x0+_x1)/2.;

#if 0
    // FIXME think about this more caerfully
    if (_y0==_y1 || _y1==_y2)
        std::cerr << "GRMinimizer: warning: tie in _prebracket(), suspect xstep is too small\n";
#endif

    _bisect_to_convergence();
    return (_x0+_x2)/2.;
}


}}}}   // namespace lsst::meas::extensions::hscpsf
