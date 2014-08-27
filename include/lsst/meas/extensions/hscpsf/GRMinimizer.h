// -*- lsst-c++ -*-
/*
 * LSST Data Management System
 * Copyright 2008-2013 LSST Corporation.
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

#ifndef LSST_MEAS_EXTENSIONS_HSCPSF_GRMINIMIZER_H
#define LSST_MEAS_EXTENSIONS_HSCPSF_GRMINIMIZER_H 1

#include "lsst/pex/exceptions.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


//
// Generic 1D function minimizer, using golden ratio method
// The function being minimized is defined by subclassing
//
class GRMinimizer {
public:
    GRMinimizer(double xmin, double xmax, double xres, double xstep);

    double minimize(double x, bool allow_endpoints=false);

    // function to minimize, specified by subclassing
    virtual double eval(double x) const = 0;

    //
    // This exception will be thrown if minimize() is called with allow_endpoints=false,
    // and an endpoint of the interval [xmin,xmax] is reached.
    //
    // FIXME: make subclass of lsst::pex::exceptions::Exception?
    //
    class OutOfRangeException : public std::exception { };

protected:
    double _xmin, _xmax, _xres, _xstep;

private:
    double _x0, _x1, _x2;
    double _y0, _y1, _y2;

    // one iteration of golden ratio minimization
    inline void _bisect_once()
    {
	if (_x1-_x0 > _x2-_x1) {
	    // trial point on low side
	    double x = _x1 - 0.38197*(_x1-_x0);
	    double y = eval(x);

	    if (y < _y1) {
		_x2 = _x1;
		_y2 = _y1;
		_x1 = x;
		_y1 = y;
	    }
	    else {
		_x0 = x;
		_y0 = y;
	    }
	}
	else {
	    // high side
	    double x = _x1 + 0.38197*(_x2-_x1);
	    double y = eval(x);

	    if (y < _y1) {
		_x0 = _x1;
		_y0 = _y1;
		_x1 = x;
		_y1 = y;
	    }
	    else {
		_x2 = x;
		_y2 = y;
	    }
	}

	assert(_x0 < _x1 && _x1 < _x2);
	assert(_y1 <= _y0 && _y1 <= _y2);
    }


    // iterate to convergence
    inline void _bisect_to_convergence()
    {
	assert(_x0 < _x1 && _x1 < _x2);
	assert(_y1 <= _y0 && _y1 <= _y2);

	while (_x2-_x0 > _xres)
	    _bisect_once();
    }


    // returns 0 if bracketing successful, -1 if _xmin was reached
    inline int _search_down()
    {
	for (;;) {
	    assert(_xmin <= _x0 && _x0 < _x1);
	    assert(_y0 < _y1);
	    
	    // special case where _xmin is reached
	    if (_x0==_xmin) {
		_x2 = _x1;
		_y2 = _y1;
		_x1 = std::min(_xmin+_xres, (_x0+_x2)/2.);
		_y1 = eval(_x1);
		return (_y0 < _y1) ? -1 : 0;
	    }

	    _x2 = _x1;
	    _y2 = _y1;
	    _x1 = _x0;
	    _y1 = _y0;

	    _x0 = std::max(_xmin, _x1 - 1.6*(_x2-_x1));
	    _y0 = eval(_x0);

	    if (_y0 >= _y1)
		return 0;
	}
    }


    // returns 0 if bracketing successful, 1 if _xmax was reached
    inline int _search_up()
    {
	for (;;) {
	    assert(_x1 < _x2 && _x2 <= _xmax);
	    assert(_y2 < _y1);
	    
	    // special case where _xmax is reached
	    if (_x2==_xmax) {
		_x0 = _x1;
		_y0 = _y1;
		_x1 = std::max(_xmax-_xres, (_x0+_x2)/2.);
		_y1 = eval(_x1);
		return (_y2 < _y1) ? 1 : 0;
	    }

	    _x0 = _x1;
	    _y0 = _y1;
	    _x1 = _x2;
	    _y1 = _y2;

	    _x2 = std::min(_xmax, _x1 + 1.6*(_x1-_x0));
	    _y2 = eval(_x2);

	    if (_y2 >= _y1)
		return 0;
	}
    }


    // returns 0 if bracketing successful, -1 if _xmin was reached, 1 if _xmax was reached
    inline int _prebracket(double xini)
    {
	assert(_xmin <= xini && xini <= _xmax);
	
	// adjust xini in corner cases
	double xmid = (_xmin+_xmax)/2.;
	xini = std::max(xini, std::min(xmid,_xmin+0.99*_xstep));
	xini = std::min(xini, std::max(xmid,_xmax-0.99*_xstep));
	
	_x0 = std::max(_xmin, xini-_xstep);
	_x1 = xini;
	_x2 = std::min(_xmax, xini+_xstep);
	
	_y0 = eval(_x0);
	_y1 = eval(_x1);
	_y2 = eval(_x2);

	if (_y0 < _y2)
	    return (_y1 <= _y0) ? 0 : _search_down();
	else
	    return (_y1 <= _y2) ? 0 : _search_up();
    }
};


}}}}   // namespace lsst::meas::extensions::hscpsf


#endif  // LSST_MEAS_EXTENSIONS_HSCPSF_GRMINIMIZER_H
