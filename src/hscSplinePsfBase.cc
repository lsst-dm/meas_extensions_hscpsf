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
#include "lsst/meas/extensions/hscpsf/GRMinimizer.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


// Helper functions for spline initialization/evaluation
namespace {

    // "spline boundary condition": used when initializing a spline, to define boundary conditions at left and right endpoints
    struct spline_bc {
        //
        // if is_natural=true, then the boundary condition is d^2y/dx^2 = 0 and the value of 'deriv' is ignored
        // if is_natural=false, then the boundary condition is dy/dx = deriv
        //
        bool is_natural;
        double deriv;
        
        // default constructor defines "natural" boundary conditions (d^2y/dx^2 = 0)
        spline_bc() : is_natural(true), deriv(0.0) { }
        
        // this constructor defines boundary conditions corresponding to fixed dy/dx
        explicit spline_bc(double deriv_) : is_natural(false), deriv(deriv_) { }
        
        spline_bc(bool is_natural_, double deriv_) : is_natural(is_natural_), deriv(deriv_) { }
    };


    void spline_init(int n, const double *x_in, const double *y_in, const spline_bc &lderiv, const spline_bc &rderiv, double *y2_out)
    {
        assert(n >= 2);
        assert(x_in != NULL);
        assert(y_in != NULL);
        assert(y2_out != NULL);
        
        std::vector<double> u(n);
        
        if (lderiv.is_natural) {
            u[0] = y2_out[0] = 0.0;
        }
        else {
            double dx = x_in[1] - x_in[0];
            double dy = y_in[1] - y_in[0];
            u[0] = 3./dx * (dy/dx - lderiv.deriv);
            y2_out[0] = -0.5;
        }
        
        for (int i = 1; i <= n-2; i++) {
            double sig = (x_in[i]-x_in[i-1]) / (x_in[i+1]-x_in[i-1]);
            double p = sig * y2_out[i-1] + 2.0;
            
            y2_out[i] = (sig-1.0)/p;
            u[i] = (y_in[i+1]-y_in[i])/(x_in[i+1]-x_in[i]) - (y_in[i]-y_in[i-1])/(x_in[i]-x_in[i-1]);
            u[i] = (6.0*u[i]/(x_in[i+1]-x_in[i-1]) - sig*u[i-1]) / p;
        }
        
        if (rderiv.is_natural) {
            y2_out[n-1] = 0.0;
        }
        else {
            double dx = x_in[n-1] - x_in[n-2];
            double dy = y_in[n-1] - y_in[n-2];
            double qn = 0.5; 
            double un = 3.0/dx * (rderiv.deriv - dy/dx);
            y2_out[n-1] = (un - qn*u[n-2]) / (qn*y2_out[n-2] + 1.0);
        }
        
        for (int i = n-2; i >= 0; i--)
            y2_out[i] = y2_out[i]*y2_out[i+1] + u[i];
    }


    inline double splinterp(double x, double x1, double x2, double f1, double f2, double ddf1, double ddf2)
    {
        double  h   = (x2-x1);
        double  dx  = (x-x1)/h;
        double  dx3 = dx*dx*dx - dx;
        double  cx3 = (1-dx)*(1-dx)*(1-dx) - (1-dx);
        
        return (1-dx)*f1 + (dx)*f2 + h*h/6*(cx3*ddf1 + dx3*ddf2);
    }
};


inline int computeNside(int nr, double dr)
{
    if ((nr <= 0) || (dr <= 0.0))
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "invalid parametrs in HscSplinePsfBase constructor");
    return int(nr*dr + 1.0e-10) + 1;
}


HscSplinePsfBase::HscSplinePsfBase(CONST_PTR(HscCandidateSet) cs, int nr, double dr)
    : HscPsfBase(cs, computeNside(nr,dr)), _nr(nr), _dr(dr)
{
    _gamma = std::vector<double>(_ncand * 2, 0.0);
    _kappa = std::vector<double>(_ncand, 0.0);
    _profile = std::vector<double>(_ncand * _nr, 0.0);
    _basis_profiles_dd = std::vector<double>(_nr * (_nr+1), 0.0);

    for (int i = 0; i < _nr; i++) {
        std::vector<double> x(_nr+1);
        std::vector<double> y(_nr+1);

        for (int j = 0; j < _nr+1; j++) {
            x[j] = j;
            y[j] = (i==j) ? 1.0 : 0.0;
        }

        spline_init(_nr+1, &x[0], &y[0], spline_bc(0.0), spline_bc(0.0), &_basis_profiles_dd[i*(_nr+1)]);
    }
}


double HscSplinePsfBase::getGamma1(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscSplinePsfBase::getGamma1()");
    return _gamma[2*icand];
}

double HscSplinePsfBase::getGamma2(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscSplinePsfBase::getGamma2()");
    return _gamma[2*icand + 1];
}

double HscSplinePsfBase::getKappa(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscSplinePsfBase::getKappa()");
    return _kappa[icand];
}

const double *HscSplinePsfBase::getProfile(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscSplinePsfBase::getProfile()");
    return &_profile[icand * _nr];
}


// -------------------------------------------------------------------------------------------------


namespace {
    struct OptimizerBase : public GRMinimizer
    {
        const HscSplinePsfBase *_p;
        int _icand;
        mutable std::vector<double> _scratch;

        OptimizerBase(const HscSplinePsfBase *p, int icand, double xmin, double xmax, double xres, double xstep)
            : GRMinimizer(xmin,xmax,xres,xstep), _p(p), _icand(icand), _scratch(p->getNx() * p->getNy())
        { }

        double eval_chi2(double x, double y, double gamma1, double gamma2, double kappa) const
        {
            _p->_fillSplineImage(&_scratch[0], _p->getProfile(_icand), _icand, x, y, gamma1, gamma2, kappa);
            return HscPsfBase::fit_basis_images(NULL, 1, _p->getNx() * _p->getNy(), _p->getIV(_icand), _p->getIm(_icand), &_scratch[0]);
        }
    };

    struct XOptimizer : public OptimizerBase
    {
        XOptimizer(const HscSplinePsfBase *p, int icand) 
            : OptimizerBase(p, icand, p->getXini(icand)-2.0, p->getXini(icand)+2.0, 1.0e-2, 1.0e-2)
        { }

        virtual double eval(double x) const 
        { 
            return this->eval_chi2(x, _p->getY(_icand), _p->getGamma1(_icand), _p->getGamma2(_icand), _p->getKappa(_icand));
        }
    };

    struct YOptimizer : public OptimizerBase
    {
        YOptimizer(const HscSplinePsfBase *p, int icand) 
            : OptimizerBase(p, icand, p->getYini(icand)-2.0, p->getYini(icand)+2.0, 1.0e-2, 1.0e-2)
        { }

        virtual double eval(double y) const 
        { 
            return this->eval_chi2(_p->getX(_icand), y, _p->getGamma1(_icand), _p->getGamma2(_icand), _p->getKappa(_icand));
        }
    };

    struct Gamma1Optimizer : public OptimizerBase
    {
        Gamma1Optimizer(const HscSplinePsfBase *p, int icand) 
            : OptimizerBase(p, icand, -0.2, 0.2, 1.0e-3, 1.0e-3)
        { }

        virtual double eval(double gamma1) const 
        { 
            return this->eval_chi2(_p->getX(_icand), _p->getY(_icand), gamma1, _p->getGamma2(_icand), _p->getKappa(_icand));
        }
    };

    struct Gamma2Optimizer : public OptimizerBase
    {
        Gamma2Optimizer(const HscSplinePsfBase *p, int icand) 
            : OptimizerBase(p, icand, -0.2, 0.2, 1.0e-3, 1.0e-3)
        { }

        virtual double eval(double gamma2) const 
        { 
            return this->eval_chi2(_p->getX(_icand), _p->getY(_icand), _p->getGamma1(_icand), gamma2, _p->getKappa(_icand));
        }
    };

    struct KappaOptimizer : public OptimizerBase
    {
        KappaOptimizer(const HscSplinePsfBase *p, int icand) 
            : OptimizerBase(p, icand, -0.2, 0.2, 1.0e-3, 1.0e-3)
        { }

        virtual double eval(double kappa) const 
        { 
            return this->eval_chi2(_p->getX(_icand), _p->getY(_icand), _p->getGamma1(_icand), _p->getGamma2(_icand), kappa);
        }
    };
};


void HscSplinePsfBase::_optimize_xy()
{
    // FIXME max pixel excursion is currently a hardcoded parameter

    for (int icand = 0; icand < _ncand; icand++) {
        _current_xy[2*icand] = XOptimizer(this,icand).minimize(_current_xy[2*icand], true);
        _current_xy[2*icand+1] = YOptimizer(this,icand).minimize(_current_xy[2*icand+1], true);
    }
}


void HscSplinePsfBase::_optimize_gamma()
{
    for (int icand = 0; icand < _ncand; icand++) {
        _gamma[2*icand] = Gamma1Optimizer(this,icand).minimize(_gamma[2*icand], true);
        _gamma[2*icand+1] = Gamma2Optimizer(this,icand).minimize(_gamma[2*icand+1], true);
    }
}


void HscSplinePsfBase::_optimize_kappa()
{
    for (int icand = 0; icand < _ncand; icand++)
        _kappa[icand] = KappaOptimizer(this,icand).minimize(_kappa[icand], true);
}


bool HscSplinePsfBase::displacement_is_maximized(int icand) const
{
    static const double max_dx = 2.0;    // FIXME hardcoded

    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "bad icand in HscSplinePsfBase::displacement_is_maximized()");

    bool dx_maximized = (std::fabs(_current_xy[2*icand] - _initial_xy[2*icand]) >= 0.999*max_dx);
    bool dy_maximized = (std::fabs(_current_xy[2*icand+1] - _initial_xy[2*icand+1]) >= 0.999*max_dx);
    return dx_maximized || dy_maximized;
}

bool HscSplinePsfBase::gamma_is_maximized(int icand) const
{
    static const double max_gamma = 0.2;   // FIXME hardcoded

    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "bad icand in HscSplinePsfBase::gamma_is_maximized()");

    bool gamma1_maximized = (std::fabs(_gamma[2*icand]) >= 0.999*max_gamma);
    bool gamma2_maximized = (std::fabs(_gamma[2*icand+1]) >= 0.999*max_gamma);
    return gamma1_maximized || gamma2_maximized;
}

bool HscSplinePsfBase::kappa_is_maximized(int icand) const
{
    static const double max_kappa = 0.2;

    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "bad icand in HscSplinePsfBase::kappa_is_maximized()");

    return (std::fabs(_kappa[icand]) >= 0.999*max_kappa);
}


void HscSplinePsfBase::_updateSplineParams()
{
    std::vector<double> psf_img(_nx * _ny);

    for (int icand = 0; icand < _ncand; icand++) {
        this->_fillSplineImage(&psf_img[0], &_profile[icand*_nr], icand, getX(icand), getY(icand), getGamma1(icand), getGamma2(icand), getKappa(icand));
        this->_updateCurrentPsfImg(icand, &psf_img[0]);
    }
}


// -------------------------------------------------------------------------------------------------


void HscSplinePsfBase::_fillSplineImage(double *out, const double *profile, int icand, double x, double y, double gamma1, double gamma2, double kappa) const
{
    std::vector<double> in(_nr+1, 0.0);
    std::vector<double> in_dd(_nr+1, 0.0);

    for (int i = 0; i < _nr; i++) {
        in[i] = profile[i];
        for (int j = 0; j < _nr+1; j++)
            in_dd[j] += profile[i] * _basis_profiles_dd[i*(_nr+1)+j];
    }

    double axx, axy, ayy;
    make_shear_matrix(axx, axy, ayy, gamma1, gamma2, kappa);

    double x0 = _xy0[2*icand];
    double y0 = _xy0[2*icand+1];

    for (int i = 0; i < _nx; i++) {
	for (int j = 0; j < _ny; j++) {
	    // pixel displacement relative to candidate
	    double dx0 = (i+x0) - x;
	    double dy0 = (j+y0) - y;

	    // sheared displacement
	    double dx = axx*dx0 + axy*dy0;
	    double dy = axy*dx0 + ayy*dy0;

	    // distance in "spline units"
	    double t = sqrt(dx*dx + dy*dy) / this->_dr;

	    int it = (int)t;
	    out[i*_ny+j] = (it < _nr) ? splinterp(t, it, it+1, in[it], in[it+1], in_dd[it], in_dd[it+1]) : 0.0;
	}
    }
}


}}}}   // namespace lsst::meas::extensions::hscpsf
