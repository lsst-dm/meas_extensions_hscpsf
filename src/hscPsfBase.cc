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
 
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/extensions/hscpsf/hscPsf.h"

namespace afwImage = lsst::afw::image;


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


HscPsfBase::HscPsfBase(CONST_PTR(HscCandidateSet) cs, int nside)
    : _cs(cs), _ncand(cs->getNcand()), 
      _nx(cs->getNx()), _ny(cs->getNy()), _nside(nside),
      _im(cs->getIm()), _iv(cs->getIV()), 
      _xy0(cs->getXY0()), _initial_xy(cs->getXY()), 
      _detection_chi2(cs->getChi2()), _ndof(cs->getNdof()),
      _current_xy(2 * cs->getNcand()),
      _current_psfimg(cs->getNcand() * cs->getNx() * cs->getNy(), 0),
      _residual_chi2(cs->getNcand()),
      _bad(cs->getNcand(), false)
{ 
    if (_ncand < 1)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "empty HscCandidateSet in HscPsfBase constructor");
    if (nside < 1)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "nside <= 1 in HscPsfBase constructor");

    for (int i = 0; i < 2*_ncand; i++)
        _current_xy[i] = _initial_xy[i];
    for (int i = 0; i < _ncand; i++)
        _residual_chi2[i] = _detection_chi2[i];

    _xmean = _ymean = 0.0;
    _xmin = _xmax = _initial_xy[0];
    _ymin = _ymax = _initial_xy[1];

    for (int i = 0; i < _ncand; i++) {
        double x = _initial_xy[2*i];
        double y = _initial_xy[2*i+1];

        _xmin = std::min(_xmin, x);
        _xmax = std::max(_xmax, x);
        _ymin = std::min(_ymin, y);
        _ymax = std::max(_ymax, y);
        _xmean += x / _ncand;
        _ymean += y / _ncand;
    }
}


const double *HscPsfBase::getIm(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getIm()");
    return &_im[icand * _nx * _ny];
}

const double *HscPsfBase::getIV(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getIV()");
    return &_iv[icand * _nx * _ny];
}

double HscPsfBase::getX(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getX()");
    return _current_xy[2*icand];
}

double HscPsfBase::getY(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getY()");
    return _current_xy[2*icand + 1];
}

double HscPsfBase::getXini(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getXini()");
    return _initial_xy[2*icand];
}

double HscPsfBase::getYini(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::getYini()");
    return _initial_xy[2*icand + 1];
}

void HscPsfBase::setBad(int icand)
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::setBad()");
    _bad[icand] = true;
}

PTR(HscCandidateSet) HscPsfBase::getTrimmedCandidateSet() const
{
    PTR(HscCandidateSet) ret = boost::make_shared<HscCandidateSet>(_cs->getMaskBits(), _nx, _ny);

    for (int icand = 0; icand < _ncand; icand++) {
        if (!_bad[icand])
            ret->add(_cs, icand, _current_xy[2*icand], _current_xy[2*icand+1]);
    }

    if (ret->getNcand() == 0)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "all candidates were marked as bad in HscPsfBase::getTrimmedCandidateSet()");

    std::cerr << "HscPsf: candidate set trimmed " << _ncand << " -> " << ret->getNcand() << std::endl;
    return ret;
}

void HscPsfBase::_updateCurrentPsfImg(int icand, const double *psf_img)
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscPsfBase::_updateCurrentPsfImg()");

    assert(psf_img != NULL);

    double chi2 = 0.0;
    for (int i = 0; i < _nx*_ny; i++)
        chi2 += _iv[icand*_nx*_ny+i] * (_im[icand*_nx*_ny+i] - psf_img[i]) * (_im[icand*_nx*_ny+i] - psf_img[i]);

    memcpy(&_current_psfimg[icand*_nx*_ny], psf_img, _nx * _ny * sizeof(double));
    _residual_chi2[icand] = chi2;
}


double HscPsfBase::get_total_reduced_chi2() const
{
    double num = 0.0;
    for (int i = 0; i < _ncand; i++)
        num += _residual_chi2[i];

    int den = 0;
    for (int i = 0; i < _ncand; i++)
        den += _ndof[i];

    return num/den;
}


void HscPsfBase::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "HscPsfBase::eval() not implemented by subclass");
}



PTR(afw::detection::Psf) HscPsfBase::clone() const
{
    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "HscPsfBase::clone() not implemented yet");   // FIXME
}

afw::geom::Point2D HscPsfBase::getAveragePosition() const
{
    return afw::geom::Point2D(_xmean, _ymean);
}

// in the same coordinate system as the pixelized image
PTR(HscPsfBase::Image) HscPsfBase::doComputeImage(afw::geom::Point2D const &position, afw::image::Color const &color) const
{
    int n = 2*_nside + 1;
    std::vector<double> img(n*n);

    // FIXME is this what I should be using?
    int x0 = (int)(position.getX() - _nside);
    int y0 = (int)(position.getY() - _nside);
    this->eval(n, n, (double)x0, (double)y0, &img[0], position.getX(), position.getY());

    PTR(Image) ret = boost::make_shared<Image> (n, n);
    ret->setXY0(x0, y0);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            (*ret)(i,j) = img[i*n+j];
    
    return ret;
}


PTR(HscPsfBase::Image) HscPsfBase::doComputeKernelImage(afw::geom::Point2D const &position, afw::image::Color const &color) const
{
    int n = 2*_nside + 1;
    std::vector<double> img(n*n);

    double x0 = position.getX() - _nside;
    double y0 = position.getY() - _nside;
    this->eval(n, n, x0, y0, &img[0], position.getX(), position.getY());

    PTR(Image) ret = boost::make_shared<Image> (n, n);
    ret->setXY0(-_nside, -_nside);
    
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            (*ret)(i,j) = img[i*n+j];
    
    return ret;
}


// static member function
void HscPsfBase::make_shear_matrix(double &axx, double &axy, double &ayy, double gamma1, double gamma2, double kappa)
{
    double ek = std::exp(kappa);
    double g = std::sqrt(gamma1*gamma1 + gamma2*gamma2);
    double c = std::cosh(g);
    double s = (g > 1.0e-10) ? (std::sinh(g)/g) : 1.0;
    
    axx = ek * (c + s*gamma1);
    axy = ek * (s*gamma2);
    ayy = ek * (c - s*gamma1);
}


// static member function
double HscPsfBase::fit_basis_images(double *out_ampl, int nbf, int nxy, const double *iv, const double *im, const double *basis_images)
{
    assert(basis_images != NULL);
    assert(iv != NULL);
    assert(im != NULL);
    assert(nbf > 0);
    assert(nxy > 0);

    Eigen::MatrixXd a(nbf,nbf);
    Eigen::VectorXd b(nbf);

    for (int i = 0; i < nbf; i++) {
	for (int j = 0; j <= i; j++) {
            double t = 0.0;
            for (int k = 0; k < nxy; k++)
                t += iv[k] * basis_images[i*nxy+k] * basis_images[j*nxy+k];
	    a(i,j) = a(j,i) = t;
        }

        double t = 0.0;
        for (int k = 0; k < nxy; k++)
            t += iv[k] * im[k] * basis_images[i*nxy+k];
        b(i) = t;
    }
    
    Eigen::VectorXd ainv_b = a.llt().solve(b);

    if (out_ampl != NULL) {
        for (int i = 0; i < nbf; i++)
            out_ampl[i] = ainv_b(i);
    }

    double residual_chi2 = 0.0;
    for (int k = 0; k < nxy; k++)
        residual_chi2 += iv[k] * im[k] * im[k];
    for (int i = 0; i < nbf; i++)
        residual_chi2 -= b[i] * ainv_b(i);

    return residual_chi2;
}


// anonymous namespace of helper functions for HscPsfBase::lanczos_offset_2d()
namespace {
    inline double lanczos_kernel(int order, double x)
    {
        if (x*x < 1.0e-20)
            return 1.0;

        if (x*x >= order*order)
            return 0.0;
        
        double t = std::sin(M_PI*x) / (M_PI*x);
        double u = std::sin(M_PI*x/order) / (M_PI*x/order);
        return t * u;
    }


    inline void set_lanczos_delimeters(int order, double x, int n, bool zero_pad, int &i0, int &i1)
    {
        assert(order >= 1);
        assert(n >= 0);
        
        double x0 = x - order + 1 + 1.0e-13;   // i0 = floor(x0)
        double x1 = x + order + 1 - 1.0e-13;   // i1 = floor(x1)
        
        if (x0 < 0.0) {
            if (!zero_pad)
                throw LSST_EXCEPT(pex::exceptions::OutOfRangeException, "lanczos interpolation out of range");
            i0 = 0;
            i1 = (x1 >= 0.0) ? std::min(int(x1),n) : 0;
            return;
        }
        
        i0 = int(x0);
        i1 = int(x1);
        
        if (i1 > n) {
            if (!zero_pad)
                throw LSST_EXCEPT(pex::exceptions::OutOfRangeException, "lanczos interpolation out of range");
            if (i0 > n)
                i0 = n;
            i1 = n;
        }
    }


    void lanczos_x_offset(int order, double x0, int nx_out, int nx_in, int ny, double *out, const double *in, bool zero_pad)
    {
        memset(out, 0, nx_out * ny * sizeof(*out));
        
        int ix0 = (x0 >= 0.0) ? int(x0) : (-int(-x0) - 1);
        double dx0 = x0 - (double)ix0;
        
        if (!zero_pad && ((ix0-order+1 < 0) || (nx_out+ix0+order > nx_in)))
            throw LSST_EXCEPT(pex::exceptions::OutOfRangeException, "lanczos interpolation out of range and zero_pad=false");
        
        std::vector<double> lv(2*order);
        for (int i = -order; i < order; i++)
            lv[i+order] = lanczos_kernel(order, dx0+i);
        
        for (int i = 0; i < nx_out; i++) {
            int j0 = std::max(0, (i+ix0) - (order-1));
            int j1 = std::min(nx_in, (i+ix0) + (order+1));
            
            for (int j = j0; j < j1; j++) {
                double t = lv[i+ix0-j+order];   // = L((i+ix0+dx0) - j)
                
                for (int k = 0; k < ny; k++)
                    out[i*ny+k] += t * in[j*ny+k];
            }
        }
    }
};   // anonymous namespace of helper functions for HscPsfBase::lanczos_offset_2d()


// static member function
void HscPsfBase::lanczos_offset_2d(int order, double x0, double y0, 
                                   int nx_out, int ny_out, double *out, int out_stride, 
                                   int nx_in, int ny_in, const double *in, int in_stride, 
                                   bool zero_pad, bool accumulate)
{
    std::vector<double> scratch(4*order);

    if ((order < 1) || (nx_out <= 0) || (ny_out <= 0) || (out == NULL) || (out_stride < ny_out) || (nx_in <= 0) || (ny_in <= 0) || (in == NULL) || (in_stride < ny_in))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters in lanczos_offset_2d()");

    std::vector<double> tmp(nx_in * ny_in);
    for (int i = 0; i < nx_in; i++)
	for (int j = 0; j < ny_in; j++)
	    tmp[i*ny_in+j] = in[i*in_stride+j];

    std::vector<double> tmp2(nx_out * ny_in);
    lanczos_x_offset(order, x0, nx_out, nx_in, ny_in, &tmp2[0], &tmp[0], zero_pad);

    // transpose
    tmp = std::vector<double>(ny_in * nx_out);
    for (int i = 0; i < nx_out; i++)
	for (int j = 0; j < ny_in; j++)
	    tmp[j*nx_out + i] = tmp2[i*ny_in + j];

    tmp2 = std::vector<double>(ny_out * nx_out);
    lanczos_x_offset(order, y0, ny_out, ny_in, nx_out, &tmp2[0], &tmp[0], zero_pad);

    for (int i = 0; i < nx_out; i++) {
	if (accumulate) {
	    for (int j = 0; j < ny_out; j++)
		out[i*out_stride + j] += tmp2[j*nx_out + i];
	}
	else {
	    for (int j = 0; j < ny_out; j++)
		out[i*out_stride + j] = tmp2[j*nx_out + i];
	}
    }
}


}}}}   // namespace lsst::meas::extensions::hscpsf
