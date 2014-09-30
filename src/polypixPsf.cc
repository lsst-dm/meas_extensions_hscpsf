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
 

#include <cmath>
#include <cstring>
#include "lsst/meas/extensions/hscpsf/hscPsf.h"

#define BIG 1.0e30

static inline double square(double x) { return x*x; }


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}   // emacs pacifier
#endif


// private helper function for constructors
void PolypixPsf::_construct(int spatialOrder, double fwhm, double backnoise2, double gain)
{
    if (spatialOrder < 0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "spatialOrder < 0 in PolypixPsf constructor");

    _spatialOrder = spatialOrder;
    _ncoeffs = ((spatialOrder+1) * (spatialOrder+2)) / 2;
    _spatialModel = boost::make_shared<HscSpatialModelPolynomial>(_spatialOrder, _xmin, _xmax, _ymin, _ymax);

    _fwhm = fwhm;
    _backnoise2 = backnoise2;
    _gain = gain;

    if (_ncand <= _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "too few spatial candidates in PolypixPsf constructor");

    _psf_nx = _nx;
    _psf_ny = _ny;

    _norm.resize(_ncand);
    _vigweight.resize(_ncand * _nx * _ny);
    _comp.resize(_ncoeffs * _psf_nx * _psf_ny, 0.0);
    _vigresi.resize(_ncand * _nx * _ny, 0.0);
    _vigchi.resize(_ncand * _nx * _ny, 0.0);    
    _chi2.resize(_ncand, 0.0);

    _psfstep = _fwhm/2.35 * 0.5;

    // FIXME understand this!
    if (_psfstep > 1.0)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "I don't currently understand what to do if pixstep > 1");

    // FIXME revisit this!
    for (int icand = 0; icand < _ncand; icand++) {
        double flux = this->getCandSet()->getFlux(icand);
        if (flux <= 0.0)
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "flux < 0 in PolypixPsf constructor");
        _norm[icand] = flux;
    }

    //
    // vigweight initialization offset modeled on psfex sample_utils.c make_weights()
    //
    const double prof_accuracy = 0.01;   // FIXME hardcoded

    for (int i = 0; i < _ncand * _nx * _ny; i++) {
        if (_iv[i] <= 0.0)   // replaces (_im[i] < -BIG) condition in psfex
            _vigweight[i] = 0.0;
        else {
            double noise2 = _backnoise2 + (prof_accuracy * prof_accuracy * _im[i] * _im[i]);
            if (_im[i] > 0.0 && _gain > 0.0)
                noise2 += _im[i] / _gain;
            _vigweight[i] = 1.0 / noise2;
        }
    }
}


PolypixPsf::PolypixPsf(CONST_PTR(HscCandidateSet) cs, int nside, int spatialOrder, double fwhm, double backnoise2, double gain)
    : HscPsfBase(cs,nside)
{ 
    this->_construct(spatialOrder, fwhm, backnoise2, gain);
}


PolypixPsf::PolypixPsf(CONST_PTR(HscCandidateSet) cs, CONST_PTR(PolypixPsf) base)
    : HscPsfBase(cs, base->_nside)
{
    this->_construct(base->_spatialOrder, base->_fwhm, base->_backnoise2, base->_gain);

    _xmin = base->_xmin;
    _xmax = base->_xmax;
    _xmean = base->_xmean;
    _ymin = base->_ymin;
    _ymax = base->_ymax;
    _ymean = base->_ymean;
}


void PolypixPsf::downsample(double *out, int nx_out, int ny_out, const double *in, double dx, double dy) const
{
    const int order = 3;
    std::vector<double> scratch(4*order);

    if ((nx_out % 2 == 0) || (ny_out % 2 == 0) || (_psf_nx % 2 == 0) || (_psf_ny % 2 == 0)) {
        //
        // The psfex behavior in the even case looks fishy to me, so I decided to throw an
        // exception to trap the even case, in order to rethink things if this case actually
        // arises in practice..
        //
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, 
                          "Image size is even in PolypixPsf::downsample(); this is currently unsupported");
    }

    memset(out, 0, nx_out * ny_out * sizeof(double));

    double x0 = (_psf_nx/2) - (nx_out/2)/_psfstep - dx/_psfstep;
    double y0 = (_psf_ny/2) - (ny_out/2)/_psfstep - dy/_psfstep;

    for (int i = 0; i < nx_out; i++) {
        for (int j = 0; j < ny_out; j++) {
	    double x = x0 + i/_psfstep;
	    double y = y0 + j/_psfstep;

	    if (x >= 0.0 && x <= _psf_nx-1 && y >= 0.0 && y <= _psf_ny-1)
		out[i*ny_out+j] = lanczos_interpolate_2d(order, x, y, _psf_nx, _psf_ny, in, _psf_ny, &scratch[0], true, true);
	}
    }
}

void PolypixPsf::upsample(double *out, const double *in, double dx, double dy) const
{
    const int order = 3;
    std::vector<double> scratch(4*order);

    if ((_nx % 2 == 0) || (_ny % 2 == 0) || (_psf_nx % 2 == 0) || (_psf_ny % 2 == 0)) {
        //
        // The psfex behavior in the even case looks fishy to me, so I decided to throw an
        // exception to trap the even case, in order to rethink things if this case actually
        // arises in practice..
        //
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, 
                          "Image size is even in PolypixPsf::upsample(); this is currently unsupported");
    }

    double x0 = (_nx/2) - (_psf_nx/2)*_psfstep + dx;
    double y0 = (_ny/2) - (_psf_ny/2)*_psfstep + dy;

    double x1 = x0 + (_psf_nx-1) * _psfstep;
    double y1 = y0 + (_psf_ny-1) * _psfstep;

    if ((x0 < 0) || (x1 > _nx-1) || (y0 < 0) || (y1 > _ny-1)) {
        //
        // Currently treated as an error (rethink if this case actually arises in practice)
        //
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException,
                          "PSF candidate image is too small in PolypixPsf::upsample()");
    }

    for (int i = 0; i < _psf_nx; i++) {
        for (int j = 0; j < _psf_ny; j++) {
	    double x = x0 + i*_psfstep;
	    double y = y0 + j*_psfstep;
            out[i*_psf_ny+j] = lanczos_interpolate_2d(order, x, y, _nx, _ny, in, _ny, &scratch[0], true, true);
	}
    }
}


void PolypixPsf::psf_make(double prof_accuracy)
{
    // FIXME rethink this
    const double regul = 1000.0;

    std::vector<double> image(_ncand * _psf_nx * _psf_ny);
    std::vector<double> weight(_ncand * _psf_nx * _psf_ny);

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        upsample(&image[icand*_psf_nx*_psf_ny], &_im[icand*_nx*_ny], dx, dy);

        for (int i = 0; i < _nx*_ny; i++)
            image[icand*_psf_nx*_psf_ny + i] /= _norm[icand];
            
        for (int i = 0; i < _nx*_ny; i++) {
            double val = image[icand*_psf_nx*_psf_ny + i];
            double norm2 = _norm[icand] * _norm[icand];
            double profaccu2 = prof_accuracy * prof_accuracy * norm2;
            double noise2 = _backnoise2 + profaccu2 * val * val;

            if (val > 0.0 && _gain > 0.0)
                noise2 += val/_gain;

            weight[icand*_nx*_ny + i] = norm2/noise2;
        }
    }

    std::vector<double> a(_ncand);
    std::vector<double> b(_ncand);
    std::vector<double> coeff(_ncoeffs);
    int npix = _psf_nx * _psf_ny;

    for (int ipix = 0; ipix < npix; ipix++) {
        for (int icand = 0; icand < _ncand; icand++) {
            a[icand] = weight[icand*npix + ipix];
            b[icand] = weight[icand*npix + ipix] * image[icand*npix + ipix];
        }

        _spatialModel->optimize(&coeff[0], _ncand, &_current_xy[0], &a[0], &b[0], regul);

        for (int icoeff = 0; icoeff < _ncoeffs; icoeff++)
            _comp[icoeff*npix + ipix] = coeff[icoeff];
    }
}


PTR(HscCandidateSet) PolypixPsf::psf_clean(double prof_accuracy)
{
    this->psf_makeresi(prof_accuracy);

    std::vector<double> chi(_ncand);
    for (int icand = 0; icand < _ncand; icand++)
        chi[icand] = sqrt(_chi2[icand]);

    std::sort(chi.begin(), chi.end());
    int lo = 0;
    int hi = _ncand;

    double chimean, chisig, chimed;
    double chisig1 = 0.0;   // used to save value from previous iteration

    //
    // Produce k-sigma-clipped statistics
    //
    for (int iter = 0; iter < 100; iter++) {
        chimean = chisig = 0.0;
        for (int i = lo; i < hi; i++) {
            chimean += chi[i];
            chisig += chi[i]*chi[i];
        }

        int n = hi-lo;
        chimean /= n;
        chisig = sqrt((chisig-chimean*chimean*n) / std::max(n-1,1));
        chimed = (n % 2) ? chi[lo+n/2] : (chi[lo+n/2-1] + chi[lo+n/2])/2.0;

        // STL binary search boilerplate
        lo = std::distance(chi.begin(), std::upper_bound(chi.begin(), chi.end(), chimed-3*chisig));
        hi = std::distance(chi.begin(), std::upper_bound(chi.begin(), chi.end(), chimed+3*chisig));
        assert(lo < hi);
    
        if (chisig < 0.1)
            break;
        if (chisig1 > 0.0 && (fabs(chisig/chisig1-1.0) <= 1.0e-4))
            break;

        chisig1 = chisig;
    }

    PTR(HscCandidateSet) ret = boost::make_shared<HscCandidateSet>(_cs->getMaskBits(), _nx, _ny);
    double chi2max = (chimed + 3*chisig) * (chimed + 3*chisig);

    for (int icand = 0; icand < _ncand; icand++)
        if (_chi2[icand] <= chi2max)
            ret->add(_cs, icand);

    return ret;
}


void PolypixPsf::psf_clip()
{
    double xc = (double)(_psf_nx/2);
    double yc = (double)(_psf_ny/2);
    double rmax2 = std::min(xc,yc) + 0.5;
    double dr2 = _fwhm / _psfstep;
    dr2 = std::max(dr2, 1.0);

    // FIXME bug in psfex here?
    if (dr2 >= rmax2)
        dr2 = rmax2/2.0;

    double rmin2 = rmax2 - dr2;
    rmin2 *= rmin2;
    rmax2 *= rmax2;
    dr2 = rmax2 - rmin2;

    for (int ic = 0; ic < _ncoeffs; ic++) {
        for (int ix = 0; ix < _psf_nx; ix++) {
            for (int iy = 0; iy < _psf_ny; iy++) {
                double r2 = (ix-xc)*(ix-xc) + (iy-yc)*(iy-yc);
            
                if (r2 >= rmax2)
                    _comp[ic*_psf_nx*_psf_ny + ix*_psf_ny + iy] = 0.0;
                else if (r2 > rmin2)
                    _comp[ic*_psf_nx*_psf_ny + ix*_psf_ny + iy] *= (rmax2-r2) / dr2;
            }
        }
    }
}


void PolypixPsf::psf_build(double *loc, const double *xy) const
{
    int npix = _psf_nx * _psf_ny;
    memset(loc, 0, npix * sizeof(double));

    // FIXME oops, _comp is transposed
    std::vector<double> _tcomp(npix * _ncoeffs);
    for (int ic = 0; ic < _ncoeffs; ic++)
        for (int ipix = 0; ipix < npix; ipix++)
            _tcomp[ipix*_ncoeffs+ic] = _comp[ic*npix+ipix];

    _spatialModel->eval(loc, 1, xy, npix, &_tcomp[0]);
}


void PolypixPsf::psf_makeresi(double prof_accuracy)
{
    const bool accuflag = (prof_accuracy > 1.0/BIG);

    std::vector<double> loc(_psf_nx * _psf_ny);
    std::vector<double> dresi(_nx*_ny, 0.0);

    memset(&_vigresi[0], 0, _vigresi.size() * sizeof(double));
    memset(&_vigchi[0], 0, _vigchi.size() * sizeof(double));

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        this->psf_build(&loc[0], &_current_xy[2*icand]);

        // temporarily set vigresi = psf in vignette coordinates (not postmultiplied by flux)
        downsample(&_vigresi[icand*_nx*_ny], _nx, _ny, &loc[0], dx, dy);

        // -------------------- fit flux --------------------

        double xi2 = 0.0;
        double xyi = 0.0;

        for (int s = icand*_nx*_ny; s < (icand+1)*_nx*_ny; s++) {
            xi2 += _vigweight[s] * _vigresi[s] * _vigresi[s];
            xyi += _vigweight[s] * _vigresi[s] * _im[s];
        }

        double norm = (xi2 > 0.0) ? (xyi/xi2) : _norm[icand];

        // -------------------- subtract psf model and compute chi^2 --------------------

        double psf_extraccu2 = prof_accuracy*prof_accuracy*norm*norm;
        double rmax2 = square(_psfstep * (double)(std::min(_nx,_ny)/2));
        double xc = _current_xy[2*icand] - _xy0[2*icand];
        double yc = _current_xy[2*icand+1] - _xy0[2*icand+1];
        double chi2 = 0.0;
        int nchi2 = 0;

        for (int ix = 0; ix < _nx; ix++) {
            for (int iy = 0; iy < _ny; iy++) {
                int s = icand*_nx*_ny + ix*_ny + iy;

                double x = ix-xc;
                double y = iy-yc;
                double wval = _vigweight[s];

                if (wval > 0.0) {
                    if (accuflag)
                        wval = 1.0 / (1.0/wval + psf_extraccu2 * _vigresi[s] * _vigresi[s]);
                    
                    // at this point in the code, vigresi = actual residual
                    _vigresi[s] = _im[s] - norm * _vigresi[s];

                    if (x*x+y*y < rmax2) {
                        _vigchi[s] = wval * square(_vigresi[s]);
                        chi2 += _vigchi[s];
                        nchi2++;
                    }
                }
            }
        }

        _chi2[icand] = (nchi2 > 1)? (chi2/(nchi2-1)) : chi2;
    }    
}


// -------------------------------------------------------------------------------------------------



// follows PsfexPsf::_doComputeImage() in meas_extensions_psfex
void PolypixPsf::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    double xy[2];
    xy[0] = x;
    xy[1] = y;

    std::vector<double> fullresIm(_psf_nx * _psf_ny);      
    psf_build(&fullresIm[0], xy);

    double dx = x - x0 - 0.5*(nx_out-1);
    double dy = y - y0 - 0.5*(ny_out-1);

    downsample(out, nx_out, ny_out, &fullresIm[0], dx, dy);

    //
    // Normalize to sum 1 
    // (FIXME put this somewhere more general?)
    //
    double acc = 0.0;
    for (int i = 0; i < nx_out * ny_out; i++)
        acc += out[i];
    for (int i = 0; i < nx_out * ny_out; i++)
        out[i] /= acc;
}


}}}}   // namespace lsst::meas::extensions::hscpsf
