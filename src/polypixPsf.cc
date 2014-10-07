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


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}   // emacs pacifier
#endif

static inline double square(double x) { return x*x; }


PolypixPsf::PolypixPsf(CONST_PTR(HscCandidateSet) cs, int nside, int psf_size, double psfstep, CONST_PTR(HscSpatialModelBase) spatialModel, double fwhm, double backnoise2, double gain)
    : HscPsfBase(cs,nside)
{ 
    this->_construct(psf_size, psfstep, spatialModel, fwhm, backnoise2, gain);
}


PolypixPsf::PolypixPsf(CONST_PTR(HscCandidateSet) cs, CONST_PTR(PolypixPsf) base)
    : HscPsfBase(cs, base->_nside)
{
    this->_construct(base->_psf_nx, base->_psfstep, base->_spatialModel, base->_fwhm, base->_backnoise2, base->_gain);

    _xmin = base->_xmin;
    _xmax = base->_xmax;
    _xmean = base->_xmean;
    _ymin = base->_ymin;
    _ymax = base->_ymax;
    _ymean = base->_ymean;
}


void PolypixPsf::psf_make(double prof_accuracy, double regul)
{
    int npix = _psf_nx * _psf_ny;
    std::vector<double> image(_ncand * npix);
    std::vector<double> weight(_ncand * npix);

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        this->_upsample(&image[icand*npix], &_im[icand*_nx*_ny], dx, dy);

        for (int ipix = 0; ipix < npix; ipix++)
            image[icand*npix + ipix] /= _flux[icand];
            
        for (int ipix = 0; ipix < npix; ipix++) {
            double val = image[icand*npix + ipix];
            double norm2 = _flux[icand] * _flux[icand];
            double profaccu2 = prof_accuracy * prof_accuracy * norm2;
            double noise2 = _backnoise2 + profaccu2 * val * val;

            if (val > 0.0 && _gain > 0.0)
                noise2 += val/_gain;

            weight[icand*npix + ipix] = norm2/noise2;
        }
    }

    std::vector<double> a(_ncand);
    std::vector<double> b(_ncand);

    for (int ipix = 0; ipix < npix; ipix++) {
        for (int icand = 0; icand < _ncand; icand++) {
            a[icand] = weight[icand*npix + ipix];
            b[icand] = weight[icand*npix + ipix] * image[icand*npix + ipix];
        }
        _spatialModel->optimize(&_tcomp[ipix*_ncoeffs], _ncand, &_current_xy[0], &a[0], &b[0], regul);
    }
}


PTR(HscCandidateSet) PolypixPsf::psf_clean(double prof_accuracy)
{
    std::vector<double> chi2 = this->_make_cleaning_chi2(prof_accuracy);

    std::vector<double> chi(_ncand);
    for (int icand = 0; icand < _ncand; icand++)
        chi[icand] = sqrt(chi2[icand]);

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
        if (chi2[icand] <= chi2max)
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

    for (int ix = 0; ix < _psf_nx; ix++) {
        for (int iy = 0; iy < _psf_ny; iy++) {
            double r2 = (ix-xc)*(ix-xc) + (iy-yc)*(iy-yc);
            double t = 1.0;            

            if (r2 >= rmax2)
                t = 0.0;
            else if (r2 > rmin2)
                t = ((rmax2-r2) / dr2);
            
            for (int ic = 0; ic < _ncoeffs; ic++)
                _tcomp[ix*_psf_ny*_ncoeffs + iy*_ncoeffs + ic] *= t;
        }
    }
}


// follows PsfexPsf::_doComputeImage() in meas_extensions_psfex
void PolypixPsf::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    double xy[2];
    xy[0] = x;
    xy[1] = y;

    std::vector<double> fullresIm(_psf_nx * _psf_ny);      
    _spatialModel->eval(&fullresIm[0], 1, xy, _psf_nx * _psf_ny, &_tcomp[0]);

    double dx = x - x0 - 0.5*(nx_out-1);
    double dy = y - y0 - 0.5*(ny_out-1);

    this->_downsample(out, nx_out, ny_out, &fullresIm[0], dx, dy);

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


// -------------------------------------------------------------------------------------------------



// private helper function for constructors
void PolypixPsf::_construct(int psf_size, double psfstep, CONST_PTR(HscSpatialModelBase) spatialModel, double fwhm, double backnoise2, double gain)
{
    _fwhm = fwhm;
    _backnoise2 = backnoise2;
    _gain = gain;

    _spatialModel = spatialModel;
    _ncoeffs = spatialModel->getNcoeffs();

    if (_ncand <= _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "too few spatial candidates in PolypixPsf constructor");

    _psf_nx = psf_size;
    _psf_ny = psf_size;

    _flux.resize(_ncand);
    _tcomp.resize(_psf_nx * _psf_ny * _ncoeffs, 0.0);
    _psfstep = (psfstep > 0.0) ? psfstep : (_fwhm/2.35 * 0.5);

    // FIXME understand this!
    if (_psfstep > 1.0)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "I don't currently understand what to do if psfstep > 1");

    // FIXME revisit this!
    for (int icand = 0; icand < _ncand; icand++) {
        double flux = this->getCandSet()->getFlux(icand);
        if (flux <= 0.0)
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "flux < 0 in PolypixPsf constructor");
        _flux[icand] = flux;
    }
}


void PolypixPsf::_downsample(double *out, int nx_out, int ny_out, const double *in, double dx, double dy) const
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

void PolypixPsf::_upsample(double *out, const double *in, double dx, double dy) const
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

    for (int i = 0; i < _psf_nx; i++) {
        for (int j = 0; j < _psf_ny; j++) {
	    double x = x0 + i*_psfstep;
	    double y = y0 + j*_psfstep;

            //
            // FIXME: following psfex, but seems fishy (better to enlarge the candidate image by
            // the interpolation kernel width?)
            //
            if (x >= 0.0 && x <= _nx-1 && y >= 0.0 && y <= _ny-1)
                out[i*_psf_ny+j] = lanczos_interpolate_2d(order, x, y, _nx, _ny, in, _ny, &scratch[0], true, true);
            else
                out[i*_psf_ny+j] = 0.0;
	}
    }
}


//
// Preserves a lot of psfex logic that I can't say I understand!
//
std::vector<double> PolypixPsf::_make_cleaning_chi2(double prof_accuracy)
{
    // -------------------- initalize vigweight --------------------

    std::vector<double> vigweight(_ncand * _nx * _ny);
    const double vw_accuracy = 0.01;   // FIXME hardcoded

    for (int i = 0; i < _ncand * _nx * _ny; i++) {
        if (_iv[i] <= 0.0)   // replaces (_im[i] < -BIG) condition in psfex
            vigweight[i] = 0.0;
        else {
            double noise2 = _backnoise2 + (vw_accuracy * vw_accuracy * _im[i] * _im[i]);
            if (_im[i] > 0.0 && _gain > 0.0)
                noise2 += _im[i] / _gain;
            vigweight[i] = 1.0 / noise2;
        }
    }

    const bool accuflag = (prof_accuracy > 1.0e-30);

    std::vector<double> loc(_psf_nx * _psf_ny);
    std::vector<double> dresi(_nx*_ny, 0.0);

    std::vector<double> vigresi(_ncand * _nx * _ny, 0.0);
    std::vector<double> ret(_ncand);

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        // psf_build
        _spatialModel->eval(&loc[0], 1, &_current_xy[2*icand], _psf_nx * _psf_ny, &_tcomp[0]);

        // temporarily set vigresi = psf in vignette coordinates (not postmultiplied by flux)
        this->_downsample(&vigresi[icand*_nx*_ny], _nx, _ny, &loc[0], dx, dy);

        // -------------------- fit flux --------------------

        double xi2 = 0.0;
        double xyi = 0.0;

        for (int s = icand*_nx*_ny; s < (icand+1)*_nx*_ny; s++) {
            xi2 += vigweight[s] * vigresi[s] * vigresi[s];
            xyi += vigweight[s] * vigresi[s] * _im[s];
        }

        double norm = (xi2 > 0.0) ? (xyi/xi2) : _flux[icand];

        // -------------------- subtract psf model and compute chi^2 --------------------

        double psf_extraccu2 = prof_accuracy*prof_accuracy*norm*norm;
        double rmax2 = square(_psfstep * (double)(std::min(_psf_nx,_psf_ny)/2));
        double xc = _current_xy[2*icand] - _xy0[2*icand];
        double yc = _current_xy[2*icand+1] - _xy0[2*icand+1];
        double chi2 = 0.0;
        int nchi2 = 0;

        for (int ix = 0; ix < _nx; ix++) {
            for (int iy = 0; iy < _ny; iy++) {
                int s = icand*_nx*_ny + ix*_ny + iy;

                double x = ix-xc;
                double y = iy-yc;
                double wval = vigweight[s];

                if (wval > 0.0) {
                    if (accuflag)
                        wval = 1.0 / (1.0/wval + psf_extraccu2 * vigresi[s] * vigresi[s]);
                    
                    // at this point in the code, vigresi = actual residual
                    vigresi[s] = _im[s] - norm * vigresi[s];

                    if (x*x+y*y < rmax2) {
                        chi2 += wval * square(vigresi[s]);
                        nchi2++;
                    }
                }
            }
        }

        ret[icand] = (nchi2 > 1)? (chi2/(nchi2-1)) : chi2;
    }    

    return ret;
}


}}}}   // namespace lsst::meas::extensions::hscpsf
