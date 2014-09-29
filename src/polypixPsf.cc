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

#define EPS 1.0e-4
#define BIG 1.0e30

static inline double square(double x) { return x*x; }


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}   // emacs pacifier
#endif


static void my_resample(const double *pix1, int w1, int h1, double *pix2, int w2, int h2,
                        double dx, double dy, double step2, double step1)
{
    assert(step1 == 1.0);

    const int order = 3;  // INTERPFAC

    std::vector<double> scratch(4*order);

    // Note integer division here!
    // In my opinion, this is a bug in psfex, but we reproduce it here...
    double x0 = (w1/2) - (w2/2)*step2 + dx;
    double y0 = (h1/2) - (h2/2)*step2 + dy;

    for (int i = 0; i < w2; i++) {
        for (int j = 0; j < h2; j++) {
	    double x = x0 + i*step2;
	    double y = y0 + j*step2;

	    if (x >= 0.0 && x <= w1-1 && y >= 0.0 && y <= h1-1)
		pix2[i*h2+j] = HscPsfBase::lanczos_interpolate_2d(order, x, y, w1, h1, pix1, h1, &scratch[0], true, true);
	    else
		pix2[i*h2+j] = 0.0;
	}
    }
}


// private helper function for constructors
void PolypixPsf::_construct(int spatialOrder, double fwhm, double backnoise2, double gain)
{
    if (spatialOrder < 0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "spatialOrder < 0 in PolypixPsf constructor");

    _spatialOrder = spatialOrder;
    _ncoeffs = ((spatialOrder+1) * (spatialOrder+2)) / 2;
    _fwhm = fwhm;
    _backnoise2 = backnoise2;
    _gain = gain;


    if (_ncand <= _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "too few spatial candidates in PolypixPsf constructor");

    std::cerr << "FIXME !!! setting (psf_nx,psf_ny) = (cand_nx,cand_ny), revisit !!!\n";
    _psf_nx = _nx;
    _psf_ny = _ny;

    _norm.resize(_ncand);
    _contextoffset.resize(2, 0.0);   // note: allocated here, but initialized in constructor body
    _contextscale.resize(2, 0.0);    // note: allocated here, but initialized in constructor body
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
    
    // FIXME revisit
    _contextoffset[0] = (_xmin + _xmax) / 2.0;
    _contextoffset[1] = (_ymin + _ymax) / 2.0;
    _contextscale[0] = (_xmax - _xmin);
    _contextscale[1] = (_ymax - _ymin);    
}


PolypixPsf::PolypixPsf(CONST_PTR(HscCandidateSet) cs, CONST_PTR(PolypixPsf) base)
    : HscPsfBase(cs, base->_nside)
{
    this->_construct(base->_spatialOrder, base->_fwhm, base->_backnoise2, base->_gain);

    _contextoffset = base->_contextoffset;
    _contextscale = base->_contextscale;

    _xmin = base->_xmin;
    _xmax = base->_xmax;
    _xmean = base->_xmean;
    _ymin = base->_ymin;
    _ymax = base->_ymax;
    _ymean = base->_ymean;
}


void PolypixPsf::psf_make(double prof_accuracy)
{
    // FIXME rethink this
    const double regul = 1000.0;

    std::vector<double> image(_ncand * _psf_nx * _psf_ny);
    std::vector<double> weight(_ncand * _psf_nx * _psf_ny);
    std::vector<double> pos(_ncand * 2);

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        my_resample(&_im[icand*_nx*_ny], _nx, _ny,
                    &image[icand*_psf_nx*_psf_ny], _psf_nx, _psf_ny,
                    dx, dy, _psfstep, 1.0);   // FIXME note in psfex code, last arg is max(_psfstep,1.0)

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

        pos[2*icand] = (_current_xy[2*icand] - _contextoffset[0]) / _contextscale[0];
        pos[2*icand+1] = (_current_xy[2*icand+1] - _contextoffset[1]) / _contextscale[1];
    }

    std::vector<double> basis(_ncand * _ncoeffs);
    this->poly_eval_basis_functions(&basis[0], &pos[0]);

    std::vector<double> pstack(_ncand);
    std::vector<double> wstack(_ncand);
    std::vector<double> coeff(_ncoeffs);
    int npix = _psf_nx * _psf_ny;

    for (int ipix = 0; ipix < npix; ipix++) {
        for (int icand = 0; icand < _ncand; icand++) {
            pstack[icand] = image[icand*npix + ipix];
            wstack[icand] = weight[icand*npix + ipix];
        }

        this->poly_fit(&coeff[0], &pstack[0], &wstack[0], &basis[0], regul);

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

    // Produce k-sigma-clipped statistics
    double locut = -BIG;
    double hicut = BIG;
    double chisig = BIG;
    int n = _ncand;

    for (int iter = 0; iter < 100; iter++) {
        double chimed = fast_median(&chi[0], n);
        double chimean = 0.0;
        double chivar = 0.0;

        int n2 = 0;
        for (int i = 0; i < n; i++) {
            if (chi[i] > locut && chi[i] < hicut) {
                chimean += chi[i];
                chivar += square(chi[i]);
                chi[n2++] = chi[i];
            }
        }

        double chisig1 = chisig;
        n = n2;

        chimean /= (double)n;
        chisig = sqrt((chivar-chimean*chimean*n) / std::max(n-1,1));
        locut = chimed - 3.0*chisig;
        hicut = chimed + 3.0*chisig;
    
        if (chisig < 0.1 || fabs(chisig/chisig1-1.0) <= EPS)
            break;
    }

    double chi2max = square(hicut);
    std::cerr << "psf_clean(" << prof_accuracy << "): chi2max=" << chi2max << std::endl;

    PTR(HscCandidateSet) ret = boost::make_shared<HscCandidateSet>(_cs->getMaskBits(), _nx, _ny);

#if 0
    //
    // Drop candidates which exceed the chi^2 threshold, the straightforward way
    //
    for (int icand = 0; icand < _ncand; icand++) {
        if (_chi2[icand] <= chi2max) {
            std::cerr << "psf_clean(" << prof_accuracy << "): dropped candidate " 
                      << icand << "/" << _ncand << " (chi2=" << _chi2[icand] << ")" << std::endl;
        }
        else
            ret->add(_cs, icand);
    }
#else
    //
    // Drop candidates which exceed the chi^2 threshold, with psfex ordering convention
    //
    std::vector<int> ix_map(_ncand);
    n = _ncand;

    for (int i = 0; i < n; i++)
        ix_map[i] = i;

    for (int i = 0; i < n; ) {
        if (_chi2[ix_map[i]] > chi2max) {
            std::cerr << "psf_clean(" << prof_accuracy << "): dropped candidate " 
                      << ix_map[i] << "/" << _ncand << " (chi2=" << _chi2[ix_map[i]] << ")" << std::endl;
            std::swap(ix_map[i], ix_map[n-1]);
            n--;
        }
        else {
            ret->add(_cs, ix_map[i]);
            i++;
        }
    }
#endif

    return ret;
}


void PolypixPsf::psf_clip()
{
    // for debug output, see below
    std::vector<double> comp0 = _comp;

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


void PolypixPsf::psf_build(double *loc, const double *pos) const
{
    int npix = _psf_nx * _psf_ny;
    memset(loc, 0, npix * sizeof(double));

    std::vector<double> basis(_ncoeffs);
    this->poly_func(&basis[0], pos);

    for (int ic = 0; ic < _ncoeffs; ic++)
        for (int ipix = 0; ipix < npix; ipix++)
            loc[ipix] += basis[ic] * _comp[ic*npix+ipix];
}


void PolypixPsf::psf_makeresi(double prof_accuracy)
{
    const bool accuflag = (prof_accuracy > 1.0/BIG);
    const double vigstep = 1.0 / _psfstep;

    std::vector<double> pos(2);
    std::vector<double> basis(_ncoeffs);
    std::vector<double> loc(_psf_nx * _psf_ny);
    std::vector<double> dresi(_nx*_ny, 0.0);

    memset(&_vigresi[0], 0, _vigresi.size() * sizeof(double));
    memset(&_vigchi[0], 0, _vigchi.size() * sizeof(double));

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        pos[0] = (_current_xy[2*icand] - _contextoffset[0]) / _contextscale[0];
        pos[1] = (_current_xy[2*icand+1] - _contextoffset[1]) / _contextscale[1];

        this->poly_func(&basis[0], &pos[0]);
        this->psf_build(&loc[0], &pos[0]);

        // temporarily set vigresi = psf in vignette coordinates (not postmultiplied by flux)
        my_resample(&loc[0], _psf_nx, _psf_ny,
                    &_vigresi[icand*_nx*_ny], _nx, _ny,
                    -dx*vigstep, -dy*vigstep, vigstep, 1.0);

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


void PolypixPsf::poly_func(double *basis, const double *pos) const
{
    double x = pos[0];
    double y = pos[1];

    std::vector<double> xpow(_spatialOrder+1);
    std::vector<double> ypow(_spatialOrder+1);

    xpow[0] = ypow[0] = 1.0;
    for (int i = 0; i < _spatialOrder; i++) {
        xpow[i+1] = xpow[i] * x;
        ypow[i+1] = ypow[i] * y;
    }

    int outpos = 0;
    for (int iy = 0; iy <= _spatialOrder; iy++)
        for (int ix = 0; ix <= _spatialOrder-iy; ix++)
            basis[outpos++] = xpow[ix] * ypow[iy];
}


// FIXME call poly_func() in loop?
void PolypixPsf::poly_eval_basis_functions(double *basis, const double *pos) const
{
    std::vector<double> xpow(_spatialOrder+1);
    std::vector<double> ypow(_spatialOrder+1);
    int outpos = 0;

    for (int icand = 0; icand < _ncand; icand++) {        
        double x = pos[2*icand];
        double y = pos[2*icand+1];

        xpow[0] = ypow[0] = 1.0;
        for (int i = 0; i < _spatialOrder; i++) {
            xpow[i+1] = xpow[i] * x;
            ypow[i+1] = ypow[i] * y;
        }

        for (int iy = 0; iy <= _spatialOrder; iy++)
            for (int ix = 0; ix <= _spatialOrder-iy; ix++)
                basis[outpos++] = xpow[ix] * ypow[iy];
    }
}


void PolypixPsf::poly_fit(double *coeffs, const double *data, const double *weights, const double *basis, double regul) const
{
    Eigen::MatrixXd alpha(_ncoeffs, _ncoeffs);
    Eigen::VectorXd beta(_ncoeffs);

    for (int i = 0; i < _ncoeffs; i++) {
        for (int j = 0; j <= i; j++) {
            double t = 0.0;
            for (int icand = 0; icand < _ncand; icand++)
                t += weights[icand] * basis[icand*_ncoeffs+i] * basis[icand*_ncoeffs+j];
            alpha(i,j) = alpha(j,i) = t;
        }

        double t = 0.0;
        for (int icand = 0; icand < _ncand; icand++)
            t += weights[icand] * basis[icand*_ncoeffs+i] * data[icand];
        beta(i) = t;
    }

    if (regul > 0.0) {
        for (int i = 0; i < _ncoeffs; i++)
            alpha(i,i) += regul;
    }

    Eigen::VectorXd ainv_b = alpha.llt().solve(beta);

    for (int i = 0; i < _ncoeffs; i++)
        coeffs[i] = ainv_b(i);
}


// follows PsfexPsf::_doComputeImage() in meas_extensions_psfex
void PolypixPsf::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    const double vigstep = 1 / _psfstep;

    double pos[2];
    pos[0] = (x - _contextoffset[0]) / _contextscale[0];
    pos[1] = (y - _contextoffset[1]) / _contextscale[1];

    std::vector<double> fullresIm(_psf_nx * _psf_ny);      
    psf_build(&fullresIm[0], pos);

    double dx = x - x0 - 0.5*(nx_out-1);
    double dy = y - y0 - 0.5*(ny_out-1);

    my_resample(&fullresIm[0], _psf_nx, _psf_ny,
                out, nx_out, ny_out,
                -dx*vigstep, -dy*vigstep, vigstep, 1.0);

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
