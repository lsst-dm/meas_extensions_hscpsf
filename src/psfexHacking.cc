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
#include "lsst/meas/extensions/hscpsf/psfexHacking.h"

#define EPS 1.0e-4
#define BIG 1.0e30

static inline double square(double x) { return x*x; }


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}   // emacs pacifier
#endif


// private helper function for constructors
void FakePsfexPsf::_construct(int spatialOrder, double fwhm, double backnoise2, double gain)
{
    if (spatialOrder < 0)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "spatialOrder < 0 in FakePsfexPsf constructor");

    _spatialOrder = spatialOrder;
    _ncoeffs = ((spatialOrder+1) * (spatialOrder+2)) / 2;
    _fwhm = fwhm;
    _backnoise2 = backnoise2;
    _gain = gain;


    if (_ncand <= _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "too few spatial candidates in FakePsfexPsf constructor");

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
            throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "flux < 0 in FakePsfexPsf constructor");
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


FakePsfexPsf::FakePsfexPsf(CONST_PTR(HscCandidateSet) cs, int nside, int spatialOrder, double fwhm, double backnoise2, double gain)
    : HscPsfBase(cs,nside)
{ 
    this->_construct(spatialOrder, fwhm, backnoise2, gain);
    
    // FIXME revisit
    _contextoffset[0] = (_xmin + _xmax) / 2.0;
    _contextoffset[1] = (_ymin + _ymax) / 2.0;
    _contextscale[0] = (_xmax - _xmin);
    _contextscale[1] = (_ymax - _ymin);    

    if (0) {
        std::cerr << "FakePsfexPsf constructor: dumping context\n";
        for (int i = 0; i < _ncand; i++)
            std::cerr << "    FakePsfexPsf constructor: n=" << i << "  context=[ " 
                      << _initial_xy[2*i] << " " << _initial_xy[2*i+1] << " ]\n";

        std::cerr << "FakePsfexPsf constructor: cmin = [ "
                  << _xmin << " " << _ymin << " ]\n";

        std::cerr << "FakePsfexPsf constructor: cmax = [ "
                  << _xmax << " " << _ymax << " ]\n";

        std::cerr << "FakePsfexPsf constructor: setting contextscale = [ "
                  << _contextscale[0] << " " << _contextscale[1] << " ]\n";

        std::cerr << "FakePsfexPsf constructor: setting contextoffset = [ "
                  << _contextoffset[0] << " " << _contextoffset[1] << " ]\n";
    }
}


FakePsfexPsf::FakePsfexPsf(CONST_PTR(HscCandidateSet) cs, CONST_PTR(FakePsfexPsf) base)
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


void FakePsfexPsf::psf_make(double prof_accuracy)
{
    // FIXME rethink this
    const double regul = 1000.0;

    std::cerr << "psf_make(" << prof_accuracy << "): nsample=" << _ncand << std::endl;
    std::cerr << "    vigsize = (" << _nx << "," << _ny << ")\n";
    std::cerr << "    psfsize = (" << _psf_nx << "," << _psf_ny << ")\n";
    std::cerr << "    psfstep=" << _psfstep << ", fwhm=" << _fwhm << "\n";
    std::cerr << "    contextoffset = [ " << _contextoffset[0] << " " << _contextoffset[1] << " ]\n";
    std::cerr << "    contextscale = [ " << _contextscale[0] << " " << _contextscale[1] << " ]\n";

    int npix = _psf_nx * _psf_ny;
    std::vector<double> image(_ncand * npix);
    std::vector<double> weight(_ncand * npix);
    std::vector<double> pos(_ncand * 2);

    for (int icand = 0; icand < _ncand; icand++) {
        // psfex sample->dx, sample->dy
        double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
        double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

        vignet_resample_xmajor(&_im[icand*_nx*_ny], _nx, _ny,
                               &image[icand*_psf_nx*_psf_ny], _psf_nx, _psf_ny,
                               dx, dy, _psfstep, 1.0);   // FIXME note in psfex code, last arg is max(_psfstep,1.0)

        for (int ipix = 0; ipix < npix; ipix++)
            image[icand*npix + ipix] /= _norm[icand];
            
        for (int ipix = 0; ipix < npix; ipix++) {
            double val = image[icand*npix + ipix];
            double norm2 = _norm[icand] * _norm[icand];
            double profaccu2 = prof_accuracy * prof_accuracy * norm2;
            double noise2 = _backnoise2 + profaccu2 * val * val;

            if (val > 0.0 && _gain > 0.0)
                noise2 += val/_gain;

            weight[icand*npix + ipix] = norm2/noise2;
        }

        pos[2*icand] = (_current_xy[2*icand] - _contextoffset[0]) / _contextscale[0];
        pos[2*icand+1] = (_current_xy[2*icand+1] - _contextoffset[1]) / _contextscale[1];
    }

    if (0) {
        std::cerr << "psf_make(" << prof_accuracy << "): dumping vig/image/weight\n";
        for (int icand = 0; icand < _ncand; icand++) {
            // psfex sample->dx, sample->dy
            double dx = _current_xy[2*icand] - _xy0[2*icand] - 0.5*(double)(_nx-1);
            double dy = _current_xy[2*icand+1] - _xy0[2*icand+1] - 0.5*(double)(_ny-1);

#if 0            
            int vi = (1013*icand+10) % _nx;
            int vj = (1031*icand+10) % _ny;
            int ii = (1013*icand+10) % _psf_nx;
            int ij = (1031*icand+10) % _psf_ny;
#else
            int vi = 11;
            int vj = 11;
            int ii = 11;
            int ij = 11;
#endif
            
            std::cerr << "  n=" << icand 
                      << "     norm=" << _norm[icand]
                      << "     gain=" << _gain
                      << "     backnoise2=" << _backnoise2
                      << "     (dx,dy)=(" << dx << "," << dy << ")"
                      << "\n            "
                      << "     vig[" << vi << "," << vj << "]=" << _im[icand*_nx*_ny + vi*_ny + vj]
                      << "     image[" << ii << "," << ij << "]=" << image[icand*_psf_nx*_psf_ny + ii*_psf_ny + ij]
                      << "     weight[" << ii << "," << ij << "]=" << weight[icand*_psf_nx*_psf_ny + ii*_psf_ny + ij]
                      << "\n            "
                      << "     context=[ " << _current_xy[2*icand] << " " << _current_xy[2*icand+1] << " ]"
                      << "     pos=[ " << pos[2*icand] << " " << pos[2*icand+1] << " ]"
                      << std::endl;
        }
    }

    std::vector<double> basis(_ncand * _ncoeffs);
    this->poly_eval_basis_functions(&basis[0], &pos[0]);

    std::vector<double> pstack(_ncand);
    std::vector<double> wstack(_ncand);
    std::vector<double> coeff(_ncoeffs);

    for (int ipix = 0; ipix < npix; ipix++) {
        for (int icand = 0; icand < _ncand; icand++) {
            pstack[icand] = image[icand*npix + ipix];
            wstack[icand] = weight[icand*npix + ipix];
        }

        this->poly_fit(&coeff[0], &pstack[0], &wstack[0], &basis[0], regul);

        for (int icoeff = 0; icoeff < _ncoeffs; icoeff++)
            _comp[icoeff*npix + ipix] = coeff[icoeff];
    }

    if (0) {
        std::cerr << "psf_make(" << prof_accuracy << "): dumping basis_functions\n";
        for (int icand = 0; icand < _ncand; icand++) {
            std::cerr << "  psf_make(" << prof_accuracy << "): ncand=" << icand << "  basis_functions = [";
            for (int icoeff = 0; icoeff < _ncoeffs; icoeff++)
                std::cerr << " " << basis[icand*_ncoeffs + icoeff];
            std::cerr << " ]\n";
        }
    }

    if (0) {
        std::cerr << "psf_make(" << prof_accuracy << "): dumping polynomial coefficients\n";
        for (int ix = 0; ix < _psf_nx; ix++) {
            for (int iy = 0; iy < _psf_ny; iy++) {
                std::cerr << "    psf_make (ix,iy)=(" << ix << "," << iy << "): comp = [";
                for (int icoeff = 0; icoeff < _ncoeffs; icoeff++)
                    std::cerr << " " << _comp[icoeff*npix + ix*_psf_ny + iy];
                std::cerr << " ]\n";
            }
        }
    }
}


PTR(HscCandidateSet) FakePsfexPsf::psf_clean(double prof_accuracy)
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
        if (_chi2[icand] > chi2max) {
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


void FakePsfexPsf::psf_clip()
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
            
    if (0) {
        for (int ix = 0; ix < _psf_nx; ix++) {
            for (int iy = 0; iy < _psf_ny; iy++) {
                std::cerr << "    psf_clip (ix,iy)=(" << ix << "," << iy << "):  [";
                for (int ic = 0; ic < _ncoeffs; ic++)
                    std::cerr << " " << comp0[ic*_psf_nx*_psf_ny + ix*_psf_ny + iy];
                std::cerr << " ] -> [";
                for (int ic = 0; ic < _ncoeffs; ic++)
                    std::cerr << " " << _comp[ic*_psf_nx*_psf_ny + ix*_psf_ny + iy];
                std::cerr << "]\n";
            }
        }
    }
}


void FakePsfexPsf::psf_build(double *loc, const double *pos) const
{
    int npix = _psf_nx * _psf_ny;
    memset(loc, 0, npix * sizeof(double));

    std::vector<double> basis(_ncoeffs);
    this->poly_func(&basis[0], pos);

    for (int ic = 0; ic < _ncoeffs; ic++)
        for (int ipix = 0; ipix < npix; ipix++)
            loc[ipix] += basis[ic] * _comp[ic*npix+ipix];
}


void FakePsfexPsf::psf_makeresi(double prof_accuracy)
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
        vignet_resample_xmajor(&loc[0], _psf_nx, _psf_ny,
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

    if (0) {
        for (int icand = 0; icand < _ncand; icand++)
            std::cerr << "    psf_makeresi(" << prof_accuracy << "): icand=" << icand << ": chi2=" << _chi2[icand] << std::endl;
    }
}


// -------------------------------------------------------------------------------------------------


void FakePsfexPsf::poly_func(double *basis, const double *pos) const
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
void FakePsfexPsf::poly_eval_basis_functions(double *basis, const double *pos) const
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


void FakePsfexPsf::poly_fit(double *coeffs, const double *data, const double *weights, const double *basis, double regul) const
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
void FakePsfexPsf::eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const
{
    const double vigstep = 1 / _psfstep;

    double pos[2];
    pos[0] = (x - _contextoffset[0]) / _contextscale[0];
    pos[1] = (y - _contextoffset[1]) / _contextscale[1];

    std::vector<double> fullresIm(_psf_nx * _psf_ny);      
    psf_build(&fullresIm[0], pos);

    double dx = x - x0 - 0.5*(nx_out-1);
    double dy = y - y0 - 0.5*(ny_out-1);

    vignet_resample_xmajor(&fullresIm[0], _psf_nx, _psf_ny,
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


// -------------------------------------------------------------------------------------------------


#define	INTERPW		6	/* Interpolation function range */
#define	INTERPFAC	3.0	/* Interpolation envelope factor */

#define	INTERPF(x)	(x<1e-5 && x>-1e-5? 1.0 \
			:(x>INTERPFAC?0.0:(x<-INTERPFAC?0.0 \
			:sin(M_PI*x)*sin(M_PI/INTERPFAC*x)/(M_PI*M_PI/INTERPFAC*x*x))))


// cut-and-paste from psfex with minor changes
void vignet_resample_ymajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                            double dx, double dy, double step2, double stepi)
{
    std::vector<double> mask, pix12;
    double *maskt, *pixout, *pixout0;
    const double *pixin, *pixin0;

    std::vector<int> nmask, start;
    int *nmaskt, *startt;

    double  mx1, mx2, my1, my2, xs1,ys1, x1, y1, x, y, dxm, dym, val, dstepi, norm;
    int	    i, j, k, n, t;
    int     ixs2, iys2, ix2, iy2, dix2, diy2, nx2, ny2, iys1a, ny1, hmw, hmh;
    int     ix, iy, ix1, iy1, interpw, interph;

    if (stepi <= 0.0)
        stepi = 1.0;
    dstepi = 1.0/stepi;


    mx1 = (double)(w1/2);		/* Im1 center x-coord*/
    mx2 = (double)(w2/2);		/* Im2 center x-coord*/
    xs1 = mx1 + dx - mx2*step2;	/* Im1 start x-coord */
    if ((int)xs1 >= w1)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "xs1 too large");
    ixs2 = 0;			/* Int part of Im2 start x-coord */
    if (xs1<0.0)
        {
            dix2 = (int)(1-xs1/step2);
            /*-- Simply leave here if the images do not overlap in x */
            if (dix2 >= w2)
                throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "dix2 too large");
            ixs2 += dix2;
            xs1 += dix2*step2;
        }
    nx2 = (int)((w1-1-xs1)/step2+1);/* nb of interpolated Im2 pixels along x */
    if (nx2>(ix2=w2-ixs2))
        nx2 = ix2;
    if (nx2<=0)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "nx2 <= 0");


    my1 = (double)(h1/2);		/* Im1 center y-coord */
    my2 = (double)(h2/2);		/* Im2 center y-coord */
    ys1 = my1 + dy - my2*step2;	/* Im1 start y-coord */
    if ((int)ys1 >= h1)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "ys1 too large");
    iys2 = 0;			/* Int part of Im2 start y-coord */
    if (ys1<0.0)
        {
            diy2 = (int)(1-ys1/step2);
            /*-- Simply leave here if the images do not overlap in y */
            if (diy2 >= h2)
                throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "diy2 too large");
            iys2 += diy2;
            ys1 += diy2*step2;
        }
    ny2 = (int)((h1-1-ys1)/step2+1);/* nb of interpolated Im2 pixels along y */
    if (ny2>(iy2=h2-iys2))
        ny2 = iy2;
    if (ny2<=0)
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "ny2 <= 0");

    /* Set the yrange for the x-resampling with some margin for interpolation */
    iys1a = (int)ys1;		/* Int part of Im1 start y-coord with margin */
    hmh = (int)((INTERPW/2)/dstepi) + 2;	/* Interpolant start */
    interph = 2*hmh;
    hmw = (int)((INTERPW/2)/dstepi) + 2;
    interpw =  2*hmw;
    if (iys1a<0 || ((iys1a -= hmh)< 0))
        iys1a = 0;
    ny1 = (int)(ys1+ny2*step2)+interpw-hmh;	/* Interpolated Im1 y size */
    if (ny1>h1)					/* with margin */
        ny1 = h1;
    /* Express everything relative to the effective Im1 start (with margin) */
    ny1 -= iys1a;
    ys1 -= (double)iys1a;
    
    /* Allocate interpolant stuff for the x direction */
    mask = std::vector<double>(nx2 * interpw, 0.0);
    nmask = std::vector<int>(nx2, 0);
    start = std::vector<int>(nx2, 0);

    /* Compute the local interpolant and data starting points in x */
    x1 = xs1;
    maskt = &mask[0];
    nmaskt = &nmask[0];
    startt = &start[0];
    for (j=nx2; j--; x1+=step2)
        {
            ix = (ix1=(int)x1) - hmw;
            dxm = (ix1 - x1 - hmw)*dstepi;/* starting point in the interp. func */
            if (ix < 0)
                {
                    n = interpw+ix;
                    dxm -= (double)ix*dstepi;
                    ix = 0;
                }
            else
                n = interpw;
            if (n>(t=w1-ix))
                n=t;
            *(startt++) = ix;
            *(nmaskt++) = n;
            norm = 0.0;
            for (x=dxm, i=n; i--; x+=dstepi)
                norm +=( *(maskt++) = INTERPF(x));
            norm = norm>0.0? 1.0/norm : dstepi;
            maskt -= n;
            for (i=n; i--;)
                *(maskt++) *= norm;
        }
    
    /* Intermediary frame-buffer */
    pix12 = std::vector<double>(nx2*ny1, 0.0);
    
    /* Make the interpolation in x (this includes transposition) */
    pixin0 = pix1+iys1a*w1;
    pixout0 = &pix12[0];
    for (k=ny1; k--; pixin0+=w1, pixout0++)
        {
            maskt = &mask[0];
            nmaskt = &nmask[0];
            startt = &start[0];
            pixout = pixout0;
            for (j=nx2; j--; pixout+=ny1)
                {
                    pixin = pixin0+*(startt++);
                    val = 0.0; 
                    for (i=*(nmaskt++); i--;)
                        val += *(maskt++)*(double)*(pixin++);
                    *pixout = (double)val;
                }
        }
    
    /* Reallocate interpolant stuff for the y direction */
    mask = std::vector<double>(ny2 * interph, 0.0);  /* Interpolation masks */
    nmask = std::vector<int>(ny2, 0);                /* Interpolation mask sizes */
    start = std::vector<int>(ny2, 0);                /* Int part of Im1 conv starts */
    
    /* Compute the local interpolant and data starting points in y */
    y1 = ys1;
    maskt = &mask[0];
    nmaskt = &nmask[0];
    startt = &start[0];
    for (j=ny2; j--; y1+=step2)
        {
            iy = (iy1=(int)y1) - hmh;
            dym = (iy1 - y1 - hmh)*dstepi;/* starting point in the interp. func */
            if (iy < 0)
                {
                    n = interph+iy;
                    dym -= (double)iy*dstepi;
                    iy = 0;
                }
            else
                n = interph;
            if (n>(t=ny1-iy))
                n=t;
            *(startt++) = iy;
            *(nmaskt++) = n;
            norm = 0.0;
            for (y=dym, i=n; i--; y+=dstepi)
                norm += (*(maskt++) = INTERPF(y));
            norm = norm>0.0? 1.0/norm : dstepi;
            maskt -= n;
            for (i=n; i--;)
                *(maskt++) *= norm;
        }
    
#if 0
    // KMS the statpix2 thing looks like a recipe for memory corruption so I removed it
    static double *statpix2;
    /* Initialize destination buffer to zero if pix2 != NULL */
    if (!pix2)
        pix2 = statpix2;
    else
        {
            memset(pix2, 0, (size_t)(w2*h2)*sizeof(double));
            statpix2 = pix2;
        }
#else
    if (!pix2)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "pix2 == NULL");
    memset(pix2, 0, (size_t)(w2*h2)*sizeof(double));    
#endif

    /* Make the interpolation in y  and transpose once again */
    pixin0 = &pix12[0];
    pixout0 = pix2+ixs2+iys2*w2;
    for (k=nx2; k--; pixin0+=ny1, pixout0++)
        {
            maskt = &mask[0];
            nmaskt = &nmask[0];
            startt = &start[0];
            pixout = pixout0;
            for (j=ny2; j--; pixout+=w2)
                {
                    pixin = pixin0+*(startt++);
                    val = 0.0; 
                    for (i=*(nmaskt++); i--;)
                        val += *(maskt++)*(double)*(pixin++);
                    *pixout = (double)val;
                }
        }
}


void vignet_resample_xmajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                            double dx, double dy, double step2, double stepi)
{
    std::vector<double> pix1_ymajor(w1*h1);
    std::vector<double> pix2_ymajor(w2*h2);
    
    for (int i = 0; i < w1; i++)
        for (int j = 0; j < h1; j++)
            pix1_ymajor[j*w1+i] = pix1[i*h1+j];

    vignet_resample_ymajor(&pix1_ymajor[0], w1, h1, &pix2_ymajor[0], w2, h2, dx, dy, step2, stepi);

    for (int i = 0; i < w2; i++)
        for (int j = 0; j < h2; j++)
            pix2[i*h2+j] = pix2_ymajor[j*w2+i];
}


// -------------------------------------------------------------------------------------------------


#define MEDIAN_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double fast_median(double *arr, int n)
{
    double  *alow, *ahigh, *amedian, *amiddle, *all, *ahh, val, valmax, valmax2;
    int      i, nless;

    if (!n)
        return 0.0;
    else if (n==1)
        return *arr;
    else if (n==2)
        return 0.5*(*arr+*(arr+1));
    
    alow = arr;
    ahigh = arr + n - 1;
    amedian = arr + n/2;
    while (ahigh > (all=alow + 1)) {
        /*-- Find median of low, middle and high items; swap into position low */
        amiddle = alow + (ahigh-alow)/2;
        if (*amiddle > *ahigh)
            MEDIAN_SWAP(*amiddle, *ahigh);
        if (*alow > *ahigh)
            MEDIAN_SWAP(*alow, *ahigh);
        if (*amiddle > *alow)
            MEDIAN_SWAP(*amiddle, *alow);
        
        /*-- Swap low item (now in position middle) into position (low+1) */
        MEDIAN_SWAP(*amiddle, *all);

        /*-- Nibble from each end towards middle, swapping items when stuck */
        ahh = ahigh;
        for (;;) {
            while (*alow > *(++all))
                ;
            while (*(--ahh) > *alow)
                ;
            if (ahh < all)
                break;
            MEDIAN_SWAP(*all, *ahh);
        }

        /*-- Swap middle item (in position low) back into correct position */
        MEDIAN_SWAP(*alow, *ahh) ;

        /*-- Re-set active partition */
        if (ahh <= amedian)
            alow = all;
        if (ahh >= amedian)
            ahigh = ahh - 1;
    }

    /* One or two elements left */
    if (ahigh == all && *alow > *ahigh)
        MEDIAN_SWAP(*alow, *ahigh);

    if (n&1)
        /*-- Odd case */
        return *amedian;
    else {
        /*-- Even case */
        valmax2 = *amedian;
        valmax = -BIG;
        alow = arr;
        nless = 0;
        for (i=n/2;i--;) {
            if ((val=*(alow++))<valmax2) {
                nless++;
                if (val > valmax)
                    valmax = val;
            }
        }
        return nless<n/2? *amedian : (*amedian+valmax)/2.0;
    }   
}



}}}}   // namespace lsst::meas::extensions::hscpsf

