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

#ifndef LSST_MEAS_EXTENSIONS_HSCPSF_PSFEX_HACKING_H
#define LSST_MEAS_EXTENSIONS_HSCPSF_PSFEX_HACKING_H 1

#include "lsst/pex/exceptions.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


class FakePsfexPsf : public HscPsfBase
{
public:
    FakePsfexPsf(CONST_PTR(HscCandidateSet) cs, int nside, int spatialOrder, double fwhm, double backnoise2, double gain);

    void psf_make(double prof_accuracy);
    void psf_makeresi(double prof_accuracy);
    void psf_clip();

    PTR(HscCandidateSet) psf_clean(double prof_accuracy);

    //
    // @out = output array of shape (_psf_nx, _psf_ny)
    // @pos = input_array of shape (2,)
    //
    void psf_build(double *loc, const double *pos);

    //
    // @basis = output array of shape (ncoeffs,)
    // @pos = input array of shape (2,)
    //
    void poly_func(double *basis, const double *pos);

    // 
    // @basis = output array of shape (ncand, ncoeffs)
    // @pos = input array of shape (ncand, 2)
    //
    void poly_eval_basis_functions(double *basis, const double *pos);

    //
    // @coeffs = output array of shape (ncoeffs,)
    // @data = input array of shape (ncand,)
    // @weights = input array of shape (ncand,)
    // @basis = input array of shape (ncand, ncoeffs), probably computed using poly_eval_basis_functions()
    //
    void poly_fit(double *coeffs, const double *data, const double *weights, const double *basis, double regul);

protected:
    int _spatialOrder;
    int _ncoeffs;

    double _fwhm;
    double _backnoise2;
    double _gain;
    
    int _psf_nx;             // psfex psf->size[0]
    int _psf_ny;             // psfex psf->size[1]
    double _psfstep;         // FIXME look at psfex source and figure out the difference between pixstep and psfstep

    std::vector<double> _norm;           // shape (_ncand); psfex sample->norm
    std::vector<double> _contextoffset;
    std::vector<double> _contextscale;

    std::vector<double> _vigweight;      // shape (_ncand, _nx, _ny)
    std::vector<double> _comp;           // shape (_ncoeffs, _psf_nx, psf_ny)

    // psf_makeresi() sets this data
    std::vector<double> _vigresi;        // shape (_ncand, _nx, _ny)
    std::vector<double> _vigchi;         // shape (_ncand, _nx, _ny)    per-pixel contribution to chi^2
    std::vector<double> _chi2;           // shape (_ncand,)             reduced chi^2
};


// -------------------------------------------------------------------------------------------------


extern double fast_median(double *arr, int n);

// you probably want this one...
extern void vignet_resample_xmajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                                   double dx, double dy, double step2, double stepi);

// ...not this one
extern void vignet_resample_ymajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                                   double dx, double dy, double step2, double stepi);


}}}}   // namespace lsst::meas::extensions::hscpsf


#endif  // LSST_MEAS_EXTENSIONS_HSCPSF_PSFEX_HACKING_H
