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
 

#include "lsst/meas/extensions/hscpsf/hscPsf.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {


template<typename T> static inline double vecdot(ssize_t len, const T *v1, const T *v2)
{
    double ret = 0.0;
    for (ssize_t i = 0; i < len; i++)
	ret += v1[i] * v2[i];
    return ret;
}

// vecma = "vector multiply-add"
template<typename T> static inline void vecma(ssize_t len, T *dst, const T *src, double alpha)
{
    for (ssize_t i = 0; i < len; i++)
	dst[i] += alpha*src[i];
}

template<typename T> static inline void vecscale(ssize_t len, T *vec, double alpha)
{
    for (ssize_t i = 0; i < len; i++)
	vec[i] *= alpha;
}


void HscSpatialModelBase::normalizePcaImages(int npca, int npix, int ncand, double *pca, double *ampl, double *sm) const
{
    const int ncoeffs = this->getNcoeffs();

    for (int i = 1; i < npca; i++) {
	for (int j = 1; j < i; j++) {
	    double t = vecdot(npix, &pca[i*npix], &pca[j*npix]);
	    vecma(npix, &pca[i*npix], &pca[j*npix], -t);
	    vecma(ncoeffs, &sm[(j-1)*ncoeffs], &sm[(i-1)*ncoeffs], t);
	}

	double t = vecdot(npix, &pca[i*npix], &pca[i*npix]);
	vecscale(npix, &pca[i*npix], 1.0/sqrt(t));
	vecscale(ncoeffs, &sm[(i-1)*ncoeffs], sqrt(t));

	t = vecdot(npix, &pca[0], &pca[i*npix]);
	vecma(npix, &pca[0], &pca[i*npix], -t);
	sm[(i-1)*ncoeffs] += t;   // shift by constant
    }

    double t = vecdot(npix, &pca[0], &pca[0]);
    vecscale(npix, &pca[0], 1.0/sqrt(t));
    vecscale(ncand, &ampl[0], sqrt(t));
    vecscale((npca-1)*ncoeffs, &sm[0], 1.0/sqrt(t));

    // FIXME ensure pca0 positive?
}


ndarray::Array<double,2,2> HscSpatialModelBase::eval(const ndarray::Array<const double,2,2> &xy, const ndarray::Array<const double,2,2> &fcoeffs) const
{
    if (xy.getSize<1>() != 2)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "xy argument to HscSpatialModelBase::eval() must be an array of shape (N,2)");
    if (fcoeffs.getSize<1>() != this->getNcoeffs())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "fcoeffs argument to HscSpatialModelBase::eval() must be an array of shape (N,Ncoeffs)");

    int ncand = xy.getSize<0>();
    int nfunc = fcoeffs.getSize<0>();

    ndarray::Array<double,2,2> out = ndarray::allocate(ncand, nfunc);
    this->eval(out.getData(), ncand, xy.getData(), nfunc, fcoeffs.getData());

    return out;
}


ndarray::Array<double,1,1> HscSpatialModelBase::optimize(const ndarray::Array<const double,2,2> &xy, const ndarray::Array<const double,1,1> &a, const ndarray::Array<const double,1,1> &b, double regul) const
{
    if (xy.getSize<1>() != 2)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "xy argument to HscSpatialModelBase::optimize() must be an array of shape (N,2)");
    if (a.getSize<0>() != xy.getSize<0>())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "xy,a arrays in HscSpatialModelBase::optimize() must have the same length");
    if (b.getSize<0>() != xy.getSize<0>())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "xy,b arrays in HscSpatialModelBase::optimize() must have the same length");

    int ncand = xy.getSize<0>();
    ndarray::Array<double,1,1> out = ndarray::allocate(this->getNcoeffs());

    this->optimize(out.getData(), ncand, xy.getData(), a.getData(), b.getData(), regul);
    return out;
}


void HscSpatialModelBase::normalizePcaImages(const ndarray::Array<double,2,2> &pcas, const ndarray::Array<double,1,1> &ampl, const ndarray::Array<double,2,2> &sm) const
{
    int npca = pcas.getSize<0>();
    int npix = pcas.getSize<1>();
    int ncand = ampl.getSize<0>();

    if (sm.getSize<0>() != (npca-1))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "sm array in HscSpatialModelBase::normalizePcaImages() must have leading dimension =(npca-1)");
    if (sm.getSize<1>() != this->getNcoeffs())
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "sm array in HscSpatialModelBase::normalizePcaImages() must have second dimension =ncoeffs");
    
    this->normalizePcaImages(npca, npix, ncand, pcas.getData(), ampl.getData(), sm.getData());
}


}}}}   // namespace lsst::meas::extensions::hscpsf
