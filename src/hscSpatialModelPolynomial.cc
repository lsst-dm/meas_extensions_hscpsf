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


//
// Note: ordering of coefficients is
//   1   x   x^2  ...  x^{N-1} x^N
//   y  xy   x^2y ... x^{N-1}y
//      ...
//   x^N
//
HscSpatialModelPolynomial::HscSpatialModelPolynomial(int order, double xmin, double xmax, double ymin, double ymax)
    : _order(order), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax)
{
    if ((order < 0) || (xmin >= xmax) || (ymin >= ymax))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters to HscSpatialModelPolynomial constructor");

    _ncoeffs = ((order+1) * (order+2)) / 2;    // integer divide
}


// virtual
int HscSpatialModelPolynomial::getNcoeffs() const 
{ 
    return _ncoeffs; 
}


//
// @out = array of shape (ncand,nfunc)
// @xy = array of shape (ncand,2)
// @fcoeffs = array of shape (nfunc,ncoeffs)
//
void HscSpatialModelPolynomial::eval(double *out, int ncand, const double *xy, int nfunc, const double *fcoeffs) const
{
    if ((ncand <= 0) || (nfunc <= 0) || (xy == NULL) || (out == NULL) || (fcoeffs == NULL))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters to HscSpatialModelPolynomial::eval()");

    std::vector<double> xypow((_order+1) * ncand * 2);    // shape (order+1,ncand,2)
    this->_eval_xypow_scaled(&xypow[0], ncand, xy);

    // FIXME (low priority) this loop could be rewritten as an Eigen matrix multiply for slightly better speed
    for (int icand = 0; icand < ncand; icand++) {
        for (int ifunc = 0; ifunc < nfunc; ifunc++) {
            double p = 0.0;

            int ic = 0;
            for (int iy = 0; iy <= _order; iy++) {
                for (int ix = 0; ix <= _order - iy; ix++) {
                    p += fcoeffs[ifunc*_ncoeffs + ic] * xypow[ix*(2*ncand)+2*icand] * xypow[iy*(2*ncand)+2*icand+1];
                    ic++;
                }
            }

            assert(ic == _ncoeffs);   // just checking!
            out[icand*nfunc + ifunc] = p;
        }
    }
}

                
//
// @out = 1D array of length ncoeffs
// @xy = array of shape (ncand,2)
// @a = 1D array of length ncand
// @b = 1D array of length ncand
//
void HscSpatialModelPolynomial::optimize(double *out, int ncand, const double *xy, const double *a, const double *b, double regul) const
{
    if ((a == NULL) || (b == NULL) || (xy == NULL) || (out == NULL))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "null pointer passed to HscSpatialModelPolynomial::optimize()");
    if (ncand < _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "number of candidates too small in HscSpatialModelPolynomial::optimize()");

    Eigen::VectorXd bmodel(_ncoeffs);
    Eigen::MatrixXd tmodel(ncand, _ncoeffs);

    std::vector<double> xypow((_order+1) * ncand * 2);    // shape (order+1,ncand,2)
    this->_eval_xypow_scaled(&xypow[0], ncand, xy);

    bmodel.setZero();

    for (int icand = 0; icand < ncand; icand++) {
        if (a[icand] <= 0.0)
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "a[i] <= 0 in HscSpatialModelPolynomial::optimize()");

        double sqrta = std::sqrt(a[icand]);

        int ic = 0;
        for (int iy = 0; iy <= _order; iy++) {
            for (int ix = 0; ix <= _order - iy; ix++) {
                double t = xypow[ix*(2*ncand)+2*icand] * xypow[iy*(2*ncand)+2*icand+1];
                tmodel(icand,ic) = sqrta * t;
                bmodel(ic) += b[icand] * t;
                ic++;
            }
        }

        assert(ic == _ncoeffs);   // just checking!
    }

    Eigen::MatrixXd amodel = tmodel.transpose() * tmodel;

    if (regul > 0.0) {
        for (int i = 0; i < _ncoeffs; i++)
            amodel(i,i) += regul;
    }

    Eigen::VectorXd ainv_b = amodel.llt().solve(bmodel);

    for (int i = 0; i < _ncoeffs; i++)
        out[i] = ainv_b(i);
}


//
// @out = array of shape (order+1, ncand, 2)
// @xy = array of length shape (ncand, 2)
//
// FIXME: this routine should include range checking
//
void HscSpatialModelPolynomial::_eval_xypow_scaled(double *out, int ncand, const double *xy) const
{
    assert((out != NULL) && (xy != NULL) && (ncand > 0));

    std::vector<double> t(2 * ncand);
    for (int i = 0; i < ncand; i++) {
	t[2*i] = (xy[2*i]-_xmin) / (_xmax-_xmin) - 0.5;
	t[2*i+1] = (xy[2*i+1]-_ymin) / (_ymax-_ymin) - 0.5;
    }

    for (int i = 0; i < 2*ncand; i++)
        out[i] = 1.0;

    for (int ic = 1; ic <= _order; ic++)
        for (int i = 0; i < 2*ncand; i++)
            out[ic*(2*ncand) + i] = t[i] * out[(ic-1)*(2*ncand) + i];
}


}}}}   // namespace lsst::meas::extensions::hscpsf
