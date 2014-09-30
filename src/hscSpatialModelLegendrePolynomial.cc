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
//   P_0(x)*P_0(y)  
//   P_1(x)*P_0(y)  P_0(x)*P_1(x)
//       ...
//   P_N(x)*P_0(y)  P_{N-1}(x)*P_1(y)  ...  P_0(x)*P_N(y)
//
HscSpatialModelLegendrePolynomial::HscSpatialModelLegendrePolynomial(int order, double xmin, double xmax, double ymin, double ymax)
    : _order(order), _xmin(xmin), _xmax(xmax), _ymin(ymin), _ymax(ymax)
{
    if ((order < 0) || (xmin >= xmax) || (ymin >= ymax))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters to HscSpatialModelLegendrePolynomial constructor");

    _ncoeffs = ((order+1) * (order+2)) / 2;    // integer divide
}


// virtual
int HscSpatialModelLegendrePolynomial::getNcoeffs() const 
{ 
    return _ncoeffs; 
}


// virtual
void HscSpatialModelLegendrePolynomial::eval(double *out, int ncand, const double *xy, int nfunc, const double *fcoeffs) const
{
    if ((ncand <= 0) || (nfunc <= 0) || (xy == NULL) || (out == NULL) || (fcoeffs == NULL))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid parameters to HscSpatialModelLegendrePolynomial::eval()");

    std::vector<double> pl((_order+1) * ncand * 2);    // shape (order+1,ncand,2)
    _eval_pl_scaled(&pl[0], ncand, xy);
    
    Eigen::VectorXd px(_order+1);
    Eigen::VectorXd py(_order+1);
    Eigen::MatrixXd c(_order+1, _order+1);

    for (int icand = 0; icand < ncand; icand++) {
        c.setZero();

	for (int i = 0; i <= _order; i++) {
	    px(i) = pl[(2*ncand)*i + 2*icand];
	    py(i) = pl[(2*ncand)*i + 2*icand + 1];
	}

	int icoeff = 0;
	for (int ifunc = 0; ifunc < nfunc; ifunc++) {
	    for (int n = 0; n <= _order; n++)
		for (int i = 0; i <= n; i++)
		    c(n-i,i) = fcoeffs[icoeff++];

	    out[icand*nfunc + ifunc] = px.dot(c * py);
	}

	assert(icoeff == nfunc * _ncoeffs);   // just checking
    }
}


// virtual 
void HscSpatialModelLegendrePolynomial::optimize(double *out, int ncand, const double *xy, const double *a, const double *b, double regul) const
{
    if ((a == NULL) || (b == NULL) || (xy == NULL) || (out == NULL))
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "null pointer passed to HscSpatialModelLegendrePolynomial::optimize()");
    if (ncand < _ncoeffs)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "number of candidates too small in HscSpatialModelLegendrePolynomial::optimize()");

    Eigen::VectorXd bmodel(_ncoeffs);
    Eigen::MatrixXd tmodel(ncand, _ncoeffs);

    std::vector<double> pl((_order+1) * ncand * 2);    // shape (order+1,ncand,2)
    _eval_pl_scaled(&pl[0], ncand, xy);

    std::vector<double> plx(_order+1);
    std::vector<double> ply(_order+1);

    bmodel.setZero();

    for (int icand = 0; icand < ncand; icand++) {
        if (a[icand] <= 0.0)
            throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "a[i] <= 0 in HscSpatialModelLegendrePolynomial::optimize()");

        double sqrta = std::sqrt(a[icand]);

	for (int i = 0; i <= _order; i++) {
	    plx[i] = pl[(2*ncand)*i + 2*icand];
	    ply[i] = pl[(2*ncand)*i + 2*icand + 1];
	}

	int icoeff = 0;
	for (int n = 0; n <= _order; n++) {
	    for (int i = 0; i <= n; i++) {
		bmodel(icoeff) += b[icand] * plx[n-i] * ply[i];
                tmodel(icand, icoeff) = sqrta * plx[n-i] * ply[i];
                icoeff++;
            }
        }

	assert(icoeff == _ncoeffs);   // just checking!
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
// @out = array of shape (order+1, ncand)
// @t = 1D array of length ncand
//
void HscSpatialModelLegendrePolynomial::_eval_pl_unscaled(double *out, int ncand, const double *t) const
{
    for (int i = 0; i < ncand; i++)
	out[i] = 1.0;

    if (_order == 0)
	return;

    for (int i = 0; i < ncand; i++)
	out[ncand+i] = t[i];

    for (int l = 2; l <= _order; l++) {
	double alpha = (double)(2*l-1) / (double)l;
	double beta = (double)(l-1) / (double)l;

	for (int i = 0; i < ncand; i++)
	    out[l*ncand+i] = alpha * t[i] * out[(l-1)*ncand+i] - beta * out[(l-2)*ncand+i];
    }
}


//
// @out = array of shape (order+1, ncand, 2)
// @xy = array of length shape (ncand, 2)
//
// FIXME: this routine should include range checking
//
void HscSpatialModelLegendrePolynomial::_eval_pl_scaled(double *out, int ncand, const double *xy) const
{
    assert((out != NULL) && (xy != NULL) && (ncand > 0));

    std::vector<double> t(2 * ncand);
    for (int i = 0; i < ncand; i++) {
	t[2*i] = 2.0 * (xy[2*i]-_xmin) / (_xmax-_xmin) - 1.0;
	t[2*i+1] = 2.0 * (xy[2*i+1]-_ymin) / (_ymax-_ymin) - 1.0;
    }

    _eval_pl_unscaled(out, 2*ncand, &t[0]);
}


}}}}   // namespace lsst::meas::extensions::hscpsf
