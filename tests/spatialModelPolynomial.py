#!/usr/bin/env python

# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
# 
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the LSST License Statement and 
# the GNU General Public License along with this program.  If not, 
# see <http://www.lsstcorp.org/LegalNotices/>.
#

"""
Run with:
   spatialModel.py
or
   python
   >>> import spatialModel; spatialModel.run()
"""

import sys
import unittest

import numpy as np
import numpy.random as rand

import lsst.utils.tests as tests
import lsst.meas.extensions.hscpsf.hscpsfLib as hscpsfLib


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class SpatialModelTestBase(unittest.TestCase):
    """
    Base class for testing spatial models.

    Subclass must define
       self.construct_sm()      -> initializes self.sm to a random(?) spatial model
       self.eval2(xy, fcoeffs)  -> same semantics as self.sm.eval(), for comparison
    """

    def setUp(self):
        pass

    def tearDown(self):
        pass


    def _eval_chi2(self, a, b, coeffs):
        fcoeffs = np.zeros((self.nfunc, self.ncoeffs))
        fcoeffs[0,:] = coeffs

        u = self.sm.eval(self.xy, fcoeffs)[:,0]
        return np.sum(0.5*a*u*u - b*u)


    def _eval_chi2_gradient(self, a, b, coeffs, epsilon=1.0e-2):
        ret = np.zeros(self.ncoeffs)

        for icoeffs in xrange(self.ncoeffs):
            dcoeffs = np.zeros(self.ncoeffs)
            dcoeffs[icoeffs] = epsilon
            ret[icoeffs] = (self._eval_chi2(a,b,coeffs+dcoeffs) - self._eval_chi2(a,b,coeffs-dcoeffs)) / (2*epsilon)

        return ret


    def _eval_pcas(self, pcas, ampl, sm):
        """Returns array of shape (ncand, npix)."""

        assert pcas.shape == (self.npca, self.npix)
        assert ampl.shape == (self.ncand,)
        assert sm.shape == (self.npca-1, self.ncoeffs)

        x = np.zeros((self.ncand, self.npca))
        x[:,0] = ampl[:]

        for ipca in xrange(1, self.npca):
            fcoeffs = np.zeros((self.nfunc, self.ncoeffs))
            fcoeffs[0,:] = sm[ipca-1,:]
            u = self.sm.eval(self.xy, fcoeffs)
            x[:,ipca] = ampl * u[:,0]

        return np.dot(x, pcas)
        

    def _test_eval(self):
        fcoeffs = rand.standard_normal(size=(self.nfunc, self.ncoeffs))
        eval1 = self.sm.eval(self.xy, fcoeffs)
        eval2 = self.eval2(self.xy, fcoeffs)

        self.assertTrue(eval1.shape == (self.ncand, self.nfunc))
        self.assertTrue(eval2.shape == (self.ncand, self.nfunc))
        self.assertLess(np.max(np.abs(eval1-eval2)), 1.0e-12)


    def _test_optimize(self):
        a = rand.uniform(low=1.0, high=2.0, size=self.ncand)
        b = rand.standard_normal(size=self.ncand)

        opt_coeffs = self.sm.optimize(self.xy, a, b, 0.0)

        grad = self._eval_chi2_gradient(a, b, opt_coeffs)
        self.assertLess(np.sum(grad**2)**0.5, 1.0e-10)


    def _test_normalize(self):
        pca = rand.standard_normal(size=(self.npca,self.npix))
        ampl = rand.uniform(size=self.ncand)
        sm = rand.standard_normal(size=(self.npca-1,self.ncoeffs))
        eval0 = self._eval_pcas(pca, ampl, sm)

        self.sm.normalizePcaImages(pca, ampl, sm)
        eval1 = self._eval_pcas(pca, ampl, sm)

        self.assertLess(np.sum((eval0-eval1)**2)**0.5, 1.0e-10)

        t = np.dot(pca, np.transpose(pca)) - np.identity(self.npca)
        self.assertLess(np.sum(t**2)**0.5, 1.0e-10)
        

    def test_spatial_model_polynomial(self):
        for n in xrange(1,50):
            self.ncand = rand.randint(150, 200)
            self.nfunc = rand.randint(1, 3)
            self.npca = rand.randint(1, 4)
            self.npix = rand.randint(100, 200)
            self.xmin = rand.uniform(0.1, 0.4)
            self.xmax = rand.uniform(0.6, 0.9)
            self.ymin = rand.uniform(0.1, 0.4)
            self.ymax = rand.uniform(0.6, 0.9)

            self.xy = np.zeros((self.ncand,2))
            self.xy[:,0] = rand.uniform(low=self.xmin, high=self.xmax, size=self.ncand)
            self.xy[:,1] = rand.uniform(low=self.ymin, high=self.ymax, size=self.ncand)

            self.construct_sm()
            self.ncoeffs = self.sm.getNcoeffs()

            self._test_eval()
            self._test_optimize()
            self._test_normalize()


class SpatialModelPolynomialTestCase(SpatialModelTestBase):
    def construct_sm(self):
        self.order = rand.randint(1,5)
        self.sm = hscpsfLib.HscSpatialModelPolynomial(self.order, self.xmin, self.xmax, self.ymin, self.ymax)


    def eval2(self, xy, fcoeffs):
        assert xy.shape == (self.ncand, 2)
        assert fcoeffs.shape == (self.nfunc, self.ncoeffs)

        ret = np.zeros((self.ncand, self.nfunc))
        tx = (xy[:,0] - self.xmin) / (self.xmax - self.xmin) - 0.5
        ty = (xy[:,1] - self.ymin) / (self.ymax - self.ymin) - 0.5

        ic = 0
        for iy in xrange(self.order+1):
            for ix in xrange(self.order-iy+1):
                ret += np.outer(tx**ix * ty**iy, fcoeffs[:,ic])
                ic += 1

        return ret


class SpatialModelLegendrePolynomialTestCase(SpatialModelTestBase):
    def construct_sm(self):
        self.order = rand.randint(1,5)
        self.sm = hscpsfLib.HscSpatialModelLegendrePolynomial(self.order, self.xmin, self.xmax, self.ymin, self.ymax)


    def eval2(self, xy, fcoeffs):
        assert xy.shape == (self.ncand, 2)
        assert fcoeffs.shape == (self.nfunc, self.ncoeffs)

        pxmat = np.zeros((self.ncand, self.order+1))
        pymat = np.zeros((self.ncand, self.order+1))

        for i in xrange(self.order+1):
            coeffs = np.zeros(i+1)
            coeffs[-1] = 1

            p = np.polynomial.legendre.Legendre(coeffs, domain=(self.xmin,self.xmax))
            pxmat[:,i] = p(xy[:,0])

            p = np.polynomial.legendre.Legendre(coeffs, domain=(self.ymin,self.ymax))
            pymat[:,i] = p(xy[:,1])

        smat = np.zeros((self.ncand, self.ncoeffs))
                
        m = 0
        for n in xrange(self.order+1):
            for i in xrange(n+1):
                smat[:,m] = pxmat[:,n-i] * pymat[:,i]
                m += 1
                
        assert m == self.ncoeffs

        return np.dot(smat, np.transpose(fcoeffs))



#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Returns a suite containing all the test cases in this module."""
    tests.init()

    suites = []
    suites += unittest.makeSuite(SpatialModelPolynomialTestCase)
    suites += unittest.makeSuite(SpatialModelLegendrePolynomialTestCase)
    return unittest.TestSuite(suites)

def run(shouldExit=False):
    """Run the tests"""
    tests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
