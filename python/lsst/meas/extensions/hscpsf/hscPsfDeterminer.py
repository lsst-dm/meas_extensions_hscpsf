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
# FIXME replace some hardcoded numbers with config params
# FIXME normalize PSF to sum =1


import os
import sys
import numpy as np

import lsst.pex.config as pexConfig
import lsst.pex.logging as pexLog
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms.psfDeterminerRegistry as PsfDeterminerRegistry
import lsst.meas.algorithms.algorithmsLib as algorithmsLib
import lsst.meas.extensions.hscpsf.hscpsfLib as hscpsfLib


class HscPsfDeterminerConfig(pexConfig.Config):
    nEigenComponents = pexConfig.Field(
        doc = "number of eigen components for PSF kernel creation",
        dtype = int,
        default = 4,
    )
    spatialOrder = pexConfig.Field(
        doc = "specify spatial order for PSF kernel creation",
        dtype = int,
        default = 2,
    )
    lanczosInterpolationOrder = pexConfig.Field(
        doc = "interpolation order used for offsetting images",
        dtype = int,
        default = 4
    )
    kernelSize = pexConfig.Field(
        doc = "radius of the kernel to create, relative to the square root of the stellar quadrupole moments",
        dtype = float,
        default = 10.0,
    )
    kernelSizeMin = pexConfig.Field(
        doc = "Minimum radius of the kernel",
        dtype = int,
        default = 25,
    )
    kernelSizeMax = pexConfig.Field(
        doc = "Maximum radius of the kernel",
        dtype = int,
        default = 45,
    )
    sizeCellX = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, column direction)",
        dtype = int,
        default = 256,
        check = lambda x: x >= 10,
    )
    sizeCellY = pexConfig.Field(
        doc = "size of cell used to determine PSF (pixels, row direction)",
        dtype = int,
        default = sizeCellX.default,
        check = lambda x: x >= 10,
    )


class HscPsfDeterminer(object):
    ConfigClass = HscPsfDeterminerConfig


    def __init__(self, config):
        self.config = config
        self.debugLog = pexLog.Debug("meas.algorithms.psfDeterminer")
        self.warnLog = pexLog.Log(pexLog.getDefaultLog(), "meas.algorithms.psfDeterminer")


    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """Determine a PCA PSF model for an exposure given a list of PSF candidates
        
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata  a home for interesting tidbits of information
        @param[in] flagKey: schema key used to mark sources actually used in PSF determination
    
        @return psf: an lsst.meas.algorithms.PcaPsf
        """
        
        if len(psfCandidateList) == 0:
            raise RuntimeError("No PSF candidates supplied.")

        # sets self.kernelSize
        self.determineKernelSize(psfCandidateList)

        # also marks bad candidates
        psf = self.fitPsf(psfCandidateList)

        psfCellSet = self.makePsfCellSet(exposure, psfCandidateList)
        self.setMetadata(psfCellSet, metadata, flagKey)

        return psf, psfCellSet


    def determineKernelSize(self, psfCandidateList):
        """Sets self.kernelSize (logic in this routine is cut-and-paste from pcaPsfDeterminer)"""

        sizes = np.zeros(len(psfCandidateList))

        for i, psfCandidate in enumerate(psfCandidateList):
            source = psfCandidate.getSource()
            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            axes = afwEll.Axes(quad)
            sizes[i] = axes.getA()

        if self.config.kernelSize >= 15:
            print >>sys.stderr, "WARNING: NOT scaling kernelSize by stellar quadrupole moment, but using absolute value"
            self.kernelSize = int(self.config.kernelSize)
        else:
            self.kernelSize = 2 * int(self.config.kernelSize * np.sqrt(np.median(sizes)) + 0.5) + 1
            self.kernelSize = max(self.kernelSize, self.config.kernelSizeMin)
            self.kernelSize = min(self.kernelSize, self.config.kernelSizeMax)

        print >>sys.stderr, 'setting kernelSize =', self.kernelSize


        
    def fitPsf(self, psfCandidateList):
        print >>sys.stderr, 'FIXME: resolve issue of proper bitmask (also in evalPsfDeterminer?)'
        mask_bits = afwImage.MaskU_getPlaneBitMask("BAD")
        mask_bits = mask_bits | afwImage.MaskU_getPlaneBitMask("CR")
        mask_bits = mask_bits | afwImage.MaskU_getPlaneBitMask("INTRP")

        cs = hscpsfLib.HscCandidateSet(mask_bits, self.kernelSize, self.kernelSize)

        for (i, cand) in enumerate(psfCandidateList):
            cs.add(cand, i, 0.0, 0.0)   # FIXME use real values of flux/size here, for completeness

        # FIXME rethink this
        for i in xrange(cs.getNcand()):
            (x0, y0) = (cs.getX0(i), cs.getY0(i))
            cs.setXY(i, x0 + 0.5*(cs.getNx()-1), y0 + 0.5*(cs.getNy()-1))

        clipPsf = hscpsfLib.HscClipPsf(cs, 15, 1.0)
        for n in xrange(5):
            clipPsf.optimize()

        cs = clipPsf.getTrimmedCandidateSet()

        gsPsf = hscpsfLib.HscGlobalSplinePsf(cs, 15, 1.0, 1.0)
        for n in xrange(5):
            gsPsf.optimize()

        for icand in xrange(gsPsf.getNcand()):
            if gsPsf.displacement_is_maximized(icand) or gsPsf.gamma_is_maximized(icand) or gsPsf.kappa_is_maximized(icand):
                gsPsf.setBad(icand)

        cs = gsPsf.getTrimmedCandidateSet()

        assert self.kernelSize % 2 == 1
        nside = (self.kernelSize-1) // 2

        pcaPsf = hscpsfLib.HscPcaPsfNoSM(cs, nside, self.config.nEigenComponents, self.config.lanczosInterpolationOrder)
        for n in xrange(5):
            pcaPsf.optimize()

        if self.config.nEigenComponents > 1:
            xvec = [ pcaPsf.getX(icand) for icand in xrange(pcaPsf.getNcand()) ]
            yvec = [ pcaPsf.getY(icand) for icand in xrange(pcaPsf.getNcand()) ]

            spatialModel = hscpsfLib.HscSpatialModelPolynomial(self.config.spatialOrder, min(xvec), max(xvec), min(yvec), max(yvec))
            pcaPsf = hscpsfLib.HscPcaPsf(pcaPsf, spatialModel)

            for n in xrange(5):
                pcaPsf.optimize()

        # The next 2 loops mark bad candidates
        bad = np.ones(len(psfCandidateList), dtype=bool)
        for icand in xrange(cs.getNcand()):
            bad[cs.getId(icand)] = False

        for i in xrange(len(psfCandidateList)):
            status = afwMath.SpatialCellCandidate.BAD if bad[i] else afwMath.SpatialCellCandidate.UNKNOWN
            psfCandidateList[i].setStatus(status)

        return pcaPsf


    def makePsfCellSet(self, exposure, psfCandidateList):
        """Construct and populate a spatial cell set (based on meas_algorithms/pcaPsf.py)"""

        mi = exposure.getMaskedImage()

        # construct and populate a spatial cell set
        bbox = mi.getBBox(afwImage.PARENT)
        psfCellSet = afwMath.SpatialCellSet(bbox, self.config.sizeCellX, self.config.sizeCellY)

        # FIXME: understand under which circumstances the try..except fails
        for (i, psfCandidate) in enumerate(psfCandidateList):
            try:
                psfCellSet.insertCandidate(psfCandidate)
            except Exception, e:
                self.debugLog.debug(2, "Skipping PSF candidate %d of %d: %s" % (i, len(psfCandidateList), e))
                continue

        return psfCellSet


    def setMetadata(self, psfCellSet, metadata, flagKey):
        """Cut-and-paste from meas_algorithms/pcaPsfDeterminer.py"""

        numGoodStars = 0
        numAvailStars = 0

        avgX = 0.0
        avgY = 0.0

        for cell in psfCellSet.getCellList():
            for cand in cell.begin(False):  # don't ignore BAD stars
                numAvailStars += 1

            for cand in cell.begin(True):  # do ignore BAD stars
                cand = algorithmsLib.cast_PsfCandidateF(cand)
                src = cand.getSource()
                if flagKey is not None:
                    src.set(flagKey, True)
                avgX += src.getX()
                avgY += src.getY()
                numGoodStars += 1

        avgX /= numGoodStars
        avgY /= numGoodStars

        if metadata != None:
            # FIXME figure out exactly what to put for spatialFitChi2
            # metadata.set("spatialFitChi2", fitChi2)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("numAvailStars", numAvailStars)
            metadata.set("avgX", avgX)
            metadata.set("avgY", avgY)
        

PsfDeterminerRegistry.register("hsc", HscPsfDeterminer)
