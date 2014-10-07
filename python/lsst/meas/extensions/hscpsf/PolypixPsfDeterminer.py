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
import os
import sys
import numpy as np

import lsst.daf.base as dafBase
import lsst.pex.config as pexConfig
import lsst.pex.exceptions as pexExceptions
import lsst.pex.logging as pexLog
import lsst.afw.cameraGeom as afwCG
import lsst.afw.geom as afwGeom
import lsst.afw.geom.ellipses as afwEll
import lsst.afw.detection as afwDetection
import lsst.afw.display.ds9 as ds9
import lsst.afw.image as afwImage
import lsst.afw.table as afwTable
import lsst.afw.math as afwMath
import lsst.meas.algorithms as measAlg
import lsst.meas.algorithms.utils as maUtils
import lsst.meas.algorithms.psfDeterminerRegistry as psfDeterminerRegistry
import lsst.meas.extensions.hscpsf.hscpsfLib as hscpsfLib


class PolypixPsfDeterminerConfig(pexConfig.Config):
    spatialOrder = pexConfig.Field(
        doc = "specify spatial order for PSF kernel creation",
        dtype = int,
        default = 2,
        check = lambda x: x >= 0,
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
    kernelSize = pexConfig.Field(
        doc = "radius of the kernel to create, relative to the square root of the stellar quadrupole moments",
        dtype = float,
        default = 81.0,
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
    samplingSize = pexConfig.Field(
        doc = "Resolution of the internal PSF model relative to the pixel size; e.g. 0.5 is equal to 2x oversampling",
        dtype = float,
        default = 0.5,
    )
    badMaskBits = pexConfig.ListField(
        doc="""List of mask bits which cause a source to be rejected as bad
N.b. INTRP is used specially in PsfCandidateSet; it means "Contaminated by neighbour"
""",
        dtype=str,
        default=["INTRP", "SAT"],
    )


class PolypixPsfDeterminer(object):
    ConfigClass = PolypixPsfDeterminerConfig

    def __init__(self, config):
        """
        @param[in] config: instance of PolypixPsfDeterminerConfig
        """
        self.config = config
        # N.b. name of component is meas.algorithms.psfDeterminer so you can turn on psf debugging
        # independent of which determiner is active
        self.debugLog = pexLog.Debug("meas.algorithms.psfDeterminer")
        self.warnLog = pexLog.Log(pexLog.getDefaultLog(), "meas.algorithms.psfDeterminer")

    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        """
        @param[in] exposure: exposure containing the psf candidates (lsst.afw.image.Exposure)
        @param[in] psfCandidateList: a sequence of PSF candidates (each an lsst.meas.algorithms.PsfCandidate);
            typically obtained by detecting sources and then running them through a star selector
        @param[in,out] metadata  a home for interesting tidbits of information
        @param[in] flagKey: schema key used to mark sources actually used in PSF determination
    
        @return psf
        """

        import lsstDebug
        display = lsstDebug.Info(__name__).display 
        displayExposure = display and \
            lsstDebug.Info(__name__).displayExposure     # display the Exposure + spatialCells 
        displayPsfCandidates = display and \
            lsstDebug.Info(__name__).displayPsfCandidates # show the viable candidates 
        displayPsfComponents = display and \
            lsstDebug.Info(__name__).displayPsfComponents # show the basis functions
        showBadCandidates = display and \
            lsstDebug.Info(__name__).showBadCandidates # Include bad candidates (meaningless, methinks)
        displayResiduals = display and \
            lsstDebug.Info(__name__).displayResiduals         # show residuals
        displayPsfMosaic = display and \
            lsstDebug.Info(__name__).displayPsfMosaic   # show mosaic of reconstructed PSF(x,y)
        matchKernelAmplitudes = lsstDebug.Info(__name__).matchKernelAmplitudes # match Kernel amplitudes for spatial plots
        normalizeResiduals = lsstDebug.Info(__name__).normalizeResiduals # Normalise residuals by object amplitude 

        mi = exposure.getMaskedImage()
        
        nCand = len(psfCandidateList)
        if nCand == 0:
            raise RuntimeError("No PSF candidates supplied.")
        #
        # How big should our PSF models be?
        #
        if display:                     # only needed for debug plots
            # construct and populate a spatial cell set
            bbox = mi.getBBox(afwImage.PARENT)
            psfCellSet = afwMath.SpatialCellSet(bbox, self.config.sizeCellX, self.config.sizeCellY)
        else:
            psfCellSet = None
        
        sizes = np.empty(nCand)
        for i, psfCandidate in enumerate(psfCandidateList):
            try:
                if psfCellSet:
                    psfCellSet.insertCandidate(psfCandidate)
            except Exception, e:
                self.debugLog.debug(2, "Skipping PSF candidate %d of %d: %s" % (i, len(psfCandidateList), e))
                continue

            source = psfCandidate.getSource()
            quad = afwEll.Quadrupole(source.getIxx(), source.getIyy(), source.getIxy())
            rmsSize = quad.getTraceRadius()
            sizes[i] = rmsSize

        if self.config.kernelSize >= 15:
            self.debugLog.debug(1, \
                "WARNING: NOT scaling kernelSize by stellar quadrupole moment, but using absolute value")
            actualKernelSize = int(self.config.kernelSize)
        else:
            actualKernelSize = 2 * int(self.config.kernelSize * np.sqrt(np.median(sizes)) + 0.5) + 1
            if actualKernelSize < self.config.kernelSizeMin:
                actualKernelSize = self.config.kernelSizeMin
            if actualKernelSize > self.config.kernelSizeMax:
                actualKernelSize = self.config.kernelSizeMax
            if display:
                rms = np.median(sizes)
                print "Median PSF RMS size=%.2f pixels (\"FWHM\"=%.2f)" % (rms, 2*np.sqrt(2*np.log(2))*rms)

        # If we manually set the resolution then we need the size in pixel units
        pixKernelSize = actualKernelSize
        if self.config.samplingSize > 0:
            pixKernelSize = int(actualKernelSize*self.config.samplingSize)
            if pixKernelSize % 2 == 0: pixKernelSize += 1 
        self.debugLog.debug(3, "PSF Kernel size=%.2f, Image Kernel Size=%.2f" % 
                            (actualKernelSize,pixKernelSize))

        nside = pixKernelSize // 2
        assert pixKernelSize == (2*nside+1)

        # Set size of image returned around candidate
        psfCandidateList[0].setHeight(pixKernelSize)
        psfCandidateList[0].setWidth(pixKernelSize)

        mask_bits = afwImage.MaskU_getPlaneBitMask(self.config.badMaskBits)            
        fluxName = 'initial.flux.sinc'    # FIXME should be in config? (meas_extensions_psfex has it in a config file)
        fluxFlagName = fluxName + ".flags"

        cs = hscpsfLib.HscCandidateSet(mask_bits, pixKernelSize, pixKernelSize)

        with ds9.Buffering():
            xpos = []; ypos = []
            for i, psfCandidate in enumerate(psfCandidateList):
                source = psfCandidate.getSource()
                xc, yc = source.getX(), source.getY()

                if fluxFlagName in source.schema and source.get(fluxFlagName):
                    continue

                flux = source.get(fluxName)
                if flux < 0 or np.isnan(flux):
                    continue

                cs.add(psfCandidate, i, flux, sizes[i])

                if flagKey is not None:
                    source.set(flagKey, True)

                xpos.append(xc); ypos.append(yc) # for QA

                if displayExposure:
                    ds9.dot("o", xc, yc, ctype=ds9.CYAN, size=4, frame=frame)

        if cs.getNcand() == 0:
            raise RuntimeError("No good PSF candidates")

        # Calculate fwhm, backnoise2 and gain
        fwhm = 2*np.sqrt(2*np.log(2))*np.median(sizes)
        backnoise2 = afwMath.makeStatistics(mi.getImage(), afwMath.VARIANCECLIP).getValue()
        ccd = afwCG.cast_Ccd(exposure.getDetector())
        if ccd:
            gain = np.mean(np.array([a.getElectronicParams().getGain() for a in ccd]))
        else:
            gain = 1.0
            self.warnLog.log(pexLog.Log.WARN, "Setting gain to %g" % gain)

        psf = hscpsfLib.PolypixPsf(cs, nside, actualKernelSize, self.config.samplingSize, self.config.spatialOrder, fwhm, backnoise2, gain)
        psf.psf_make(0.2, 1000.0)
        cs = psf.psf_clean(0.2)

        psf = hscpsfLib.PolypixPsf(cs, psf)
        psf.psf_make(0.1, 1000.0)
        cs = psf.psf_clean(0.1)

        psf = hscpsfLib.PolypixPsf(cs, psf)
        psf.psf_make(0.05, 1000.0)
        cs = psf.psf_clean(0.05)

        psf = hscpsfLib.PolypixPsf(cs, psf)
        psf.psf_make(0.01, 1000.0)
        psf.psf_clip()
        cs = psf.psf_clean(0.01)

        # The next 2 loops mark bad candidates
        bad = np.ones(len(psfCandidateList), dtype=bool)
        for icand in xrange(cs.getNcand()):
            bad[cs.getId(icand)] = False

        for i in xrange(len(psfCandidateList)):
            status = afwMath.SpatialCellCandidate.BAD if bad[i] else afwMath.SpatialCellCandidate.UNKNOWN
            psfCandidateList[i].setStatus(status)

        xpos = np.array(xpos); ypos = np.array(ypos)
        numGoodStars = len(xpos)
        avgX, avgY = np.mean(xpos), np.mean(ypos)

        #
        # Generate some QA information
        #
        # Count PSF stars
        #
        if metadata != None:
            metadata.set("spatialFitChi2", np.nan)
            metadata.set("numAvailStars", nCand)
            metadata.set("numGoodStars", numGoodStars)
            metadata.set("avgX", avgX)
            metadata.set("avgY", avgY)

        psfCellSet = None
        return psf, psfCellSet


#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-    

psfDeterminerRegistry.register("polypix", PolypixPsfDeterminer)
