"""
Example config file:

    import comp2PsfDeterminer

    root.calibrate.measurePsf.psfDeterminer.name = 'comp2'
    root.calibrate.measurePsf.psfDeterminer['comp2'].psfName1 = 'psfex'
    root.calibrate.measurePsf.psfDeterminer['comp2'].psfName2 = 'pixpoly'
    root.calibrate.measurePsf.psfDeterminer['psfex'].spatialOrder = 2        # for example
    root.calibrate.measurePsf.psfDeterminer['pixpoly'].spatialOrder = 2      # for example

    comp2PsfDeterminer.configHack(root)     # always required
"""

import os
import sys
import numpy as np

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms.psfDeterminerRegistry as PsfDeterminerRegistry


class Comp2PsfDeterminerConfig(pexConfig.Config):
    psfName1 = pexConfig.Field(
        doc = "name of first PSF being compared",
        dtype = str
    )
    psfName2 = pexConfig.Field(
        doc = "name of second PSF being compared",
        dtype = str
    )


class Comp2PsfDeterminer(object):
    ConfigClass = Comp2PsfDeterminerConfig
    
    # set in configHack(), see end of file
    determiner_config1 = None
    determiner_config2 = None

    def __init__(self, config):
        if (self.determiner_config1 is None) or (self.determiner_config2 is None):
            raise RuntimeError('config file needs to call comp2PsfDeterminer.configHack(root) at the end')

        self.config = config
        self.determiner_class1 = PsfDeterminerRegistry[self.config.psfName1]
        self.determiner_class2 = PsfDeterminerRegistry[self.config.psfName2]

        assert isinstance(self.determiner_config1, self.determiner_class1.ConfigClass)
        assert isinstance(self.determiner_config2, self.determiner_class2.ConfigClass)


    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):
        xsrc = np.array([ c.getSource().getX() for c in psfCandidateList ])
        ysrc = np.array([ c.getSource().getY() for c in psfCandidateList ])
        status0 = [ c.getStatus() for c in psfCandidateList ]
        
        # These are (x,y) values where the psfs will be compared
        xvec = np.linspace(np.min(xsrc)+16.32, np.max(xsrc)-19.77, 5)
        yvec = np.linspace(np.min(ysrc)+20.45, np.max(ysrc)-23.89, 5)

        print '\nComp2PsfDeterminer: starting first psf run (%s)\n' % self.config.psfName1

        determiner1 = self.determiner_class1(self.determiner_config1)
        (psf1, psfCellSet1) = determiner1.determinePsf(exposure, psfCandidateList, metadata, flagKey=None)
        bad1 = [ i for (i,c) in enumerate(psfCandidateList) if (c.getStatus() == afwMath.SpatialCellCandidate.BAD) ]

        print '\nComp2PsfDeterminer: starting second psf run (%s)\n' % self.config.psfName2

        for (c,s) in zip(psfCandidateList, status0):
            c.setStatus(s)

        determiner2 = self.determiner_class2(self.determiner_config2)
        (psf2, psfCellSet2) = determiner2.determinePsf(exposure, psfCandidateList, metadata, flagKey=None)
        bad2 = [ i for (i,c) in enumerate(psfCandidateList) if (c.getStatus() == afwMath.SpatialCellCandidate.BAD) ]

        print '\nComp2PsfDeterminer: table of comparison values follows'

        for x in xvec:
            for y in yvec:
                print '(x,y) =', (x,y)

                self.compare_images('doComputeImage', 
                                    psf1.doComputeImage(afwGeom.Point2D(x,y), afwImage.Color()),
                                    psf2.doComputeImage(afwGeom.Point2D(x,y), afwImage.Color()))

                self.compare_images('doComputeKernelImage', 
                                    psf1.doComputeKernelImage(afwGeom.Point2D(x,y), afwImage.Color()),
                                    psf2.doComputeKernelImage(afwGeom.Point2D(x,y), afwImage.Color()))

        print '\nComp2PsfDeterminer: the following candidates were marked bad by %s: %s' % (self.config.psfName1, bad1)
        print 'Comp2PsfDeterminer: the following candidates were marked bad by %s: %s\n' % (self.config.psfName2, bad2)

        if bad1 != bad2:
            print '   !!! WARNING the two sets of bad candidates differ !!!\n'

        return (psf2, psfCellSet2)


    def compare_images(self, label, im1, im2):
        (w1, h1, x1, y1) = (im1.getWidth(), im1.getHeight(), im1.getXY0().getX(), im1.getXY0().getY())
        (w2, h2, x2, y2) = (im2.getWidth(), im2.getHeight(), im2.getXY0().getX(), im2.getXY0().getY())

        if (w1,h1,x1,y1) != (w2,h2,x2,y2):
            print "   %s: couldn't compare images, (w,h,x,y) values differ:  %s  %s" % (label,(w1,h1,x1,y1),(w2,h2,x2,y2))
            return

        arr1 = im1.getArray()
        arr2 = im2.getArray()
        num = np.sum((arr1-arr2)**2)
        den = (np.sum(arr1**2) * np.sum(arr2**2))**0.5
        print "   %s: %s" % (label, num/den)


def configHack(root):
    k1 = root.calibrate.measurePsf.psfDeterminer['comp2'].psfName1
    k2 = root.calibrate.measurePsf.psfDeterminer['comp2'].psfName2

    if k1 == k2:
        raise RuntimeError('comp2PsfDeterminer: fatal: psfName1 and psfName2 are not allowed to be identical')

    Comp2PsfDeterminer.determiner_config1 = root.calibrate.measurePsf.psfDeterminer[k1]
    Comp2PsfDeterminer.determiner_config2 = root.calibrate.measurePsf.psfDeterminer[k2]

PsfDeterminerRegistry.register("comp2", Comp2PsfDeterminer)
