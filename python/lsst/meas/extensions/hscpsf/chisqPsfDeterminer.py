import os
import sys
import errno
import numpy as np

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import lsst.meas.algorithms.psfDeterminerRegistry as PsfDeterminerRegistry


class ChisqPsfDeterminerConfig(pexConfig.Config):
    psfName = pexConfig.Field(
        doc = "name of PSF being evaluated",
        dtype = str
    )
    outDir = pexConfig.Field(
        doc = "output directory, should be different for each psf determiner being evaluated (FIXME: this is a hack)",
        dtype = str
    )
    iterative = pexConfig.Field(
        doc = "if True, will loop over candidates and fit PSF with the N-th candidate excluded.  A better test (penalizes overfitting rather than rewarding) but much slower!",
        dtype = bool
    )


class ChisqPsfDeterminer(object):
    ConfigClass = ChisqPsfDeterminerConfig
    ConfigHack = None


    def __init__(self, config):
        if self.ConfigHack is None:
            # FIXME get rid of this hack
            raise RuntimeError('config file needs to call chisqPsfDeterminer.configHack(root) at the end')

        self.config = config
        self.real_determiner_class = PsfDeterminerRegistry[self.config.psfName]
        assert isinstance(self.ConfigHack, self.real_determiner_class.ConfigClass)


    def determinePsf(self, exposure, psfCandidateList, metadata=None, flagKey=None):            
        assert len(psfCandidateList) >= 2
        ncand = len(psfCandidateList)

        # create output directory
        try:
            os.makedirs(self.config.outDir)
        except OSError as exc:
            # OK if directory already exists
            if (exc.errno != errno.EEXIST) or not os.path.isdir(self.config.outDir):
                raise

        outid = self.exposure_id(exposure)
        filename = os.path.join(self.config.outDir, ('%s.txt' % outid))
        print >>sys.stderr, 'ChisqPsfDeterminer: output filename =', filename

        if os.path.exists(filename):
            print >>sys.stderr, 'ChisqPsfDeterminer: old version already exists, renaming %s -> %s.old' % (filename,filename)
            os.rename(filename, filename+'.old')
    
        self.orig_status = [ cand.getStatus() for cand in psfCandidateList ]
        self.status_a = [ afwMath.SpatialCellCandidate.BAD for cand in psfCandidateList ]
        self.status_b = [ afwMath.SpatialCellCandidate.BAD for cand in psfCandidateList ]
        self.chi2_a = np.array([ -1.0e9 for cand in psfCandidateList ])
        self.chi2_b = np.array([ -1.0e9 for cand in psfCandidateList ])

        #
        # Run once with all candidates.
        # Note: we set flagKey=None here
        #
        if self.config.iterative:
            # the iterative case generates so much log output that a few grep-hooks are useful
            print
            print '**** chisqPsfDeterminer: first run starts here ****'
            print

        flagKey0 = flagKey if not self.config.iterative else None
        (psf0, psfCellSet0) = self.run_psf_determiner(exposure, psfCandidateList, metadata, flagKey0)
        # save_xcenter = [ cand.getXCenter() for cand in psfCandidateList ]
        # save_ycenter = [ cand.getYCenter() for cand in psfCandidateList ]

        for (i,cand) in enumerate(psfCandidateList):
            try:
                self.chi2_a[i] = self.best_fit_residual_chi2(psf0, cand, cand.getXCenter(), cand.getYCenter())
            except:
                continue

            self.status_a[i] = cand.getStatus()


        if self.config.iterative:
            for (i,cand) in enumerate(psfCandidateList):
                for (c,s) in zip(psfCandidateList, self.orig_status):
                    c.setStatus(s)

                # for (cand, status, xcenter, ycenter) in zip(psfCandidateList, self.orig_status, save_xcenter, save_ycenter):
                #     cand.setStatus(status)
                #     cand.setXCenter(xcenter)
                #     cand.setYCenter(ycenter)
                #
                # FIXME how to get this config?  
                #   root.calibrate.measurePsf.psfDeterminer[self.config.psfName]

                print
                print '**** chisqPsfDeterminer: iterative run %d/%d starts here ****' % (i, ncand)
                print

                (psf1, psfCellSet1) = self.run_psf_determiner(exposure, psfCandidateList[:i] + psfCandidateList[(i+1):], metadata, flagKey=None)

                try:
                    self.chi2_b[i] = self.best_fit_residual_chi2(psf1, cand, cand.getXCenter(), cand.getYCenter())
                except:
                    continue

                self.status_b[i] = cand.getStatus()

            for (c,s) in zip(psfCandidateList, self.orig_status):
                c.setStatus(s)

            # for (cand, status, xcenter, ycenter) in zip(psfCandidateList, self.orig_status, save_xcenter, save_ycenter):
            #     cand.setStatus(status)
            #     cand.setXCenter(xcenter)
            #     cand.setYCenter(ycenter)

            #
            # Final run with all candidates, to set metadata and flagKey
            #
            print
            print '**** chisqPsfDeterminer: last run starts here ****'
            print

            (psf0, psfCellSet0) = self.run_psf_determiner(exposure, psfCandidateList, metadata, flagKey)


        f_out = open(filename, 'w')

        for i in xrange(ncand):
            print >>f_out, i,
            print >>f_out, self.status_str(self.status_a[i]),
            print >>f_out, self.chi2_a[i],
            print >>f_out, self.status_str(self.status_b[i]),
            print >>f_out, self.chi2_b[i]

        f_out.close()
        print >>sys.stderr, 'wrote', filename

        return (psf0, psfCellSet0)


    def exposure_id(self, exposure):
        """FIXME how to do this properly?"""

        x = exposure.getInfo().getMetadata()
        det_id = x.get('DET-ID')

        # FIXME reverse engineer the ccd from the det_id; there must be a better way to do this!                                                                                                               h = { 107:105, 112:106, 113:107, 114:108, 115:109, 108:110,111:111 }
        if 0 <= det_id < 105:
            ccd = det_id
        elif h.has_key(det_id):
            ccd = h[det_id]
        else:
            raise RuntimeError('unrecognized det_id: %d' % det_id)

        outid = '%s_%s_%d' % (x.get('OBJECT'), x.get('EXP-ID'), ccd)
        outid = '_'.join(outid.split())  # replace spaces with underscores
        return outid


    def status_str(self, status):
        if status == afwMath.SpatialCellCandidate.GOOD:
            return 'GOOD'
        elif status == afwMath.SpatialCellCandidate.BAD:
            return 'BAD'
        return 'UNKNOWN'


    def run_psf_determiner(self, exposure, psfCandidateList, metadata, flagKey):
        real_determiner = self.real_determiner_class(self.ConfigHack)
        (psf, psfCellSet) = real_determiner.determinePsf(exposure, psfCandidateList, metadata, flagKey)
        return (psf, psfCellSet)


    def best_fit_residual_chi2(self, psf, cand, x, y):
        """The point (x,y) is just used as a starting point for chi^2 minimization."""

        for n in xrange(5):
            m = gr_1d_minimizer(lambda xp: self.residual_chi2_at_xy(psf, cand, xp, y),
                                xmin=x-2, xmax=x+2, xres=0.01, xstep=0.2)

            x = m.minimize(x, allow_endpoints=False)

            m = gr_1d_minimizer(lambda yp: self.residual_chi2_at_xy(psf, cand, x, yp),
                                xmin=y-2, xmax=y+2, xres=0.01, xstep=0.2)

            y = m.minimize(y, allow_endpoints=False)
            
        return self.residual_chi2_at_xy(psf, cand, x, y)


    def residual_chi2_at_xy(self, psf, cand, x, y):
        ctr = afwGeom.Point2D(x, y)
        psf_image = psf.computeImage(ctr)

        #
        # FIXME: getMaskedImage() will fail for candidates which are near the edge of the
        # Exposure.  The best fix for this would be to modify Psf::getMaskedImage() so that
        # it returns missing pixels with some sort of extrapolated value + BAD mask bits,
        # but for now we just skip the candidate if getMaskedImage() fails.
        #
        # We deliberately return a crazy value (-1.0e9), so that any postprocessing code which
        # fails to test for a negative value and skip it will fail in a way which is obvious.
        #
        try:
            masked_image = cand.getMaskedImage(psf_image.getWidth(), psf_image.getHeight())
        except:
            print >>sys.stderr, 'chisqPsfDeterminer: skipped candidate at (x,y)=(%s,%s), probably near the edge' % (cand.getXCenter(), cand.getYCenter())
            return -1.0e-9
        return self.compute_residual_chi2(masked_image.getImage(), self.get_inverse_variance(masked_image), psf.computeImage(ctr))


    def get_inverse_variance(self, masked_image):
        mask_bits = afwImage.MaskU_getPlaneBitMask("BAD")
        mask_bits = mask_bits | afwImage.MaskU_getPlaneBitMask("CR")
        mask_bits = mask_bits | afwImage.MaskU_getPlaneBitMask("INTRP")

        iv_img = afwImage.ImageD(masked_image.getWidth(), masked_image.getHeight())
        iv_img.setXY0(masked_image.getXY0())

        iv_img.getArray()[:,:] = np.where((masked_image.getMask().getArray() & mask_bits) == 0,
                                          1.0 / masked_image.getVariance().getArray(),
                                          0.0)

        return iv_img

    
    def overlap_images(self, im_a, im_b):
        """Returns (subimage_array_a, subimage_array_b).  NOTE: subimages are transposed"""

        (x0_a, y0_a) = (im_a.getX0(), im_a.getY0())
        (x0_b, y0_b) = (im_b.getX0(), im_b.getY0())
            
        (x1_a, y1_a) = (x0_a + im_a.getWidth(), y0_a + im_a.getHeight())
        (x1_b, y1_b) = (x0_b + im_b.getWidth(), y0_b + im_b.getHeight())

        x0_overlap = max(x0_a, x0_b)
        y0_overlap = max(y0_a, y0_b)
        x1_overlap = min(x1_a, x1_b)
        y1_overlap = min(y1_a, y1_b)

        assert (x0_overlap < x1_overlap) and (y0_overlap < y1_overlap)

        subimage_array_a = im_a.getArray()[(y0_overlap-y0_a):(y1_overlap-y0_a),(x0_overlap-x0_a):(x1_overlap-x0_a)]
        subimage_array_b = im_b.getArray()[(y0_overlap-y0_b):(y1_overlap-y0_b),(x0_overlap-x0_b):(x1_overlap-x0_b)]
        return (subimage_array_a, subimage_array_b)


    def compute_residual_chi2(self, cand_img, iv_img, psf_img):
        assert cand_img.getXY0() == iv_img.getXY0()
        assert cand_img.getWidth() == iv_img.getWidth()
        assert cand_img.getHeight() == cand_img.getHeight()

        (psf_subarr, cand_subarr) = self.overlap_images(psf_img, cand_img)
        (psf_subarr, iv_subarr) = self.overlap_images(psf_img, iv_img)
        
        a = np.sum(iv_subarr * psf_subarr**2)
        b = np.sum(iv_subarr * psf_subarr * cand_subarr)
        c = np.sum(iv_img.getArray() * cand_img.getArray()**2)
        ndof = np.sum(iv_img.getArray() > 0)
        
        # reduced chi^2
        return (c - b**2/a) / ndof



class gr_1d_minimizer:
    """
    Note that this class is essentially a Python rewrite of its C++ counterpart in HscPsfUtil.h,
    so changes made to either should be reflected in the counterpart.
    """

    warn_count = 0

    def __init__(self, f, xmin, xmax, xres, xstep):
        assert xmin < xmax
        assert xres > 0.0
        assert xstep > 0.0

        self.f = f
        self.xmin = xmin
        self.xmax = xmax
        self.xres = xres
        self.xstep = xstep


    def _bisect_once(self):
        """One iteration of golden ratio minimization"""

        if (self.x1 - self.x0) > (self.x2 - self.x1):
            # trial point on low side
            x = self.x1 - 0.38197 * (self.x1 - self.x0)
            y = self.f(x)
            if y < self.y1:
                (self.x1, self.y1, self.x2, self.y2) = (x, y, self.x1, self.y1)
            else:
                (self.x0, self.y0) = (x, y)
        else:
            # high side
            x = self.x1 + 0.38197 * (self.x2 - self.x1)
            y = self.f(x)
            if y < self.y1:
                (self.x0, self.y0, self.x1, self.y1) = (self.x1, self.y1, x, y)
            else:
                (self.x2, self.y2) = (x, y)

        assert self.x0 < self.x1 < self.x2
        assert (self.y1 <= self.y0) and (self.y1 <= self.y2)


    def _bisect_to_convergence(self):
        assert self.x0 < self.x1 < self.x2
        assert (self.y1 <= self.y0) and (self.y1 <= self.y2)

        while (self.x2 - self.x0) > self.xres:
            self._bisect_once()


    def _search_down(self):
        """Returns 0 if bracketing successful, -1 if xmin was reached."""

        while True:
            assert self.xmin <= self.x0 < self.x1
            assert self.y0 < self.y1
	    
	    # special case where xmin is reached
            if self.x0 == self.xmin:
                (self.x2, self.y2) = (self.x1, self.y1)
                self.x1 = min(self.xmin+self.xres, (self.x0+self.x2)/2.)
                self.y1 = self.f(self.x1)
                return -1 if (self.y0 < self.y1) else 0

            (self.x1, self.y1, self.x2, self.y2) = (self.x0, self.y0, self.x1, self.y1)
            self.x0 = max(self.xmin, self.x1 - 1.6*(self.x2-self.x1))
            self.y0 = self.f(self.x0)

	    if self.y0 >= self.y1:
		return 0


    def _search_up(self):
        """Returns 0 if bracketing successful, 1 if xmax was reached."""

        while True:
            assert self.x1 < self.x2 <= self.xmax
            assert self.y2 < self.y1

	    # special case where xmax is reached
	    if self.x2 == self.xmax:
                (self.x0, self.y0) = (self.x1, self.y1)
                self.x1 = max(self.xmax-self.xres, (self.x0+self.x2)/2.)
                self.y1 = self.f(self.x1)
                return 1 if (self.y2 < self.y1) else 0

            (self.x0, self.y0, self.x1, self.y1) = (self.x1, self.y1, self.x2, self.y2)
            self.x2 = min(self.xmax, self.x1 + 1.6*(self.x1-self.x0))
            self.y2 = self.f(self.x2)

	    if self.y2 >= self.y1:
                return 0

    
    def _prebracket(self, xini):
        assert self.xmin <= xini <= self.xmax

	# adjust xini in corner cases
	xmid = (self.xmin + self.xmax) / 2.
	xini = max(xini, min(xmid,self.xmin+0.99*self.xstep))
	xini = min(xini, max(xmid,self.xmax-0.99*self.xstep))

	self.x0 = max(self.xmin, xini - self.xstep)
	self.x1 = xini;
	self.x2 = min(self.xmax, xini + self.xstep)
	
        (self.y0, self.y1, self.y2) = (self.f(self.x0), self.f(self.x1), self.f(self.x2))

	if self.y0 < self.y2:
            return 0 if (self.y1 <= self.y0) else self._search_down()
        else:
            return 0 if (self.y1 <= self.y2) else self._search_up()


    def minimize(self, x, allow_endpoints=False):
        assert self.xmin <= x <= self.xmax
        
        res = self._prebracket(x)
        
        if (res != 0) and not allow_endpoints:
            raise RuntimeError('gr_1d_minimizer: called with allow_endpoints=false and endpoint was reached')

	if res > 0:
	    return 0.5 * (self.x1 + self.x2)
	if res < 0:
            return 0.5 * (self.x0 + self.x1)

        if (self.y0 == self.y1) or (self.y1 == self.y2):
            if self.warn_count <= 3:
                print >>sys.stderr, 'gr_1d_minimizer: warning: tie in _prebracket(), suspect xstep is too small'
            self.warn_count += 1

        self._bisect_to_convergence();
	return 0.5 * (self.x0 + self.x2)



def configHack(root):
    """FIXME make this happen automatically, and get rid of the need to call this in the config file!"""

    k = root.calibrate.measurePsf.psfDeterminer['chisq'].psfName    
    ChisqPsfDeterminer.ConfigHack = root.calibrate.measurePsf.psfDeterminer[k]


PsfDeterminerRegistry.register("chisq", ChisqPsfDeterminer)
