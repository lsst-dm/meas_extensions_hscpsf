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

#ifndef LSST_MEAS_EXTENSIONS_HSCPSF_H
#define LSST_MEAS_EXTENSIONS_HSCPSF_H 1

#include "lsst/meas/algorithms/ImagePsf.h"
#include "lsst/meas/algorithms/PsfCandidate.h"


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


class HscSpatialModelBase;


class HscCandidateSet {
public:
    HscCandidateSet(afw::image::MaskPixel mask_bits, int nx, int ny);

    void add(const meas::algorithms::PsfCandidate<float> &cand, int id, double flux, double size);
    void add(CONST_PTR(HscCandidateSet) cs, int index, double x, double y);
    void add(CONST_PTR(HscCandidateSet) cs, int index);

    int getNx() const     { return _nx; }
    int getNy() const     { return _ny; }
    int getNcand() const  { return _ncand; }
    afw::image::MaskPixel getMaskBits() const { return _mask_bits; }

    const double *getIm() const;
    const double *getIV() const;
    const double *getXY() const;
    const int *getXY0() const;
    const double *getChi2() const;
    const int *getNdof() const;

    double getX(int icand) const;
    double getY(int icand) const;
    double getFlux(int icand) const;
    double getSize(int icand) const;

    int getX0(int icand) const;
    int getY0(int icand) const;
    int getId(int icand) const;

    void setXY(int icand, double x, double y);

private:
    const afw::image::MaskPixel _mask_bits;
    const int _nx;
    const int _ny;

    int _ncand;
    int _ncand_alloc;

    void _reserve(int cand_new);

    std::vector<double> _im;    // shape (ncand,nx,ny) array 
    std::vector<double> _iv;    // shape (ncand,nx,ny) array 
    std::vector<double> _xy;    // shape (ncand,2) array    [ XXX is this relative to xy0? ]
    std::vector<int>    _xy0;   // shape (ncand,2) array
    std::vector<double> _chi2;  // shape (ncand,) array; note that this is the chi^2 to zero, not the chi^2 to the PSF model!
    std::vector<int>    _ndof;  // shape (ncand,) array; number of unmasked pixels
    std::vector<int>    _id;    // shape (ncand,) array
    std::vector<double> _flux;  // shape (ncand,) array
    std::vector<double> _size;  // shape (ncand,) array
};

    
class HscPsfBase : public meas::algorithms::ImagePsf {
public:
    // @nside determines size of PSF: returned image is (2*nside+1)-by-(2*nside+1)
    HscPsfBase(CONST_PTR(HscCandidateSet) cs, int nside);

    // called when PSF is evaluated e.g. at star or galaxy position (FIXME make pure virtual)
    virtual void eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const;

    int getNcand() const { return _ncand; }
    int getNx() const { return _nx; }
    int getNy() const { return _ny; }
    int getNside() const { return _nside; }

    CONST_PTR(HscCandidateSet) getCandSet() const { return _cs; }
    
    const double *getIm(int icand) const;
    const double *getIV(int icand) const;

    double getX(int icand) const;
    double getY(int icand) const;
    double getXini(int icand) const;
    double getYini(int icand) const;

    // note: has no effect except to determine which candidates are returned by getTrimmedCandidateSet()
    void setBad(int icand);

    // returns new HscCandidateSet with updated (x,y) values and with bad candidates dropped
    PTR(HscCandidateSet) getTrimmedCandidateSet() const;

    double get_total_reduced_chi2() const;

    // devirtualize ImagePsf base class
    virtual PTR(afw::detection::Psf)  clone() const;
    virtual afw::geom::Point2D        getAveragePosition() const;
    virtual PTR(Image)                doComputeImage(afw::geom::Point2D const &position, afw::image::Color const &color) const;
    virtual PTR(Image)                doComputeKernelImage(afw::geom::Point2D const &position, afw::image::Color const &color) const;

    // static helper functions
    static void make_shear_matrix(double &axx, double &axy, double &ayy, double gamma1, double gamma2, double kappa);
    static double fit_basis_images(double *out_ampl, int nbf, int nxy, const double *iv, const double *im, const double *basis_images);   // returns residual chi^2

    //
    // Static helper function: Given a grid of values defined at integer points (i,j), where (0 <= i < nx_in) and (0 <= j < ny_in)
    // this routine interpolates onto a grid defined at points (x0+i,y0+j), where (0 <= i < nx_out) and (0 <= j < ny_out)
    //
    // Note: the transpose of 
    //    lanczos_offset_2d(x0, y0, nx_out, ny_out, nx_in, ny_in)
    // is
    //    lanczos_offset_2d(-x0, -y0, nx_in, ny_in, nx_out, ny_out)
    //
    // For this reason, we don't define a separate extirpolation-type routine (but note that it may be appropriate to
    // set the 'accumulate' flag in an extirpolation context)
    //
    static void lanczos_offset_2d(int order, double x0, double y0,
                                  int nx_out, int ny_out, double *out, int out_stride,
                                  int nx_in, int ny_in, const double *in, int in_stride, 
                                  bool zero_pad=false, bool accumulate=false);
    
protected:
    CONST_PTR(HscCandidateSet) _cs;
    const int _ncand;
    const int _nx;
    const int _ny;
    const int _nside;
    const double *_im;
    const double *_iv;
    const int *_xy0;
    const double *_initial_xy;
    const double *_detection_chi2;
    const int *_ndof;

    std::vector<double> _current_xy;       // shape (ncand, 2)
    std::vector<double> _current_psfimg;   // shape (ncand, nx, ny)
    std::vector<double> _residual_chi2;    // shape (ncand,)
    std::vector<bool>   _bad;              // shape (ncand,)

    // computed at construction
    double _xmin, _xmax, _xmean;
    double _ymin, _ymax, _ymean;
    
    //
    // Called by subclass when PSF model parameters change, to udpate _current_psfimg.  
    // Also preserves invariants by updating _residual_chi2
    // 
    void _updateCurrentPsfImg(int icand, const double *img);
};
    

class PolypixPsf : public HscPsfBase
{
public:
    PolypixPsf(CONST_PTR(HscCandidateSet) cs, int nside, int spatialOrder, double fwhm, double backnoise2, double gain);
    PolypixPsf(CONST_PTR(HscCandidateSet) cs, CONST_PTR(PolypixPsf) base);

    void psf_make(double prof_accuracy);
    void psf_makeresi(double prof_accuracy);
    void psf_clip();

    PTR(HscCandidateSet) psf_clean(double prof_accuracy);

    //
    // @out = output array of shape (_psf_nx, _psf_ny)
    // @pos = input_array of shape (2,)
    //
    void psf_build(double *loc, const double *pos) const;

    //
    // @basis = output array of shape (ncoeffs,)
    // @pos = input array of shape (2,)
    //
    void poly_func(double *basis, const double *pos) const;

    // 
    // @basis = output array of shape (ncand, ncoeffs)
    // @pos = input array of shape (ncand, 2)
    //
    void poly_eval_basis_functions(double *basis, const double *pos) const;

    //
    // @coeffs = output array of shape (ncoeffs,)
    // @data = input array of shape (ncand,)
    // @weights = input array of shape (ncand,)
    // @basis = input array of shape (ncand, ncoeffs), probably computed using poly_eval_basis_functions()
    //
    void poly_fit(double *coeffs, const double *data, const double *weights, const double *basis, double regul) const;

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

    virtual void eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const;

private:
    void _construct(int spatialOrder, double fwhm, double backnoise2, double gain);
};


extern double fast_median(double *arr, int n);

// you probably want this one...
extern void vignet_resample_xmajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                                   double dx, double dy, double step2, double stepi);

// ...not this one
extern void vignet_resample_ymajor(const double *pix1, int w1, int h1, double *pix2, int w2, int h2, 
                                   double dx, double dy, double step2, double stepi);


class HscSplinePsfBase : public HscPsfBase {
public:
    HscSplinePsfBase(CONST_PTR(HscCandidateSet) cs, int nr, double dr);

    double getGamma1(int icand) const;
    double getGamma2(int icand) const;
    double getKappa(int icand) const;
    const double *getProfile(int icand) const;    

    bool displacement_is_maximized(int icand) const;
    bool gamma_is_maximized(int icand) const;
    bool kappa_is_maximized(int icand) const;

    // Warning: optimize() should call _updateSplineParams() at the end
    virtual void optimize() = 0;

    //
    // @out = array of shape (nx,ny)
    // @profile = array of shape (nr,)
    //
    // Note: @icand is only needed for (x0,y0)
    //
    void _fillSplineImage(double *out, const double *profile, int icand, double x, double y, double gamma1, double gamma2, double kappa) const;

protected:
    int _nr;
    double _dr;

    std::vector<double>  _gamma;              // shape (ncand,2)
    std::vector<double>  _kappa;              // shape (ncand,)
    std::vector<double>  _profile;            // shape (ncand, nr)
    std::vector<double>  _basis_profiles_dd;  // shape (nr, nr+1)
    
    // Warning: these routines do not call _updateSplineParams()!
    void _optimize_xy();
    void _optimize_gamma();
    void _optimize_kappa();

    //
    // Called by subclass whenever {x,y,gamma,kappa,profile} are updated.
    //
    // Preseves invariants by computing the PSF image at each candidate and 
    // calling HscPsfBase::_updateCurrentPsfImg().
    // 
    void _updateSplineParams();
};

    
class HscClipPsf : public HscSplinePsfBase {
public:
    HscClipPsf(CONST_PTR(HscCandidateSet) cs, int nr, double dr);

    virtual void optimize();

protected:
    void _optimize_profile();
    void _dump(const char *msg, int level) const;
};

    
class HscGlobalSplinePsf : public HscSplinePsfBase {
public:
    HscGlobalSplinePsf(CONST_PTR(HscCandidateSet) cs, int nr, double dr, double sigma);

    virtual void optimize();

protected:
    std::vector<double> _global_profile;   // shape (nr,)
    std::vector<double> _candidate_ampl;   // shape (ncand,)

    // Warning: these routines do not call _updateGlobalSplineParams()
    void _optimize_global_profile();
    void _optimize_ampl();

    //
    // Called when one of { x, y, gamma, kappa, global_profile, candidate_ampl } are updated
    //
    // Preserves invariants by
    //   - normalizing the global profile
    //   - bringing HscSplinePsfBase::_profile into sync with {_global_profile, _candidate_ampl}
    //   - calling HscSplinePsfBase::_updateSplineParams()
    //
    // HscSplinePsfBase::_updateSplineParams() should always(?) be called via this routine.
    //
    void _updateGlobalSplineParams();

    void _dump(const char *msg, int level=1) const;
};


class HscPcaPsfBase : public HscPsfBase
{
public:
    // Note: PCA images are shape (2*nside+1, 2*nside+1)
    HscPcaPsfBase(CONST_PTR(HscCandidateSet) cs, int nside, int npca, int lanczos_order);

    inline int getNpca() const { return _npca; }
    inline int getNside() const { return _nside; }
    inline int getLanczosOrder() const { return _lanczos_order; }
    inline const std::vector<double> &getPcas() const { return _pca; }

protected:
    int _npca;
    int _lanczos_order;
    int _nn;       // (2*nside+1)**2

    std::vector<double>  _pca;       // shape (npca, 2*nside+1, 2*nside+1)

    void _dump(const char *msg, int level) const;

    //
    // @out = array of shape (nx,ny)
    // @in = array of shape (2*nside+1,2*nside+1)
    // 
    void _interpolate(double *out, const double *in, int icand, double x, double y) const;
    void _interpolate(double *out, const double *in, int icand) const;

    //
    // @out = array of shape (2*nside+1,2*nside+1)
    // @in = array of shape (nx,ny)
    //
    // Note: out array is accumulated; caller should initialize to zero if necessary
    //
    void _extirpolate(double *out, const double *in, int icand, double x, double y) const;
    void _extirpolate(double *out, const double *in, int icand) const;

    // 2-line helper functions to make code prettier
    void _scale_pca(double *out, double t) const;
    double _dot_pca(const double *in1, const double *in2) const;
    void _shift_pca(double *out, const double *in, double t) const;
};


class HscPcaPsfNoSM : public HscPcaPsfBase
{
public:
    HscPcaPsfNoSM(CONST_PTR(HscCandidateSet) cs, int nside, int npca, int lanczos_order);
    
    virtual void optimize();

    // callback for XY optimization
    double eval_shifted_chi2(int icand, double x, double y) const;

protected:
    std::vector<double> _pca_coeffs;      // shape (ncand,npca)

    void _init_pcas();
    void _normalize_pcas();
    void _optimize_pcas_and_update();
    void _optimize_pca_coeffs_and_update();
    void _optimize_xy();
};


class HscPcaPsf : public HscPcaPsfBase
{
public:
    HscPcaPsf(CONST_PTR(HscPcaPsfBase) psf0, CONST_PTR(HscSpatialModelBase) spatialModel);

    virtual void optimize();
    virtual void eval(int nx_out, int ny_out, double x0, double y0, double *out, double x, double y) const;

    // callback for XY optimization
    double eval_shifted_chi2(int icand, double x, double y) const;

protected:
    int _sm_ncoeffs;
    std::vector<double> _sm_coeffs;        // shape (npca-1, sm_ncoeffs)
    std::vector<double> _candidate_ampl;   // shape (ncand,)

    CONST_PTR(HscSpatialModelBase) _spatial_model;

    void _optimize_pcas_and_update();
    void _optimize_spatial_model_and_update();
    void _optimize_xy();
    void _optimize_ampl_and_update();
};


class HscSpatialModelBase {
public:
    // number of coefficents needed to represent a spatially varying scalar quantity (assumed fixed at construction)
    virtual int getNcoeffs() const = 0;

#if !defined(SWIG)
    //
    // @out = array of shape (ncand,nfunc)
    // @xy = array of shape (ncand,2)
    // @fcoeffs = array of shape (nfunc,ncoeffs)
    //
    virtual void eval(double *out, int ncand, const double *xy, int nfunc, const double *fcoeffs) const = 0;
    
    //
    // @out = 1D array of length ncoeffs
    // @xy = array of shape (ncand,2)
    // @a = 1D array of length ncand
    // @b = 1D array of length ncand
    //
    // Note: this interfaces assumes no prior on the coefficients, and will need to change in order to implement kriging
    //
    virtual void optimize(double *out, int ncand, const double *xy, const double *a, const double *b) const = 0;

    //
    // @pcas = array of shape (npca, npix)
    // @ampl = array of shape (ncand,)
    // @sm = array of shape (npca-1, ncoeffs)
    //
    // Note: not pure virtual; there is a reasonable default defined in hscSpatialModelBase.cpp
    //
    virtual void normalizePcaImages(int npca, int npix, int ncand, double *pca, double *ampl, double *sm) const;
#endif  // !defined(SWIG)

    virtual ~HscSpatialModelBase() { }

    //
    // ndarray-based versions of eval(), optimize() and normalizePcaImages()
    // Currently these are only used in python unit tests, but maybe they will become the primary interface eventually...
    //
    ndarray::Array<double,2,2> eval(const ndarray::Array<const double,2,2> &xy, const ndarray::Array<const double,2,2> &fcoeffs) const;

    ndarray::Array<double,1,1> optimize(const ndarray::Array<const double,2,2> &xy, const ndarray::Array<const double,1,1> &a, const ndarray::Array<const double,1,1> &b) const;

    void normalizePcaImages(const ndarray::Array<double,2,2> &pcas, const ndarray::Array<double,1,1> &ampl, const ndarray::Array<double,2,2> &sm) const;
};


class HscSpatialModelPolynomial : public HscSpatialModelBase {
public:
    HscSpatialModelPolynomial(int order, double xmin, double xmax, double ymin, double ymax);
    virtual ~HscSpatialModelPolynomial() { }

    // Devirtualize HscSpatialModelBase
    virtual int getNcoeffs() const;
#if !defined(SWIG)
    virtual void eval(double *out, int ncand, const double *xy, int nfunc, const double *fcoeffs) const;
    virtual void optimize(double *out, int ncand, const double *xy, const double *a, const double *b) const;
#endif

protected:
    int _order;
    int _ncoeffs;
    double _xmin, _xmax;
    double _ymin, _ymax;

    void _eval_pl_unscaled(double *out, int ncand, const double *t) const;
    void _eval_pl_scaled(double *out, int ncand, const double *xy) const;
};


}}}}   // namespace lsst::meas::extensions::hscpsf


#endif  // LSST_MEAS_EXTENSIONS_HSCPSF_H
