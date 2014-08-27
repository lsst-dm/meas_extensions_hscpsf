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
 
#include "lsst/afw/image/MaskedImage.h"
#include "lsst/meas/extensions/hscpsf/hscPsf.h"

namespace afwImage = lsst::afw::image;


namespace lsst { namespace meas { namespace extensions { namespace hscpsf {
#if 0
}}}}
#endif


HscCandidateSet::HscCandidateSet(afwImage::MaskPixel mask_bits, int nx, int ny) 
    : _mask_bits(mask_bits), _nx(nx), _ny(ny), _ncand(0), _ncand_alloc(0)
{ 
    if (nx <= 0 || ny <= 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "non-positive (nx,ny) in HscCandidateSet constructor");
}


void HscCandidateSet::add(const meas::algorithms::PsfCandidate<float> &cand, int id, double flux, double size)
{
    CONST_PTR(afwImage::MaskedImage<float>) masked_image;

    try {
        masked_image = cand.getMaskedImage(_nx, _ny);
    }
    catch (...) {
        //
        // FIXME: getMaskedImage() will fail for candidates which are near the edge of the
        // Exposure.  The best fix for this would be to modify Psf::getMaskedImage() so that
        // it returns missing pixels with some sort of extrapolated value + BAD mask bits,
        // but for now we just skip the candidate if getMaskedImage() fails
        //
        std::cerr << "HscPsf: skipped candidate at (x,y)=(" << cand.getXCenter() << "," << cand.getYCenter() << "), probably near the edge\n";
        return;
    }
    
    CONST_PTR(afwImage::Image<float>) image = masked_image->getImage();
    CONST_PTR(afwImage::Image<afwImage::VariancePixel>) variance = masked_image->getVariance();
    CONST_PTR(afwImage::Mask<afwImage::MaskPixel>) mask = masked_image->getMask();

    double chi2 = 0.0;
    int nmasked = 0;
    int nunmasked = 0;

    _reserve(_ncand+1);

    for (int i = 0; i < _nx; i++) {
        for (int j = 0; j < _ny; j++) {
            if ((*mask)(i,j) & _mask_bits) {
                // masked pixels are represented by zero inverse variance, i.e. infinitely noisy
                _im[_ncand*_nx*_ny + i*_ny + j] = 0.0;
                _iv[_ncand*_nx*_ny + i*_ny + j] = 0.0;
                nmasked++;
            }
            else {
                double im = (*image)(i,j);
                double var = (*variance)(i,j);

                if (var <= 0.0)
                    throw LSST_EXCEPT(lsst::pex::exceptions::RuntimeErrorException, "non-positive variance in PsfCandidate");

                _im[_ncand*_nx*_ny + i*_ny + j] = im;
                _iv[_ncand*_nx*_ny + i*_ny + j] = 1.0 / var;
                chi2 += im*im / var;
                nunmasked++;
            }
        }
    }

    if (nunmasked == 0) {
        // I don't think this ever happens, but just to avoid weird corner cases...
        std::cerr << "HscPsf: all pixels in image are masked, skipping candidate\n";
        return;
    }

    if (0) {
        std::cerr << "HscCandidateSet: added image " << id << "/" << _ncand << ": " 
                  << nmasked << "/" << (nmasked+nunmasked) << " pixels masked\n";
    }

    _xy[_ncand*2] = cand.getXCenter();
    _xy[_ncand*2+1] = cand.getYCenter();
    _xy0[_ncand*2] = masked_image->getX0();
    _xy0[_ncand*2+1] = masked_image->getY0();
    _chi2[_ncand] = chi2;
    _ndof[_ncand] = nunmasked;
    _id[_ncand] = id;
    _flux[_ncand] = flux;
    _size[_ncand] = size;

    if (0) {
        std::cerr << "image\n";
        for (int i = 0; i < _nx; i++) {
            for (int j = 0; j < _ny; j++)
                std::cerr << " " << _im[_ncand*_nx*_ny + i*_ny + j];
            std::cerr << "\n";
        }
        
        std::cerr << "iv\n";
        for (int i = 0; i < _nx; i++) {
            for (int j = 0; j < _ny; j++)
                std::cerr << " " << _iv[_ncand*_nx*_ny + i*_ny + j];
            std::cerr << "\n";
        }

        std::cerr << "xy " << _xy[_ncand*2] << " " << _xy[_ncand*2+1] << "\n";
        std::cerr << "xy0 " << _xy0[_ncand*2] << " " << _xy0[_ncand*2+1] << "\n";
        std::cerr << "id " << _id[_ncand] << "\n";
    }

    _ncand++;
}


void HscCandidateSet::add(CONST_PTR(HscCandidateSet) cs, int index, double x, double y)
{
    if (index < 0 || index >= cs->_ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid index in HscCandidateSet::add()");
    if ((cs->_nx != _nx) || (cs->_ny != _ny))
        throw LSST_EXCEPT(pex::exceptions::RuntimeErrorException, "(nx,ny) mismatch in HscCandidateSet::add()");

    // FIXME should sanity check (x,y)

    _reserve(_ncand+1);

    memcpy(&_im[_ncand*_nx*_ny], &cs->_im[index*_nx*_ny], _nx*_ny * sizeof(double));
    memcpy(&_iv[_ncand*_nx*_ny], &cs->_iv[index*_nx*_ny], _nx*_ny * sizeof(double));

    _xy[2*_ncand] = x;
    _xy[2*_ncand+1] = y;
    _xy0[2*_ncand] = cs->_xy0[2*index];
    _xy0[2*_ncand+1] = cs->_xy0[2*index+1];
    _chi2[_ncand] = cs->_chi2[index];
    _ndof[_ncand] = cs->_ndof[index];
    _id[_ncand] = cs->_id[index];
    _flux[_ncand] = cs->_flux[index];
    _size[_ncand] = cs->_size[index];
    _ncand++;
}


void HscCandidateSet::add(CONST_PTR(HscCandidateSet) cs, int index)
{
    if (index < 0 || index >= cs->_ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid index in HscCandidateSet::add()");
    this->add(cs, index, cs->_xy[2*index], cs->_xy[2*index+1]);
}


const double *HscCandidateSet::getIm() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getIm() called on empty CandidateSet");
    return &_im[0];
}

const double *HscCandidateSet::getIV() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getIV() called on empty CandidateSet");
    return &_iv[0];
}

const double *HscCandidateSet::getXY() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getXY() called on empty CandidateSet");
    return &_xy[0];
}

const int *HscCandidateSet::getXY0() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getXY0() called on empty CandidateSet");
    return &_xy0[0];
}

const double *HscCandidateSet::getChi2() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getChi2() called on empty CandidateSet");
    return &_chi2[0];
}

const int *HscCandidateSet::getNdof() const
{
    if (_ncand == 0)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getNdof() called on empty CandidateSet");
    return &_ndof[0];
}

double HscCandidateSet::getX(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getX() called with icand out of range");
    return _xy[2*icand];
}

double HscCandidateSet::getY(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getY() called with icand out of range");
    return _xy[2*icand + 1];
}

int HscCandidateSet::getX0(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getX0() called with icand out of range");
    return _xy0[2*icand];
}

int HscCandidateSet::getY0(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getY0() called with icand out of range");
    return _xy0[2*icand + 1];
}

int HscCandidateSet::getId(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::getId() called with icand out of range");
    return _id[icand];
}

double HscCandidateSet::getFlux(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscCandidateSet::getFlux()");
    return _flux[icand];
}

double HscCandidateSet::getSize(int icand) const
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(pex::exceptions::InvalidParameterException, "invalid candidate index in HscCandidateSet::getSize()");
    return _size[icand];
}

void HscCandidateSet::setXY(int icand, double x, double y)
{
    if (icand < 0 || icand >= _ncand)
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::setXY() called with icand out of range");
    if ((x < _xy0[2*icand]-0.5) || (x > _xy0[2*icand]+_nx-0.5) || (y < _xy0[2*icand+1]-0.5) || (y > _xy0[2*icand+1]+_ny-0.5))
        throw LSST_EXCEPT(lsst::pex::exceptions::InvalidParameterException, "HscCandidateSet::setXY() called with (x,y) outside of image");

    _xy[2*icand] = x;
    _xy[2*icand+1] = y;
}


void HscCandidateSet::_reserve(int new_ncand)
{
    assert(new_ncand > _ncand);

    if (_ncand_alloc < new_ncand) {
        _ncand_alloc = std::max(new_ncand+16, 2*_ncand);
        _im.resize(_ncand_alloc * _nx * _ny, 0);
        _iv.resize(_ncand_alloc * _nx * _ny, 0);
        _xy.resize(_ncand_alloc * 2, 0);
        _xy0.resize(_ncand_alloc * 2, 0);
        _chi2.resize(_ncand_alloc, 0);
        _ndof.resize(_ncand_alloc, 0);
        _id.resize(_ncand_alloc, -1);
        _flux.resize(_ncand_alloc, 0.0);
        _size.resize(_ncand_alloc, 0.0);
    }
}



}}}}   // namespace lsst::meas::extensions::hscpsf
