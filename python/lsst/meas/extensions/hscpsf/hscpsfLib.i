// -*- lsst-c++ -*-

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
 
%module(package="lsst.meas.extensions.hscpsf") hscpsfLib

namespace lsst { namespace meas { namespace extensions { namespace hscpsf { } } } }

%include "lsst/p_lsstSwig.i"
%include "lsst/base.h"     // PTR()

%{
#include "lsst/afw/detection.h"
#include "lsst/meas/algorithms.h"
#include "lsst/pex/logging.h"
#include "lsst/meas/extensions/hscpsf/hscPsf.h"

#define PY_ARRAY_UNIQUE_SYMBOL LSST_MEAS_EXTENSIONS_HSCPSF_NUMPY_ARRAY_API
#include "numpy/arrayobject.h"
#include "ndarray/swig.h"
%}

%init %{
     import_array();
%}

%lsst_exceptions();

%import "lsst/afw/detection/detectionLib.i"
%import "lsst/meas/algorithms/algorithmsLib.i"

%shared_ptr(lsst::meas::extensions::hscpsf::HscCandidateSet);

%shared_ptr(lsst::meas::extensions::hscpsf::HscPsfBase);
%shared_ptr(lsst::meas::extensions::hscpsf::HscSplinePsfBase);
%shared_ptr(lsst::meas::extensions::hscpsf::HscClipPsf);
%shared_ptr(lsst::meas::extensions::hscpsf::HscGlobalSplinePsf);
%shared_ptr(lsst::meas::extensions::hscpsf::HscPcaPsfBase);
%shared_ptr(lsst::meas::extensions::hscpsf::HscPcaPsfNoSM);
%shared_ptr(lsst::meas::extensions::hscpsf::HscPcaPsf);

%shared_ptr(lsst::meas::extensions::hscpsf::HscSpatialModelBase)
%shared_ptr(lsst::meas::extensions::hscpsf::HscSpatialModelPolynomial)

%include "lsst/meas/extensions/hscpsf/hscPsf.h"

