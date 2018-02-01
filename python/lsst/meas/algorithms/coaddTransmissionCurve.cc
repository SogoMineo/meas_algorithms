/*
 * LSST Data Management System
 * Copyright 2018 LSST/AURA.
 *
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 * See the COPYRIGHT file
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
 * see <https://www.lsstcorp.org/LegalNotices/>.
 */
#include "pybind11/pybind11.h"

#include "lsst/meas/algorithms/CoaddTransmissionCurve.h"
#include "lsst/afw/image/Wcs.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst { namespace meas { namespace algorithms {

PYBIND11_PLUGIN(coaddTransmissionCurve) {
    py::module mod("coaddTransmissionCurve");
    py::module::import("lsst.afw.image");
    py::module::import("lsst.afw.table");
    mod.def("makeCoaddTransmissionCurve", &makeCoaddTransmissionCurve, "coaddWcs"_a, "inputSensors"_a);
    return mod.ptr();
}

}}} // lsst::meas::algorithms
