/*
 * LSST Data Management System
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

#include "lsst/geom/Point.h"
#include "lsst/afw/table/io/python.h"
#include "lsst/meas/algorithms/PcaPsf.h"

namespace py = pybind11;
using namespace pybind11::literals;

namespace lsst {
namespace meas {
namespace algorithms {
namespace {

PYBIND11_MODULE(pcaPsf, mod) {
    afw::table::io::python::declarePersistableFacade<PcaPsf>(mod, "PcaPsf");

    py::class_<PcaPsf, std::shared_ptr<PcaPsf>, afw::table::io::PersistableFacade<PcaPsf>, KernelPsf>
            clsPcaPsf(mod, "PcaPsf");

    clsPcaPsf.def(py::init<std::shared_ptr<afw::math::LinearCombinationKernel>, geom::Point2D const &>(),
                  "kernel"_a, "averagePosition"_a = geom::Point2D());

    clsPcaPsf.def("clone", &PcaPsf::clone);
    clsPcaPsf.def("getKernel", &PcaPsf::getKernel);
}

}  // namespace
}  // namespace algorithms
}  // namespace meas
}  // namespace lsst
