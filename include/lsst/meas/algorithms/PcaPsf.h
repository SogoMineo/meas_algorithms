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

#ifndef LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED
#define LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED

#include "lsst/geom/Point.h"
#include "lsst/meas/algorithms/KernelPsf.h"

namespace lsst {
namespace meas {
namespace algorithms {

/**
 * @brief Represent a PSF as a linear combination of PCA (== Karhunen-Loeve) basis functions
 */
class PcaPsf : public afw::table::io::PersistableFacade<PcaPsf>, public KernelPsf {
public:
    /**
     *  @brief Constructor for a PcaPsf
     *
     *  @param[in] kernel           Kernel that defines the Psf.
     *  @param[in] averagePosition  Average position of stars used to construct the Psf.
     */
    explicit PcaPsf(std::shared_ptr<afw::math::LinearCombinationKernel> kernel,
                    geom::Point2D const& averagePosition = geom::Point2D());

    /// Polymorphic deep copy; should usually be unnecessary as Psfs are immutable.x
    std::shared_ptr<afw::detection::Psf> clone() const override;

    /// Return a clone with specified kernel dimensions
    std::shared_ptr<afw::detection::Psf> resized(int width, int height) const override;

    /// PcaPsf always has a LinearCombinationKernel, so we can override getKernel to make it more useful.
    std::shared_ptr<afw::math::LinearCombinationKernel const> getKernel() const;

private:
    // Name used in table persistence; the rest of is implemented by KernelPsf.
    std::string getPersistenceName() const override { return "PcaPsf"; }
};

}  // namespace algorithms
}  // namespace meas
}  // namespace lsst

#endif  // !LSST_MEAS_ALGORITHMS_PcaPsf_h_INCLUDED
