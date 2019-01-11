#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE PsfCaching
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#include "boost/test/unit_test.hpp"
#pragma clang diagnostic pop

#include <iostream>
#include <iterator>
#include <algorithm>
#include <map>
#include <cmath>

#include "boost/filesystem.hpp"

#include "ndarray/eigen.h"
#include "lsst/geom/Point.h"
#include "lsst/afw/math/Kernel.h"
#include "lsst/afw/math/FunctionLibrary.h"
#include "lsst/afw/detection/Psf.h"
#include "lsst/meas/algorithms/KernelPsf.h"
#include "lsst/meas/algorithms/DoubleGaussianPsf.h"

BOOST_AUTO_TEST_CASE(FixedPsfCaching) {
    using namespace lsst::afw::detection;
    using namespace lsst::afw::geom;
    using namespace lsst::afw::image;
    using namespace lsst::meas::algorithms;
    DoubleGaussianPsf psf(7, 7, 1.5, 3.0, 0.2);
    std::shared_ptr<Psf::Image> im1 = psf.computeKernelImage(lsst::geom::Point2D(0, 0), Color(), Psf::INTERNAL);
    std::shared_ptr<Psf::Image> im2 = psf.computeImage(lsst::geom::Point2D(0, 0), Color(), Psf::INTERNAL);
    BOOST_CHECK_CLOSE(ndarray::asEigenMatrix(im1->getArray()).sum(), 1.0, 1E-8);
    BOOST_CHECK_EQUAL(ndarray::asEigenMatrix(im1->getArray()), ndarray::asEigenMatrix(im2->getArray()));
    std::shared_ptr<Psf::Image> im3 = psf.computeKernelImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im1 == im3);
    std::shared_ptr<Psf::Image> im4 = psf.computeKernelImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im3 == im4);
    std::shared_ptr<Psf::Image> im5 = psf.computeImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im2 != im5);
    std::shared_ptr<Psf::Image> im6 = psf.computeImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im5 == im6);
}

BOOST_AUTO_TEST_CASE(VariablePsfCaching) {
    using namespace lsst::afw::detection;
    using namespace lsst::afw::geom;
    using namespace lsst::afw::math;
    using namespace lsst::afw::image;
    using namespace lsst::meas::algorithms;
    std::vector<std::shared_ptr<Kernel::SpatialFunction>> spatialFuncs;
    spatialFuncs.push_back(std::make_shared< PolynomialFunction2<double> >(1));
    spatialFuncs.push_back(std::make_shared< PolynomialFunction2<double> >(1));
    spatialFuncs.push_back(std::make_shared< PolynomialFunction2<double> >(0));
    spatialFuncs[0]->setParameter(0, 1.0);
    spatialFuncs[0]->setParameter(1, 0.5);
    spatialFuncs[0]->setParameter(2, 0.5);
    spatialFuncs[1]->setParameter(0, 1.0);
    spatialFuncs[1]->setParameter(1, 0.5);
    spatialFuncs[1]->setParameter(2, 0.5);
    GaussianFunction2<double> kernelFunc(1.0, 1.0);
    AnalyticKernel kernel(7, 7, kernelFunc, spatialFuncs);
    KernelPsf psf(kernel);
    std::shared_ptr<Psf::Image> im1 = psf.computeKernelImage(lsst::geom::Point2D(0, 0), Color(), Psf::INTERNAL);
    std::shared_ptr<Psf::Image> im2 = psf.computeImage(lsst::geom::Point2D(0, 0), Color(), Psf::INTERNAL);
    BOOST_CHECK_CLOSE(ndarray::asEigenMatrix(im1->getArray()).sum(), 1.0, 1E-8);
    BOOST_CHECK_EQUAL(ndarray::asEigenMatrix(im1->getArray()), ndarray::asEigenMatrix(im2->getArray()));
    std::shared_ptr<Psf::Image> im3 = psf.computeKernelImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im1 != im3);
    std::shared_ptr<Psf::Image> im4 = psf.computeKernelImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im3 == im4);
    std::shared_ptr<Psf::Image> im5 = psf.computeImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im2 != im5);
    std::shared_ptr<Psf::Image> im6 = psf.computeImage(lsst::geom::Point2D(5, 6), Color(), Psf::INTERNAL);
    BOOST_CHECK(im5 == im6);
}
