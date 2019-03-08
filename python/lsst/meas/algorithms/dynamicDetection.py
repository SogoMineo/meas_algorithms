
__all__ = ["DynamicDetectionConfig", "DynamicDetectionTask"]

import numpy as np

from lsst.pex.config import Field, ConfigurableField
from lsst.pipe.base import Struct

from .detection import SourceDetectionConfig, SourceDetectionTask
from .skyObjects import SkyObjectsTask

from lsst.afw.detection import FootprintSet
from lsst.afw.table import SourceCatalog, SourceTable
from lsst.meas.base import ForcedMeasurementTask

import lsst.afw.image
import lsst.afw.math


class DynamicDetectionConfig(SourceDetectionConfig):
    """Configuration for DynamicDetectionTask"""
    prelimThresholdFactor = Field(dtype=float, default=0.5,
                                  doc="Fraction of the threshold to use for first pass (to find sky objects)")
    skyObjects = ConfigurableField(target=SkyObjectsTask, doc="Generate sky objects")
    doBackgroundTweak = Field(dtype=bool, default=True,
                              doc="Tweak background level so median PSF flux of sky objects is zero?")
    minNumSources = Field(dtype=int, default=10,
                          doc="Minimum number of sky sources in statistical sample; "
                              "if below this number, we refuse to modify the threshold.")

    def setDefaults(self):
        SourceDetectionConfig.setDefaults(self)
        self.skyObjects.nSources = 1000  # For good statistics


class DynamicDetectionTask(SourceDetectionTask):
    """Detection of sources on an image with a dynamic threshold

    We first detect sources using a lower threshold than normal (see config
    parameter ``prelimThresholdFactor``) in order to identify good sky regions
    (configurable ``skyObjects``). Then we perform forced PSF photometry on
    those sky regions. Using those PSF flux measurements and estimated errors,
    we set the threshold so that the stdev of the measurements matches the
    median estimated error.
    """
    ConfigClass = DynamicDetectionConfig
    _DefaultName = "dynamicDetection"

    def __init__(self, *args, **kwargs):
        """Constructor

        Besides the usual initialisation of configurables, we also set up
        the forced measurement which is deliberately not represented in
        this Task's configuration parameters because we're using it as part
        of the algorithm and we don't want to allow it to be modified.
        """
        SourceDetectionTask.__init__(self, *args, **kwargs)
        self.makeSubtask("skyObjects")

        # Set up forced measurement.
        config = ForcedMeasurementTask.ConfigClass()
        config.plugins.names = ['base_TransformedCentroid', 'base_PsfFlux', 'base_LocalBackground']
        # We'll need the "centroid" and "psfFlux" slots
        for slot in ("shape", "psfShape", "apFlux", "modelFlux", "gaussianFlux", "calibFlux"):
            setattr(config.slots, slot, None)
        config.copyColumns = {}
        self.skySchema = SourceTable.makeMinimalSchema()
        self.skyMeasurement = ForcedMeasurementTask(config=config, name="skyMeasurement", parentTask=self,
                                                    refSchema=self.skySchema)

    def calculateThreshold(self, exposure, seed, sigma=None):
        """Calculate new threshold

        This is the main functional addition to the vanilla
        `SourceDetectionTask`.

        We identify sky objects and perform forced PSF photometry on
        them. Using those PSF flux measurements and estimated errors,
        we set the threshold so that the stdev of the measurements
        matches the median estimated error.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure on which we're detecting sources.
        seed : `int`
            RNG seed to use for finding sky objects.
        sigma : `float`, optional
            Gaussian sigma of smoothing kernel; if not provided,
            will be deduced from the exposure's PSF.

        Returns
        -------
        result : `lsst.pipe.base.Struct`
            Result struct with components:

            - ``multiplicative``: multiplicative factor to be applied to the
                configured detection threshold (`float`).
            - ``additive``: additive factor to be applied to the background
                level (`float`).
        """
        # Make a catalog of sky objects
        fp = self.skyObjects.run(exposure.maskedImage.mask, seed)
        skyFootprints = FootprintSet(exposure.getBBox())
        skyFootprints.setFootprints(fp)
        table = SourceTable.make(self.skyMeasurement.schema)
        catalog = SourceCatalog(table)
        catalog.reserve(len(skyFootprints.getFootprints()))
        skyFootprints.makeSources(catalog)
        key = catalog.getCentroidKey()
        for source in catalog:
            peaks = source.getFootprint().getPeaks()
            assert len(peaks) == 1
            source.set(key, peaks[0].getF())
            source.updateCoord(exposure.getWcs())

        # Forced photometry on sky objects
        self.skyMeasurement.run(catalog, exposure, catalog, exposure.getWcs())

        # Calculate new threshold
        fluxes = catalog["base_PsfFlux_instFlux"]
        area = catalog["base_PsfFlux_area"]
        bg = catalog["base_LocalBackground_instFlux"]

        good = (~catalog["base_PsfFlux_flag"] & ~catalog["base_LocalBackground_flag"] &
                np.isfinite(fluxes) & np.isfinite(area) & np.isfinite(bg))

        if good.sum() < self.config.minNumSources:
            self.log.warn("Insufficient good flux measurements (%d < %d) for dynamic threshold calculation",
                          good.sum(), self.config.minNumSources)
            return Struct(multiplicative=1.0, additive=0.0)

        bgMedian = np.median((fluxes/area)[good])

        lq, uq = np.percentile((fluxes - bg*area)[good], [25.0, 75.0])
        stdevMeas = 0.741*(uq - lq)
        medianError = np.median(catalog["base_PsfFlux_instFluxErr"][good])
        return Struct(multiplicative=medianError/stdevMeas, additive=bgMedian)

    def detectFootprints(self, exposure, doSmooth=True, sigma=None, clearMask=True, expId=None):
        """Detect footprints with a dynamic threshold

        This varies from the vanilla ``detectFootprints`` method because we
        do detection twice: one with a low threshold so that we can find
        sky uncontaminated by objects, then one more with the new calculated
        threshold.

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure to process; DETECTED{,_NEGATIVE} mask plane will be
            set in-place.
        doSmooth : `bool`, optional
            If True, smooth the image before detection using a Gaussian
            of width ``sigma``.
        sigma : `float`, optional
            Gaussian Sigma of PSF (pixels); used for smoothing and to grow
            detections; if `None` then measure the sigma of the PSF of the
            ``exposure``.
        clearMask : `bool`, optional
            Clear both DETECTED and DETECTED_NEGATIVE planes before running
            detection.
        expId : `int`, optional
            Exposure identifier, used as a seed for the random number
            generator. If absent, the seed will be the sum of the image.

        Return Struct contents
        ----------------------
        positive : `lsst.afw.detection.FootprintSet`
            Positive polarity footprints (may be `None`)
        negative : `lsst.afw.detection.FootprintSet`
            Negative polarity footprints (may be `None`)
        numPos : `int`
            Number of footprints in positive or 0 if detection polarity was
            negative.
        numNeg : `int`
            Number of footprints in negative or 0 if detection polarity was
            positive.
        background : `lsst.afw.math.BackgroundList`
            Re-estimated background.  `None` if
            ``reEstimateBackground==False``.
        factor : `float`
            Multiplication factor applied to the configured detection
            threshold.
        prelim : `lsst.pipe.base.Struct`
            Results from preliminary detection pass.
        """
        maskedImage = exposure.maskedImage
        psf = self.getPsf(exposure, sigma=sigma)

        # random seed needs to fit in a C++ 'int' so pybind doesn't choke on it
        seed = (expId if expId is not None else int(maskedImage.image.array.sum())) % (2**31 - 1)

        # We want to keep any large-scale background (e.g., scattered light from bright stars)
        # from being selected for sky objects in the calculation,
        # so do a preliminary detection pass without either the local or wide
        # temporary background subtraction; the DETECTED pixels will mark
        # the area to ignore.
        originalMask = maskedImage.mask.array.copy()
        try:
            self.clearMask(exposure.mask)
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
            tweakDetResults = self.applyThreshold(convolveResults.middle, maskedImage.getBBox(), self.config.prelimThresholdFactor)
            self.finalizeFootprints(maskedImage.mask, tweakDetResults, convolveResults.sigma, self.config.prelimThresholdFactor)
            # Calculate the proper threshold
            threshResults = self.calculateThreshold(exposure, seed, sigma=convolveResults.sigma)
            bgLevel = threshResults.additive
            factor = threshResults.multiplicative
            self.log.info("Modifying configured detection threshold by factor %f to %f",
                          factor, factor*self.config.thresholdValue)

        finally:
            maskedImage.mask.array[:] = originalMask

        background = lsst.afw.math.BackgroundList()
        if self.config.doBackgroundTweak:
            self.tweakBackground(exposure, bgLevel, background)

        if clearMask:
            self.clearMask(maskedImage.mask)
        else:
            oldDetected = maskedImage.mask.array & maskedImage.mask.getPlaneBitMask(["DETECTED",
                                                                                     "DETECTED_NEGATIVE"])

        with self.tempWideBackgroundContext(exposure):
            # Could potentially smooth with a wider kernel than the PSF in order to better pick up the
            # wings of stars and galaxies, but for now sticking with the PSF as that's more simple.
            convolveResults = self.convolveImage(maskedImage, psf, doSmooth=doSmooth)
            middle = convolveResults.middle
            sigma = convolveResults.sigma
            results = self.applyThreshold(middle, maskedImage.getBBox(), factor)
            results.background = background
            if self.config.doTempLocalBackground:
                self.applyTempLocalBackground(exposure, middle, results)
            self.finalizeFootprints(maskedImage.mask, results, sigma, factor)

            self.clearUnwantedResults(maskedImage.mask, results)

        if self.config.reEstimateBackground:
            self.reEstimateBackground(maskedImage, results.background)

        self.display(exposure, results, middle)

        return results

    def tweakBackground(self, exposure, bgLevel, bgList=None):
        """Modify the background by a constant value

        Parameters
        ----------
        exposure : `lsst.afw.image.Exposure`
            Exposure for which to tweak background.
        bgLevel : `float`
            Background level to remove
        bgList : `lsst.afw.math.BackgroundList`, optional
            List of backgrounds to append to.

        Returns
        -------
        bg : `lsst.afw.math.BackgroundMI`
            Constant background model.
        """
        self.log.info("Tweaking background by %f to match sky photometry", bgLevel)
        exposure.image -= bgLevel
        bgStats = lsst.afw.image.MaskedImageF(1, 1)
        bgStats.set(bgLevel, 0, bgLevel)
        bg = lsst.afw.math.BackgroundMI(exposure.getBBox(), bgStats)
        bgData = (bg, lsst.afw.math.Interpolate.LINEAR, lsst.afw.math.REDUCE_INTERP_ORDER,
                  lsst.afw.math.ApproximateControl.UNKNOWN, 0, 0, False)
        if bgList is not None:
            bgList.append(bgData)
        return bg
