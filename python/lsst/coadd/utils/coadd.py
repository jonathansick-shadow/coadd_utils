#
# LSST Data Management System
# Copyright 2008, 2009, 2010, 2011, 2012 LSST Corporation.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#
import math

import lsst.pex.config as pexConfig
from lsst.pex.logging import Log
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
from .utilsLib import addToCoadd, setCoaddEdgeBits

__all__ = ["Coadd"]


class CoaddConfig(pexConfig.Config):
    """Config for Coadd
    """
    badMaskPlanes = pexConfig.ListField(
        dtype = str,
        doc = "mask planes that, if set, the associated pixel should not be included in the coadd",
        default = ("NO_DATA", "SAT"),
    )


class Coadd(object):
    """Coadd by weighted addition

    This class may be subclassed to implement other coadd techniques.
    Typically this is done by overriding addExposure.
    """
    ConfigClass = CoaddConfig

    def __init__(self, bbox, wcs, badMaskPlanes, logName="coadd.utils.Coadd"):
        """Create a coadd

        @param[in] bbox: bounding box of coadd Exposure with respect to parent (lsst.afw.geom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param[in] wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param[in] badMaskPlanes: mask planes to pay attention to when rejecting masked pixels.
            Specify as a collection of names.
            badMaskPlanes should always include "NO_DATA".
        @param[in] logName: name by which messages are logged
        """
        self._log = Log(Log.getDefaultLog(), logName)
        self._bbox = bbox
        self._wcs = wcs
        self._badPixelMask = afwImage.MaskU.getPlaneBitMask(badMaskPlanes)
        self._coadd = afwImage.ExposureF(bbox, wcs)
        self._weightMap = afwImage.ImageF(bbox)
        self._filterDict = dict()  # dict of filter name: filter object for all filters seen so far

        self._statsControl = afwMath.StatisticsControl()
        self._statsControl.setNumSigmaClip(3.0)
        self._statsControl.setNumIter(2)
        self._statsControl.setAndMask(self._badPixelMask)

    @classmethod
    def fromConfig(cls, bbox, wcs, config, logName="coadd.utils.Coadd"):
        """Create a coadd

        @param[in] bbox: bounding box of coadd Exposure with respect to parent (lsst.afw.geom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param[in] wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param[in] config: coadd config; an instance of CoaddConfig
        @param[in] logName: name by which messages are logged
        """
        return cls(
            bbox = bbox,
            wcs = wcs,
            badMaskPlanes = config.badMaskPlanes,
            logName = logName,
        )

    def addExposure(self, exposure, weightFactor=1.0):
        """Add an Exposure to the coadd

        @param[in] exposure: Exposure to add to coadd; this should be:
            - background-subtracted or background-matched to the other images being coadded
            - psf-matched to the desired PSF model (optional)
            - warped to match the coadd
            - photometrically scaled to the desired flux magnitude
        @param[in] weightFactor: extra weight factor for this exposure

        @return
        - overlapBBox: region of overlap between exposure and coadd in parent coordinates (afwGeom.Box2I)
        - weight: weight with which exposure was added to coadd; weight = weightFactor / clipped mean variance

        Subclasses may override to preprocess the exposure or change the way it is added to the coadd.
        """
        maskedImage = exposure.getMaskedImage()

        # compute the weight
        statObj = afwMath.makeStatistics(maskedImage.getVariance(), maskedImage.getMask(),
                                         afwMath.MEANCLIP, self._statsControl)
        meanVar = statObj.getResult(afwMath.MEANCLIP)[0]
        weight = weightFactor / float(meanVar)
        if math.isnan(weight):
            raise RuntimeError("Weight is NaN (weightFactor=%s; mean variance=%s)" % (weightFactor, meanVar))

        # save filter info
        filter = exposure.getFilter()
        self._filterDict.setdefault(filter.getName(), filter)

        self._log.log(Log.INFO, "Add exposure to coadd with weight=%0.3g" % (weight,))

        overlapBBox = addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
                                 maskedImage, self._badPixelMask, weight)

        return overlapBBox, weight

    def getCoadd(self):
        """Get the coadd exposure for all exposures you have coadded so far

        If all exposures in this coadd have the same-named filter then that filter is set in the coadd.
        Otherwise the coadd will have the default unknown filter.

        @warning: the Calib is not be set.
        """
        # make a deep copy so I can scale it
        coaddMaskedImage = self._coadd.getMaskedImage()
        scaledMaskedImage = coaddMaskedImage.Factory(coaddMaskedImage, True)

        # set the edge pixels
        setCoaddEdgeBits(scaledMaskedImage.getMask(), self._weightMap)

        # scale non-edge pixels by weight map
        scaledMaskedImage /= self._weightMap

        scaledExposure = afwImage.makeExposure(scaledMaskedImage, self._wcs)
        if len(self._filterDict) == 1:
            scaledExposure.setFilter(self._filterDict.values()[0])
        return scaledExposure

    def getFilters(self):
        """Return a collection of all the filters seen so far in in addExposure
        """
        return self._filterDict.values()

    def getBadPixelMask(self):
        """Return the bad pixel mask
        """
        return self._badPixelMask

    def getBBox(self):
        """Return the bounding box of the coadd
        """
        return self._bbox

    def getWcs(self):
        """Return the wcs of the coadd
        """
        return self._wcs

    def getWeightMap(self):
        """Return the weight map for all exposures you have coadded so far

        The weight map is a float Image of the same dimensions as the coadd; the value of each pixel
        is the sum of the weights of all exposures that contributed to that pixel.
        """
        return self._weightMap
