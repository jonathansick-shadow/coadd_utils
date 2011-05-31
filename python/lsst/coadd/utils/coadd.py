# 
# LSST Data Management System
# Copyright 2008, 2009, 2010 LSST Corporation.
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
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath
import makeBitMask
import utilsLib

__all__ = ["Coadd"]

class Coadd(object):
    """Coadd by weighted addition
    
    This class may be subclassed to implement other coadd techniques.
    Typically this is done by overriding addExposure.
    """
    def __init__(self, bbox, wcs, badMaskPlanes, coaddZeroPoint, logName="coadd.utils.Coadd"):
        """Create a coadd
        
        @param bbox: bounding box of coadd Exposure with respect to parent (lsst.afw.geom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param badMaskPlanes: mask planes to pay attention to when rejecting masked pixels.
            Specify as a collection of names.
            badMaskPlanes should always include "EDGE".
        @param coaddZeroPoint: photometric zero point of coadd (mag)
        @param logName: name by which messages are logged
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), logName)

        self._badPixelMask = makeBitMask.makeBitMask(badMaskPlanes)
        self._coaddZeroPoint = float(coaddZeroPoint)

        self._bbox = bbox
        self._wcs = wcs
        self._coadd = afwImage.ExposureF(bbox, wcs)

        coddFluxMag0 = 10**(0.4 * coaddZeroPoint)
        calib = afwImage.Calib()
        calib.setFluxMag0(coddFluxMag0)
        if abs(calib.getMagnitude(1.0) - self._coaddZeroPoint) > 1.0e-4:
            raise RuntimeError("Bug: calib.getMagnitude(1.0) = %0.4f != %0.4f = coaddZeroPoint" % \
                (calib.getMagnitude(1.0), self._coaddZeroPoint))
        self._coadd.setCalib(calib)

        self._weightMap = afwImage.ImageF(bbox, afwImage.PARENT)
        
        self._statsControl = afwMath.StatisticsControl()
        self._statsControl.setNumSigmaClip(3.0)
        self._statsControl.setNumIter(2)
        self._statsControl.setAndMask(self._badPixelMask)
    
    @classmethod
    def fromPolicy(cls, bbox, wcs, policy, logName="coadd.utils.Coadd"):
        """Create a coadd
        
        @param bbox: bounding box of coadd Exposure with respect to parent (lsst.afw.geom.Box2I):
            coadd dimensions = bbox.getDimensions(); xy0 = bbox.getMin()
        @param wcs: WCS of coadd exposure (lsst.afw.math.Wcs)
        @param policy: coadd policy; see policy/CoaddPolicyDictionary.paf
        @param logName: name by which messages are logged
        """
        return cls(
            bbox = bbox,
            wcs = wcs,
            badMaskPlanes = policy.getArray("badMaskPlanes"),
            coaddZeroPoint = policy.get("coaddZeroPoint"),
            logName = logName,
        )

    def addExposure(self, exposure, weightFactor=1.0):
        """Add an Exposure to the coadd
        
        @param exposure: Exposure to add to coadd; this must be:
            - background-subtracted
            - warped to match the coadd
            - photometrically calibrated (have a Calib object with nonzero fluxMag0)
        @param weightFactor: extra weight factor for this exposure

        @return
        - overlapBBox: region of overlap between exposure and coadd in parent coordinates (afwGeom.Box2I)
        - weight: weight with which exposure was added to coadd; weight = weightFactor / clipped mean variance
        
        Subclasses may override to preprocess the exposure or change the way it is added to the coadd.
        """
        # normalize a deep copy of the masked image so flux is 1 at the coadd zero point;
        # use a deep copy to avoid altering the input exposure
        fluxAtZeropoint = exposure.getCalib().getFlux(self._coaddZeroPoint)
        scaleFac = 1.0 / fluxAtZeropoint
        maskedImage = exposure.getMaskedImage()
        maskedImage = maskedImage.Factory(maskedImage, True)
        maskedImage *= scaleFac
        
        # compute the weight
        statObj = afwMath.makeStatistics(maskedImage.getVariance(), maskedImage.getMask(),
            afwMath.MEANCLIP, self._statsControl)
        meanVar, meanVarErr = statObj.getResult(afwMath.MEANCLIP);
        weight = weightFactor / float(meanVar)

        self._log.log(pexLog.Log.INFO, "add masked image to coadd; scaled by %0.3g; weight=%0.3g" % \
            (scaleFac, weight,))

        overlapBBox = utilsLib.addToCoadd(self._coadd.getMaskedImage(), self._weightMap,
            maskedImage, self._badPixelMask, weight)

        return overlapBBox, weight

    def getCoadd(self):
        """Get the coadd exposure for all exposures you have coadded so far
        """
        # make a deep copy so I can scale it
        coaddMaskedImage = self._coadd.getMaskedImage()
        scaledMaskedImage = coaddMaskedImage.Factory(coaddMaskedImage, True)

        # set the edge pixels
        utilsLib.setCoaddEdgeBits(scaledMaskedImage.getMask(), self._weightMap)
        
        # scale non-edge pixels by weight map
        scaledMaskedImage /= self._weightMap
        
        scaledExposure = afwImage.makeExposure(scaledMaskedImage, self._wcs)
        scaledExposure.setCalib(self._coadd.getCalib())
        return scaledExposure
    
    def getCoaddZeroPoint(self):
        """Return the coadd photometric zero point.
        
        getCoaddFluxMag0 gives the same information in different units.
        """
        return self._coaddZeroPoint
    
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
        
