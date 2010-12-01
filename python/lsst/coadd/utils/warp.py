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
import sys
import lsst.pex.logging as pexLog
import lsst.afw.image as afwImage
import lsst.afw.math as afwMath

__all__ = ["Warp"]

class Warp(object):
    """Warp exposures
    """
    def __init__(self, warpingKernelName, logName="coadd.utils.WarpExposure"):
        """Create a Warp
        
        Inputs:
        - warpingKernelName: name of warping kernel
        - logName: name by which messages are logged
        """
        self._log = pexLog.Log(pexLog.Log.getDefaultLog(), logName)
        self._warpingKernel = afwMath.makeWarpingKernel(warpingKernelName)

    def warpExposure(self, dimensions, wcs, exposure):
        """Warp an exposure
        
        Inputs:
        - dimensions: dimensions of warped exposure; must be the type of object returned by
            exposure.getMaskedImage().getDimensions() (presently std::pair<int, int>)
        - wcs: WCS of warped exposure
        - exposure: Exposure to warp
            
        Returns:
        - warpedExposure: warped exposure
        """
        self._log.log(pexLog.Log.INFO, "warp exposure")
        maskedImage = exposure.getMaskedImage()
        blankMaskedImage = maskedImage.Factory(dimensions)
        warpedExposure = afwImage.makeExposure(blankMaskedImage, wcs)
        afwMath.warpExposure(warpedExposure, exposure, self._warpingKernel)
        return warpedExposure