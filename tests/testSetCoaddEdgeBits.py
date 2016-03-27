#!/usr/bin/env python

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

"""Test lsst.coadd.utils.setCoaddEdgeBits
"""
import unittest

import numpy

import lsst.afw.geom as afwGeom
import lsst.afw.image as afwImage
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.coadd.utils as coaddUtils

Verbosity = 0  # increase to see trace
pexLog.Trace_setVerbosity("lsst.coadd.utils", Verbosity)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class SetCoaddEdgeBitsTestCase(unittest.TestCase):
    """A test case for setCoaddEdgeBits
    """

    def testRandomMap(self):
        """Test setCoaddEdgeBits using a random depth map
        """
        imDim = afwGeom.Extent2I(50, 55)
        coaddMask = afwImage.MaskU(imDim)

        numpy.random.seed(12345)
        depthMapArray = numpy.random.randint(0, 3, list((imDim[1], imDim[0]))).astype(numpy.uint16)
        depthMap = afwImage.makeImageFromArray(depthMapArray)

        refCoaddMaskArray = coaddMask.getArray()
        edgeMask = afwImage.MaskU.getPlaneBitMask("NO_DATA")
        refCoaddMaskArray |= numpy.array(numpy.where(depthMapArray > 0, 0, edgeMask),
                                         dtype=refCoaddMaskArray.dtype)

        coaddUtils.setCoaddEdgeBits(coaddMask, depthMap)
        coaddMaskArray = coaddMask.getArray()
        if numpy.any(refCoaddMaskArray != coaddMaskArray):
            errMsgList = (
                "Coadd mask does not match reference:",
                "computed=  %s" % (coaddMaskArray,),
                "reference= %s" % (refCoaddMaskArray,),
            )
            errMsg = "\n".join(errMsgList)
            self.fail(errMsg)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SetCoaddEdgeBitsTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
