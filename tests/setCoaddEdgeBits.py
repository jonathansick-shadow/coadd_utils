#!/usr/bin/env python
"""
Test lsst.coadd.utils.setCoaddEdgeBits
"""

import os
import math
import pdb # we may want to say pdb.set_trace()
import unittest

import numpy

import eups
import lsst.afw.image as afwImage
import lsst.afw.image.testUtils as imTestUtils
import lsst.afw.math as afwMath
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.coadd.utils as coaddUtils

Verbosity = 0 # increase to see trace
pexLog.Trace_setVerbosity("lsst.coadd.utils", Verbosity)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


class SetCoaddEdgeBitsTestCase(unittest.TestCase):
    """
    A test case for setCoaddEdgeBits
    """
    def testRandomMap(self):
        """Test setCoaddEdgeBits using a random depth map
        """
        imShape = (50, 50)
        coaddMask = afwImage.MaskU(imShape[0], imShape[1], 0)
        

        numpy.random.seed(12345)
        depthMapArray = numpy.random.randint(0, 3, list(imShape))
        depthMap = imTestUtils.imageFromArray(depthMapArray, afwImage.ImageU)
        
        refCoaddMaskArray = imTestUtils.arrayFromMask(coaddMask)
        edgeMask = afwImage.MaskU.getPlaneBitMask("EDGE")
        refCoaddMaskArray |= numpy.where(depthMapArray > 0, 0, edgeMask)
        
        coaddUtils.setCoaddEdgeBits(coaddMask, depthMap)
        coaddMaskArray = imTestUtils.arrayFromMask(coaddMask)
        if numpy.any(refCoaddMaskArray != coaddMaskArray):
            errMsgList = (
                "Coadd mask does not match reference=%s:" % (badPixelMask,),
                "computed=  %s" % (coaddMaskArray,),
                "reference= %s" % (refCoaddMaskArray,),
            )
            errMsg = "\n".join(errMsgList)
            self.fail(errMsg)

#-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

def suite():
    """
    Returns a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(SetCoaddEdgeBitsTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)

if __name__ == "__main__":
    utilsTests.run(suite())
