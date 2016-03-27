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

"""Test lsst.coadd.utils.setMaskBits
"""
import unittest

import lsst.afw.image as afwImage
import lsst.utils.tests as utilsTests
import lsst.pex.logging as pexLog
import lsst.coadd.utils as coaddUtils

Verbosity = 0

pexLog.Trace_setVerbosity("lsst.coadd.utils", Verbosity)


def countBits(val):
    nBits = 0
    while val > 0:
        if val & 1 != 0:
            nBits += 1
        val >>= 1
    return nBits

MaxBitMask = (2**afwImage.MaskU_getNumPlanesMax() - 1)

MaskPlaneNameIDDict = dict(afwImage.MaskU().getMaskPlaneDict())


class AddToCoaddTestCase(unittest.TestCase):
    """A test case for setMaskBits
    """

    def testBits(self):
        """Test that the bits set are correct"""

        fullPlaneNameList = MaskPlaneNameIDDict.keys()
        totNumBits = countBits(MaxBitMask)
        for i in range(len(fullPlaneNameList)):
            numPlanes = i + 1
            setPlaneNameList = fullPlaneNameList[0:numPlanes]

            bitMask = coaddUtils.makeBitMask(setPlaneNameList)
            self.assertEqual(countBits(bitMask), numPlanes)
            for planeName, planeId in MaskPlaneNameIDDict.iteritems():
                self.assertEqual((bitMask & (1 << planeId)) > 0, planeName in setPlaneNameList)

            invBitMask = coaddUtils.makeBitMask(setPlaneNameList, doInvert=True)
            self.assertEqual(countBits(invBitMask), totNumBits - numPlanes)
            for planeName, planeId in MaskPlaneNameIDDict.iteritems():
                self.assertEqual((invBitMask & (1 << planeId)) > 0, planeName not in setPlaneNameList)


def suite():
    """Return a suite containing all the test cases in this module.
    """
    utilsTests.init()

    suites = []
    suites += unittest.makeSuite(AddToCoaddTestCase)
    suites += unittest.makeSuite(utilsTests.MemoryTestCase)

    return unittest.TestSuite(suites)


def run(shouldExit=False):
    """Run the tests"""
    utilsTests.run(suite(), shouldExit)

if __name__ == "__main__":
    run(True)
