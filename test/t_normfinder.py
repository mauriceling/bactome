"""
Test script for src/normfinder.py

This script only works with Python version < 3.0 due to print 
statements. 

Date created: 23rd November 2011

Licence: Python Software Foundation License version 2
"""

import sys
import os
import unittest

sys.path.append(os.path.dirname(os.getcwd()))

import gnormfinder as g

testdata = open('genormfinder_790452.csv', 'r').readlines()
testdata = [x[:-1].split(',') for x in testdata]

class testGO(unittest.TestCase):
    def testnormfinder_1(self):
        result = g.normfinder(testdata, True)
        self.assertTrue(len(result), 932)
        
    def testnormfinder_2(self):
        result = g.normfinder(testdata, True)
        self.assertTrue(result['2'][0], 'SPCS1')
        self.assertAlmostEqual(result['2'][1], 0.32639493)
        self.assertTrue(result['66'][0], 'HADHB')
        self.assertAlmostEqual(result['66'][1], 0.35463989)
        