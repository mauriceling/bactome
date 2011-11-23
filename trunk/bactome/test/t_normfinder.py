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

sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'src'))

import gnormfinder as g

testdata = open('genormfinder_790452.csv', 'r').readlines()
testdata = [x[:-1].split(',') for x in testdata]

class testGO(unittest.TestCase):
    def testnormfinder_1(self):
        result = g.normfinder(testdata, True)
        self.assertTrue(len(result), 932)
        
    def testnormfinder_2(self):
        result = g.normfinder(testdata, True)
        self.assertTrue(result['SPCS1'], [0.32639493])
        self.assertTrue(result['HADHB'], [0.35463989])
        
    def testnormfinder_3(self):
        result = g.normfinder(testdata, True)
        testvalues = [0.851, 0.852, 0.854, 0.860, 0.868, 0.877, 
                      0.882, 0.885, 0.888]
        print 'Generated results for ACTG1: ' + str(result['ACTG1'])
        print 'Test values for ACTG1: ' + str(testvalues)
        self.assertTrue(len(result['ACTG1']), 9)
