"""
Test script for src/rflp.py

Copyright 2009 Maurice Ling <mauriceling@acm.org> on behalf of
Maurice Ling, Chin-How Lee, Jack Si-Hao Oon, Kun-Cheng Lee.

Licence: Python Software Foundation License version 2

@see: Lee, CH, Lee, KC, Oon, JSH, Ling, MHT. 2010. Bactome, I: 
Python in DNA Fingerprinting. In: Peer-Reviewed Articles from 
PyCon Asia-Pacific 2010. The Python Papers 5(3): 6.
"""

import sys
import os
import unittest

sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'src'))

from rflp import Nei_Li

profile = [[3.8, 4.0, 4.1, 4.2, 4.4, 5.0, 5.4, 5.9, 6.0, 6.4, 6.5, 6.8, 7.2, 7.3],
           [4.0, 4.1, 4.2, 5.0, 5.4, 5.9, 6.0, 6.4, 6.8, 7.2, 7.3],
           [3.8, 4.2, 5.0, 5.4, 5.9, 6.0, 6.4, 6.5, 6.8, 7.2, 7.3],
           [3.8, 4.0, 4.2, 4.4, 5.0, 5.4, 5.9, 6.0, 6.4, 6.5, 6.8, 7.3],
           [4.2, 5.4, 5.9, 6.0, 6.4, 6.8, 7.3],
           [3.8, 4.2, 4.4, 5.0, 5.4, 5.9, 6.0, 6.4, 6.8, 7.2, 7.3],
           [3.8, 4.1, 4.2, 4.4, 5.0, 5.4, 5.9, 6.0,6.4, 6.5, 6.8, 7.3],
           [3.8, 4.2, 5.0, 5.4, 5.9, 6.0,6.4, 6.5, 6.8, 7.3]]

distance = [(p1, p2, Nei_Li(profile[p1], profile[p2]))
            for p1 in range(len(profile))
            for p2 in range(len(profile))]

print(distance)
