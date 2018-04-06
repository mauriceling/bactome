'''
Bactome package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from distutils.core import setup
import py2exe

# Make Windows executable:python setup.py py2exe

setup(name='oliver',
      version='1.0',
      description='OLIgonucleotide Variable Expression Ranker (OLIVER)',
      long_description='''
OLIgonucleotide Variable Expression Ranker (OLIVER): A set of tools for 
identifying suitable reference (invariant) genes from large transcriptome 
datasets. 

Part of Bactome project (https://github.com/mauriceling/bactome)

References:
- Chan, OYW, Keng, BMH, Ling, MHT. 2014. Correlation and Variation Based 
Method for Reference Genes Identification from Large Datasets. Electronic 
Physician 6(1): 719-727.

Authors:
- Oliver YW Chan, Raffles Institution, Singapore
- Bryan MH Keng, Raffles Institution, Singapore
- Maurice HT Ling (mauriceling@acm.org), Department of Zoology, The 
University of Melbourne, Australia

Copyright: Copyright (c) 2013-2014, Maurice Ling (on behalf of all 
authors)
''',
      author='Maurice HT Ling',
      author_email='mauriceling@acm.org',
      license = 'General Public License version 3',
      console=['oliver.py'])
