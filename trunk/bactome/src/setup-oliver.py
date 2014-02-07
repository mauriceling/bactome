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

References:
- Chan, OYW, Keng, BMH, Ling, MHT. 2014. Correlation and Variation Based 
Method for Reference Genes Identification from Large Datasets. Electronic 
Physician 6(1): 719-727.
''',
      author='Maurice HT Ling',
      author_email='mauriceling@acm.org',
      license = 'General Public License version 3',
      console=['oliver.py'])
