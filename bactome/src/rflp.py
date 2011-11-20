"""
Simulate and analyze genetic distances from RFLP results

Copyright 2009 Maurice Ling <mauriceling@acm.org> on behalf of
Maurice Ling, Chin-How Lee, Jack Si-Hao Oon, Kun-Cheng Lee.

Licence: Python Software Foundation License version 2

@see: Lee, CH, Lee, KC, Oon, JSH, Ling, MHT. 2010. Bactome, I: 
Python in DNA Fingerprinting. In: Peer-Reviewed Articles from 
PyCon Asia-Pacific 2010. The Python Papers 5(3): 6.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import math

def restriction_digest(seq, enzyme, max_band=23130, min_band=2000, 
                       linear=False, p='yes'):
    """
    Performs restriction endonuclease digestion on a sequence and group the
    resulting fragments into 3 groups so simulate different agarose gel
    electrophoresis:
    1. fragment length more than the maximum size
    2. fragment length between the maximum and minimum size
    3. fragment length less than the minimum size
    
    Parameters:
        seq = DNA sequence for restriction endonuclease digestion
        enzyme = Restriction endonuclease object from Bio.Restriction
            package
        max_band = size of maximum band in basepairs. Default = 23130
        min_band = size of minimum band in basepairs. Default = 2000
        linear = flag to define if DNA sequence is linear.
            Default = False (DNA is circular)
        p = flag to determine if the data is to be printed. Default = yes
        
    Result:
        (Number of fragments after digestion, 
        List of fragments with molecular size above max_band, 
        List of fragments with molecular size between max_band and min_band, 
        List of fragments with molecular size below min_band)
    """
    digest = enzyme.search(seq, linear=linear)
    digest.sort()
    fragment = [digest[x+1] - digest[x]
                for x in range(len(digest) - 1)]
    fragment.sort()
    ogel = [x for x in fragment if x > max_band]
    gel = [x for x in fragment if x <= max_band and x >= min_band]
    ugel = [x for x in fragment if x < min_band]
    ogel.sort()
    gel.sort()
    ugel.sort()
    if p == 'yes':
        print 'Enzyme: ' + str(enzyme)
        print 'Restriction site: ' + enzyme.site
        print 'Number of fragments: ' + str(len(fragment))
        print 'Number of fragments (x > ' + str(max_band) + '): ' + \
            str(len(ogel))
        print 'Number of fragments (' + str(max_band) + ' < x < ' + \
            str(min_band) + '): ' + str(len(gel))
        print 'Number of fragments (x < ' + str(min_band) + '): ' + \
            str(len(ugel))
        print
    return (len(fragment), ogel, gel, ugel)

def restriction_supplier(seq, max_band=23130, min_band=2000,
                         suppliers='ACEGFIHKJMONQPSRUVX',
                         linear=False, p='yes'):
    """
    Performs restriction endonuclease analysis by batch, based on supplier.
    
    Parameters:
        seq = DNA sequence for restriction endonuclease digestion
        max_band = size of maximum band in basepairs. Default = 23130
        mmin_band = size of minimum band in basepairs. Default = 2000
        suppliers = restriction enzyme supplier. Default = ACEGFIHKJMONQPSRUVX
            where
                A = Amersham Pharmacia Biotech
                C = Minotech Biotechnology
                E = Stratagene
                G = Qbiogene
                F = Fermentas AB
                I = SibEnzyme Ltd.
                H = American Allied Biochemical, Inc.
                K = Takara Shuzo Co. Ltd.
                J = Nippon Gene Co., Ltd.
                M = Roche Applied Science
                O = Toyobo Biochemicals
                N = New England Biolabs
                Q = CHIMERx
                P = Megabase Research Products
                S = Sigma Chemical Corporation
                R = Promega Corporation
                U = Bangalore Genei
                V = MRC-Holland
                X = EURx Ltd.
        linear = flag to define if DNA sequence is linear.
            Default = False (DNA is circular)
        p = flag to determine if the data is to be printed. Default = yes

    Returns:
    {Restriction endonuclease :
        (Total number of fragments after digestion,
        Number of fragments with molecular size above max_band, 
        Number of fragments with molecular size between max_band and min_band, 
        Number of fragments with molecular size below min_band)}
    """
    from Bio.Restriction import RestrictionBatch
    count = 0
    result = {}
    for enzyme in RestrictionBatch(first=[], 
                                   suppliers=[x.upper() for x in suppliers]):
        try:
            digest = restriction_digest(seq, enzyme, max_band, min_band,
                                        linear, p)
        except MemoryError:
            print 'Memory Error during ' + str(enzyme) + ' digestion'
        result[str(enzyme)] = (digest[0], len(digest[1]),
                               len(digest[2]), len(digest[3]))
        count = count + 1
        if p != 'yes':
            if count % 10 == 0:
                print str(count) + ' restriction endonuclease processed'
    return result
    
def setCompare(original, test, absent):
    """
    Used for processing set-based (unordered or nominal) distance of 
    categorical data.
    
    Parameters:
        original: list of original data
        test: list of data to test against original
        absent: indicator to define absent data
    """
    original_only = float(len([x for x in original
                               if x not in test]))
    test_only = float(len([x for x in test if x not in original]))
    both = float(len([x for x in original if x in test]))
    return (original_only, test_only, both)

def listCompare(original, test, absent):
    """
    Used for processing list-based (ordered or ordinal) distance of 
    categorical data.
    
    Parameters:
        original: list of original data
        test: list of data to test against original
        absent: indicator to define absent data
    """
    original = list(original)
    test = list(test)
    original_only = 0.0
    test_only = 0.0
    both = 0.0
    for i in range(len(original)):
        if original[i] == absent and test[i] == absent: pass
        elif original[i] == test[i]: both = both + 1
        elif original[i] <> absent and test[i] == absent:
            original_only = original_only + 1
        elif original[i] == absent and test[i] <> absent: 
            test_only = test_only + 1
        else: pass
    return (original_only, test_only, both)

def Nei_Li(original, test, absent=0, type='Set'):
    """
    Nei and Li Distance is distance measure for nominal or ordinal data.
    
    Given 2 lists (original and test), calculates the Nei and Li Distance 
    based on the formula,
    
    1 - [2 x (number of regions where both species are present)/
    [(2 x (number of regions where both species are present)) + 
    (number of regions where only one species is present)]]
        
    @see: Nei M, Li WH (1979) Mathematical models for studying genetic
    variation in terms of restriction endonucleases.
    Proc Natl Acad Sci USA 76:5269-5273
    
    @param original: list of original data
    @param test: list of data to test against original
    @param absent: user-defined identifier for absent of region, default = 0
    @param type: {Set | List}, define whether use Set comparison (unordered) or
        list comparison (ordered), default = Set
    """
    if type == 'Set':
        (original_only, test_only, both) = setCompare(original, test, absent)
    else:
        (original_only, test_only, both) = listCompare(original, test, absent)
    return 1-((2*both)/((2*both)+original_only+test_only))