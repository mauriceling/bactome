"""
Generate potential primers for RAPD/RFLP from reference genomic sequence

Copyright 2009 Maurice Ling <mauriceling@acm.org> on behalf of
Maurice Ling, Chin-How Lee, Jack Si-Hao Oon, Kun-Cheng Lee.

Licence: Python Software Foundation License version 2

@see: Lee, CH, Lee, KC, Oon, JSH, Ling, MHT. 2010. Bactome, I: 
Python in DNA Fingerprinting. In: Peer-Reviewed Articles from 
PyCon Asia-Pacific 2010. The Python Papers 5(3): 6.
"""

import anydbm
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
import math

from marshaldbm import marshaldbm

def fasta_seq(fasta='ATCC8739.fasta'):
    """
    Open fasta format file and returns the sequence.
    """
    f = SeqIO.parse(open(fasta, 'rU'), 'fasta').next()
    return f.seq

def generate_fragments(seq, tablefile, fragment_size=15):
    """
    Slice an entire sequence into fragments of equal lengths
    
    Parameters:
        seq = sequence to be sliced (output from fasta_seq function)
        tablefile = name of file to store the sliced sequence
        fragment_size = length of each sequence to be sliced into. This will
            be the maximum primer/probe length.
            Default value is 15 bases.
    """
    genome_fragment = anydbm.open(tablefile, 'c')
    last_start= len(seq) - fragment_size
    count = 0
    for start in range(0, last_start, fragment_size):
        genome_fragment[str(count)] = seq[start:15+start].data
        count = count + 1
        if (count % 1000) == 0:
            print str(count) + ' fragments inserted into ' + tablefile
    print str(count) + ' fragments inserted into ' + tablefile
    genome_fragment.close()

def LCS(seq, test):
    """
    Finds the longest common substring between the 2 inputs (seq and test).
    
    Adapted from http://en.wikibooks.org/wiki/Algorithm_Implementation/
    Strings/Longest_common_substring
    """
    m = len(seq)
    n = len(test)
    L = [[0] * (n+1) for i in xrange(m+1)]
    LCS = set()
    longest = 0
    for i in xrange(m):
        for j in xrange(n):
            if seq[i] == test[j]:
                v = L[i][j] + 1
                L[i+1][j+1] = v
                if v > longest:
                    longest = v
                    LCS = set()
                if v == longest:
                    LCS.add(seq[i-v+1:i+1])
    return LCS

def generate_all_LCS(fragfile, LCSfile, 
                     min_primer=7, max_interval=198, min_interval=6):
    """
    Generates primers where forward primer sequence = reverse primer sequence
    given the size of amplicon and minimum length of primers.
    
    Parameters:
        fragfile = file of sliced sequence (tablefile parameter of 
            generate_fragments function)
        LCSfile = name of file to store least common substrings' (primers) data 
        min_primer = smallest acceptable length for primers
        max_interval = number of maximum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 198 as 3000bp amplicon is the maximum normal
            Taq polymerase can amplify. 198 intervals of 15 bases per 
            fragment + 2 flanking fragment of maximum of 15 bases (primer
            length) gives a total of 200 fragments of 15 bases = 3000 bases.
        min_interval = number of minimum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 6 as 120bp amplicon is the minimum length visible
            on 1.5% agarose gel that is also suitable for 3000 bp. 6 intervals 
            of 15 bases per fragment + 2 flanking fragment of maximum of 15 
            bases (primer length) gives a total of 8 fragments of 15 bases = 
            120 bases.
    """
    genome_fragment = anydbm.open(fragfile, 'c')
    num_of_frag = len(genome_fragment)
    num_analyzed = 0
    result = {}
    fragment_size = len(genome_fragment['0'])
    for size in range(max_interval, min_interval, -1):
        com_seq = anydbm.open(LCSfile, 'c')
        print 'Analysing for ' + str(fragment_size * size) + 'bp amplicons'
        num_com = 0
        for frag in range(0, num_of_frag - size):
            p1 = genome_fragment[str(frag)]
            p2 = genome_fragment[str(frag+size)]
            p2 = Seq(p2, generic_dna).reverse_complement().data
            try: common_seq = LCS(p1, p2).pop()
            except : common_seq = ''
            num_analyzed = num_analyzed + 1
            if len(common_seq) > min_primer:
                com_seq['|'.join([str(frag), str(frag+size),
                                  str(size)])] = common_seq
                num_com = num_com + 1
            if (num_analyzed % 10000) == 0:
                print str(num_analyzed) + ' pairs analyzed'
                print '        ' + str(num_com) + ' LCS more than ' + \
                      str(min_primer) + 'bp for ' + \
                      str(15 * size) + 'bp amplicons'
        com_seq.close()
        result[size] = str(num_com)
    print str(num_analyzed) + ' pairs analyzed'
    print '        ' + str(num_com) + ' LCS more than ' + str(min_primer) + \
        'bp for ' + str(15 * size) + 'bp amplicons'
    genome_fragment.close()
    return result

def get_LCS_by_amplicon(LCSfile, amplicon_size):
    """
    Extract all the primers and amplicon positions from the LCSfile given the
    amplicon length.
    
    Parameters:
        LCSfile = file of least common substrings' (primers) data
        amplicon_size = size of amplified product in terms of number of
            fragments intervening 2 primers. For example, if the fragment 
            size is 15 and the amplicon size is 60, we are then looking
            for primers (forward primer sequence = reverse primer
            sequence) that gives amplicons of 900 bp (60*15) excluding
            the flanking primers.
    
    Returns:
        (List of primers,
        List of amplicons' data)
        where
        1. each data element in "List of amplicons' data" is in the format of
           '<start fragment position>|<end fragment position>|<amplicon size>'
            <start fragment position> and <end fragment position> multiplied 
            by the fragment size will give the corresponding position on the 
            fasta sequence. For example, '89572|89632|60' represents the 
            amplicon from 1343580 bp (889572*15) to 1344480 bp (89632*15) of 
            the fasta sequence, corresponding to an amplicon size of 900bp 
            (60*15).
        2. the number of primers = number of primers's data. This means that in
            primers = get_LCS_by_amplicon('primer.table', 60)
            primers[0][10] sequence will amplify primers[1][10]
    """
    com_seq = anydbm.open(LCSfile, 'c')
    keys = com_seq.keys()
    return ([com_seq[x]
             for x in keys
                 if x.split('|')[2] == str(amplicon_size)],
            [x
             for x in keys
                 if x.split('|')[2] == str(amplicon_size)])
    
def LCS_uniqueness(LCSfile, amplicon_size):
    """
    Extract all the primers and amplicon positions from the LCSfile given the
    amplicon length and groups the primers into 2 groups: those that only 
    yield one amplicon and those that yield more than one amplicons.
    
    Parameters:
        LCSfile = file of least common substrings' (primers) data
        amplicon_size = size of amplified product in terms of number of
            fragments intervening 2 primers. For example, if the fragment 
            size is 15 and the amplicon size is 60, we are then looking
            for primers (forward primer sequence = reverse primer
            sequence) that gives amplicons of 900 bp (60*15) excluding
            the flanking primers.
    
    Returns:
        (List of unique primers,
        List of non-unique primers)
        where 'List of unique primers' are primers amplifying only one 
        amplicon, and 'List of non-unique primers' is a tuple of 'number
        of amplicons' and primer sequence.
    """
    LCS = get_LCS_by_amplicon(LCSfile, amplicon_size)[0]
    unique_LCS = [x
                  for x in LCS
                      if LCS.count(x) == 1]
    non_unique_LCS = list(set([(LCS.count(x), x)
                               for x in LCS
                                   if LCS.count(x) > 1]))
    return (unique_LCS, non_unique_LCS)
    
def LCS_tabulate(LCSfile, max_interval=198, min_interval=6, 
                 fragment_size=15, p='yes'):
    """
    Provides a tabulation of the given LCSfile, grouped by amplicon size.
    
    Parameters:
        LCSfile = file of least common substrings' (primers) data
        max_interval = number of maximum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 198 as 3000bp amplicon is the maximum normal
            Taq polymerase can amplify. 198 intervals of 15 bases per 
            fragment + 2 flanking fragment of maximum of 15 bases (primer
            length) gives a total of 200 fragments of 15 bases = 3000 bases.
        min_interval = number of minimum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 6 as 120bp amplicon is the minimum length visible
            on 1.5% agarose gel that is also suitable for 3000 bp. 6 intervals 
            of 15 bases per fragment + 2 flanking fragment of maximum of 15 
            bases (primer length) gives a total of 8 fragments of 15 bases = 
            120 bases.
        fragment_size = length of each sequence to be sliced into. This will
            be the maximum primer/probe length.
            Default value is 15 bases.
        p = flag to determine if the tabulation data is to be printed. 
            Default = yes
    """
    result = {}
    for amplicon_size in range(max_interval, min_interval, -1):
        (unique_LCS, non_unique_LCS) = LCS_uniqueness(LCSfile, amplicon_size)
        if p == 'yes':
            print 'Longest Common Substring (primer) for amplicon size of ' + \
                str(amplicon_size*fragment_size)
            print 'Unique Primers (only one amplicon). N = ' + \
                str(len(unique_LCS))
            print '  '.join(unique_LCS)
            print
            print 'Non-Unique Primers (more than one amplicon)'
            print '\t'.join(['# amplicons', 'Primer sequence'])
            for lcs in non_unique_LCS:
                print '\t\t'.join([str(lcs[0]), lcs[1]])
            print
            print
        result[str(amplicon_size)] = (unique_LCS, non_unique_LCS)
    return result

def inverse_LCS(LCSfile, inverseLCSfile, max_interval=198, min_interval=6, 
                fragment_size=15):
    """
    Generates the index file of LCSfile. 
    
    Format of LCSfile: <start fragment position>|
                       <end fragment position>|
                       <amplicon size> = <primer sequence>
    Format of inverseLCSfile: <primer sequence> =
                              [<start fragment position>|
                               <end fragment position>|
                               <amplicon size>, ...]
    
    Parameters:
        LCSfile = file of least common substrings' (primers) data
        inverseLCSfile = name of file to store the index of LCSfile
        max_interval = number of maximum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 198 as 3000bp amplicon is the maximum normal
            Taq polymerase can amplify. 198 intervals of 15 bases per 
            fragment + 2 flanking fragment of maximum of 15 bases (primer
            length) gives a total of 200 fragments of 15 bases = 3000 bases.
        min_interval = number of minimum fragment intervals between 2 primers. 
            This value is determines the maximum length of amplicon generated 
            and is determined by the length of each sliced fragment.
            Default value is 6 as 120bp amplicon is the minimum length visible
            on 1.5% agarose gel that is also suitable for 3000 bp. 6 intervals 
            of 15 bases per fragment + 2 flanking fragment of maximum of 15 
            bases (primer length) gives a total of 8 fragments of 15 bases = 
            120 bases.
        fragment_size = length of each sequence to be sliced into. This will
            be the maximum primer/probe length.
            Default value is 15 bases.
    """
    inverseLCSfile = marshaldbm(inverseLCSfile, 'c')
    count = 0
    for amplicon_size in range(max_interval, min_interval, -1):
        seq, pos = get_LCS_by_amplicon(LCSfile, amplicon_size)
        print 'Processing for amplicon size of ' + \
            str(amplicon_size*fragment_size)
        for index in range(len(seq)):
            count = count + 1
            try:
                temp = inverseLCSfile[seq[index]]
                inverseLCSfile[seq[index]] = temp + [pos[index]]
            except KeyError:
                inverseLCSfile[seq[index]] = [pos[index]]
            if (count % 1000) == 0:
                print str(count) + ' LCS processed'
    print str(count) + ' LCS processed'
    inverseLCSfile.close()

def look_for_primers(inverseLCSfile, num_of_amplicons, min_tm, p='yes'):
    """
    Scans the LCSfile index for primer(s) that amplifies a given number of
    amplicons. The primer must not end with adenosine and thymidine, and
    must meet minimum annealing temperature.
    
    Parameters:
        inverseLCSfile = name of LCSfile index file
        num_of_amplicons = number of amplicons generated by required primer(s)
        min_tm = minimum annealing temperature (in degrees centigrade) of 
            primers. Calculated as 2AT + 4GC.
        p = flag to determine if data is to be printed. Default = yes
    
    Returns:
        List of primers meeting the criteria
    """
    f = marshaldbm(inverseLCSfile, 'c')
    keys = f.keys()
    primers = []
    for k in keys:
        temperature = 4*(k.count('G')+k.count('C')) + \
                      2*(k.count('A')+k.count('T'))
        if (len(f[k])) == num_of_amplicons and \
           temperature > min_tm-1 and \
           not k.endswith('A') and not k.endswith('T'):
            primers.append(k)
            if p == 'yes':
                print k
                print f[k]
                print
    return primers
            
def generate_primers_from_fasta(fasta='ATCC8739.fasta',
                                genome_fragment='genome_fragment.table',
                                max_primer_length=15,
                                min_primer_length=7,
                                max_amplicon_length=3100,
                                min_amplicon_length=300,
                                primer_file='primer.table',
                                inverse_primer_file='invprimer.table'):
    """
    Generate a file of primers with amplicon size and genomic position from
    the Fasta sequence file.
    
    Parameters:
        fasta = name of fasta file to process
        genome_fragment = name of file to store the sliced sequence.
            Default = 'genome_fragment.table'
        max_primer_length = maximum length of primer (fragment size of genome
            slices). Default = 15 bp
        min_primer_length = minimum length of primer. Default = 7 bp
        max_amplicon_length = maximum size of amplicon. Default = 3100 bp
        min_amplicon_length = minimum size of amplicon. Default = 300 bp
        primer_file = file of least common substrings' (primers) data.
            Default = 'primer.table'
        inverse_primer_file = name of primers index file. 
            Default = 'invprimer.table'
    
    Output files:
        1. genome fragment file from 'generate_fragments' function
        2. primer file from 'generate_all_LCS' function
        3. inverse primer file from 'inverse_LCS' function
    """
    seq = fasta_seq(fasta)
    generate_fragments(seq, genome_fragment, fragment_size=max_primer_length)
    generate_all_LCS(genome_fragment, primer_file, 
                     min_primer_length, 
                     int(max_amplicon_length/max_primer_length), 
                     int(min_amplicon_length/max_primer_length))
    inverse_LCS(primer_file, inverse_primer_file, 
                int(max_amplicon_length/max_primer_length), 
                int(min_amplicon_length/max_primer_length), 
                max_primer_length)