'''!
Random Amino Acid and Nucleotide Sequence Generator

Date created: 8th May 2018

License: GNU General Public License version 3 for academic or 
not-for-profit use only


Bactome package is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as 
published by the Free Software Foundation, either version 3 of the 
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
'''
import random

import fire

class RandomSequence(object):
    '''!
    Class to generate random sequences based on input parameters. 
    Although the initial intention is to generate DNA/RNA sequences, 
    it can be used to generate other random sequences; such as, 
    amino acid sequences or even non-biological sequences.

    Basic example:

        >>> s = RandomSequence()
        >>> selection = {'A':25, 'U':25, 'G':25, 'C':25}
        >>> s.initiateNucleotide(selection)
        >>> s.generateSequence(50, True, True)
        'AGAUCCCCCGGACAGGUUAUGGAUAAGUGUCCACCUCUAACUAUUUCUGC'

    Example without start and stop:

        >>> s = RandomSequence()
        >>> selection = {'A':25, 'T':25, 'G':25, 'C':25}
        >>> s.initiateNucleotide(selection)
        >>> s.initiateStart('TTG,CTG,ATG')
        >>> s.initiateStop('TAA,TAG,TGA')
        >>> s.generateSequence(50, False, False)
        'AATCCGGACAAATCAAGCACAGGCAGTTTATTCAAGGGGTACCCCAGAAC'
    '''

    def __init__(self):
        '''!
        Constructor method.
        '''
        self.start_codons = []
        self.stop_codons = []
        self.sseq = []

    def initiateNucleotide(self, ndict):
        '''!
        Initialize the object with a bag of items (atomic sequence), 
        defined as a dictionary, for random sequence generation. The 
        atomic sequences will be randomized.

        @param ndict Dictionary: Represents a bag of items where the 
        key is the atomic sequence and the value is the relative 
        occurrence. For example, {'A':20, 'T':20, 'G':30, 'C':30} 
        will initialize for the generation of a random sequence with 
        60% GC.
        '''
        for k in ndict:
            self.sseq = self.sseq + [k] * int(ndict[k])
            random.shuffle(self.sseq)
        for i in range(10):
            random.shuffle(self.sseq)

    def initiateStart(self, start_codons='TTG,CTG,ATG'):
        '''!
        Initialize a list of start codons. This is used for two 
        purposes. Firstly, it can be used to check the randomly 
        generated sequence to ensure that there is no start codons 
        within the sequence. Secondly, it can be used to cap the 
        start of a newly generated sequence.

        @param start_codons String: A comma-delimited list to 
        represent start codons. Default = 'TTG,CTG,ATG'.
        '''
        start_codons = start_codons.strip()
        start_codons = [x.strip() for x in start_codons.split(',')]
        self.start_codons = start_codons

    def initiateStop(self, stop_codons='TAA,TAG,TGA'):
        '''!
        Initialize a list of stop codons. This is used for two 
        purposes. Firstly, it can be used to check the randomly 
        generated sequence to ensure that there is no stop codons 
        within the sequence. Secondly, it can be used to cap the 
        end of a newly generated sequence.

        @param stop_codons String: A comma-delimited list to 
        represent start codons. Default = 'TAA,TAG,TGA'.
        '''
        stop_codons = stop_codons.strip()
        stop_codons = [x.strip() for x in stop_codons.split(',')]
        self.stop_codons = stop_codons

    def _cleanStop(self, seq):
        '''!
        Private method - to remove all occurrences of stop codons 
        in a sequence.

        @param seq String: Sequence to clean of all stop codons.

        @return Cleaned sequence.
        '''
        for stop in self.stop_codons:
            seq = seq.replace(stop, '')
        return seq

    def _cleanStart(self, seq):
        '''!
        Private method - to remove all occurrences of start codons 
        in a sequence.

        @param seq String: Sequence to clean of all start codons.

        @return Cleaned sequence.
        '''
        for start in self.start_codons:
            seq = seq.replace(start, '')
        return seq

    def generateSequence(self, length, allow_start=False,
                         allow_stop=False):
        '''!
        Generate a random sequence.

        @param length Integer: Length of random sequence to generate. 
        @param allow_start Boolean: Flag to allow occurrence(s) of 
        start codon(s) within the randomly generated sequence. Default 
        = False.
        @param allow_stop Boolean: Flag to allow occurrence(s) of 
        stop codon(s) within the randomly generated sequence. Default 
        = False.

        @return Randomly generated sequence of specified length.
        '''
        length = int(length)
        sequence = ''
        while len(sequence) < length:
            sequence = sequence + random.choice(self.sseq)
            if not allow_start:
                sequence = self._cleanStart(sequence)
            if not allow_stop:
                sequence = self._cleanStop(sequence)
        return sequence[:length]

def selectionGenerator(selection):
    '''!
    Function to convert a string-based atomic sequence(s) definition 
    into a dictionary, suitable as parameter for 
    RandomSequence.initiateNucleotide() function.

    @param selection String: Definition of the atomic sequence(s), 
    defined as a semi-colon delimited definition of comma-delimited 
    definition, used for random sequence generation.
    '''
    selection = str(selection)
    selection = [pair.strip() for pair in selection.split(';')]
    selection = [[pair.split(',')[0], pair.split(',')[1]] 
                 for pair in selection]
    selection = [[pair[0].strip(), pair[1].strip()] 
                 for pair in selection]
    ndict = {}
    for k in selection:
        ndict[str(k[0])] = int(k[1])
    return ndict

def _initial_RandomSequence(start_codons, stop_codons,
                            selection, source_seq):
    '''!
    Private method - Instantiate a RandomSequence object based on 
    parameters. Main requirement is either selection or source_seq. 
    If both selection and source_seq are given (not empty) strings, 
    only selection will be used.

    @param start_codons String: A comma-delimited list to represent 
    start codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no start codons within the sequence. Secondly, it can be used 
    to cap the start of a newly generated sequence. 
    @param stop_codons String: A comma-delimited list to represent 
    stop codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no stop codons within the sequence. Secondly, it can be used 
    to cap the stop of a newly generated sequence.Default = 
    'TAA,TAG,TGA'.
    @param selection String: Definition of the atomic sequence(s), 
    defined as a semi-colon delimited definition of comma-delimited 
    definition, used for random sequence generation.
    @param source_seq String: A sequence to be used to generate the 
    atomic sequence(s), used for random sequence generation. The 
    atomic sequence(s) is/are essentially an single-character 
    alphanumeric of the sequence, where relative abundance of each 
    character is determined by the occurrence in the sequence.

    @return RandomSequence object.
    '''
    o = RandomSequence()
    o.initiateStart(start_codons)
    o.initiateStop(stop_codons)
    if source_seq == '':
        ndict = selectionGenerator(selection)
    else:
        source_seq = [x for x in source_seq]
        ndict = {}
        for k in list(set(source_seq)):
            ndict[k] = source_seq.count(k)
    o.initiateNucleotide(ndict)
    return o

def _generate_sequence(o, min_length, max_length, 
                       allow_start, allow_stop,
                       cap_start, cap_stop):
    '''!
    Private method - Generates a fixed or variable length random 
    sequence based on given RandomSequence object. For fixed length 
    random sequence generation, min_length == max_length. If 
    min_length != max_length, variable length random sequence will 
    be generated with length uniformly distributed between min_length 
    and max_length.

    @param o Object: RandomSequence object.
    @param min_length Integer: Minimum length of randomly generated 
    sequence.
    @param max_length Integer: Maximum length of randomly generated 
    sequence.
    @param allow_start Boolean: Flag to allow occurrence(s) of 
    start codon(s) within the randomly generated sequence. 
    @param allow_stop Boolean: Flag to allow occurrence(s) of 
    stop codon(s) within the randomly generated sequence. 
    @param cap_start Boolean: Flag to determine whether to cap the 
    start of the randomly generated sequence with a start codon. 
    @param cap_stop Boolean: Flag to determine whether to cap the 
    end of the randomly generated sequence with a stop codon. 

    @return Generated random sequence.
    '''
    min_length = int(min_length)
    max_length = int(max_length)
    if min_length == max_length:
        sequence = o.generateSequence(max_length, allow_start, 
                                      allow_stop)
    else:
        sequence = o.generateSequence(max_length, allow_start, 
                                      allow_stop)
        length = min_length + \
                 int((max_length - min_length) * random.random())
        sequence = sequence[:length+1]
    if cap_start:
        sequence = random.choice(o.start_codons) + sequence
    if cap_stop:
        sequence = sequence + random.choice(o.stop_codons)
    return sequence

def gFixedLength(length, n, allow_start=False, allow_stop=False, 
                 start_codons='TTG,CTG,ATG', cap_start=False,
                 stop_codons='TAA,TAG,TGA', cap_stop=False,
                 selection='A,250;T,250;G,250;C,250', 
                 source_seq='', fasta=True, prefix='Test'):
    '''!
    Function to generate one or more fixed length random sequences as 
    per specification.
    
    Usage:

        python randomseq.py FLS --length=100 --n=10 --allow_start=False --allow_stop=False --start_codons='TTG,CTG,ATG' --stop_codons='TAA,TAG,TGA' --cap_start=True --cap_stop=True --selection=A,250;T,250;G,250;C,250 --fasta=True, --prefix='Test'

    @param length Integer: Length of random sequence to generate. 
    @param n Integer: Number of random sequence(s) to generate.
    @param allow_start Boolean: Flag to allow occurrence(s) of start 
    codon(s) within the randomly generated sequence. Default = False.
    @param allow_stop Boolean: Flag to allow occurrence(s) of stop 
    codon(s) within the randomly generated sequence. Default = False.
    @param start_codons String: A comma-delimited list to represent 
    start codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no start codons within the sequence. Secondly, it can be used 
    to cap the start of a newly generated sequence.Default = 
    'TTG,CTG,ATG'.
    @param cap_start Boolean: Flag to determine whether to cap the 
    start of the randomly generated sequence with a start codon. 
    Default = False.
    @param stop_codons String: A comma-delimited list to represent 
    stop codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no stop codons within the sequence. Secondly, it can be used 
    to cap the stop of a newly generated sequence.Default = 
    'TAA,TAG,TGA'.
    @param cap_stop Boolean: Flag to determine whether to cap the 
    end of the randomly generated sequence with a stop codon. 
    Default = False.
    @param selection String: Definition of the atomic sequence(s), 
    defined as a semi-colon delimited definition of comma-delimited 
    definition, used for random sequence generation. Default = 
    'A,250;T,250;G,250;C,250' which means that the list of atomic 
    sequences is defined as 250 'A', 250 'T', 250 'G', and 250 'C'.
    @param source_seq String: A sequence to be used to generate the 
    atomic sequence(s), used for random sequence generation. The 
    atomic sequence(s) is/are essentially an single-character 
    alphanumeric of the sequence, where relative abundance of each 
    character is determined by the occurrence in the sequence.
    @param fasta Boolean: Flag to determine whether the output is to 
    be formatted as a FASTA file. Default = True.
    @param prefix String: Prefix to the running number for the 
    sequence name, if the output is to be a FASTA format. Default = 
    'Test'.
    '''
    o = _initial_RandomSequence(start_codons, stop_codons,
                                selection, source_seq)
    for i in range(int(n)):
        sequence = _generate_sequence(o, length, length, allow_start, 
                                      allow_stop, cap_start, cap_stop)
        if fasta:
            title = '_'.join([str(prefix), str(i+1)])
            print('> %s' % title)
        print(sequence)

def gVariableLength(min_length, max_length, n, 
                    allow_start=False, allow_stop=False, 
                    start_codons='TTG,CTG,ATG', cap_start=False,
                    stop_codons='TAA,TAG,TGA', cap_stop=False,
                    selection='A,250;T,250;G,250;C,250', 
                    source_seq='', fasta=True, prefix='Test'):
    '''!
    Function to generate one or more variable length random sequences 
    as per specification.
    
    Usage:

        python randomseq.py VLS --min_length=90 --max_length=110 --n=10 --allow_start=False --allow_stop=False --start_codons='TTG,CTG,ATG' --stop_codons='TAA,TAG,TGA' --cap_start=True --cap_stop=True --selection=A,250;T,250;G,250;C,250 --fasta=True, --prefix='Test'

    @param min_length Integer: Minimum length of randomly generated 
    sequence.
    @param max_length Integer: Maximum length of randomly generated 
    sequence.
    @param n Integer: Number of random sequence(s) to generate.
    @param allow_start Boolean: Flag to allow occurrence(s) of start 
    codon(s) within the randomly generated sequence. Default = False.
    @param allow_stop Boolean: Flag to allow occurrence(s) of stop 
    codon(s) within the randomly generated sequence. Default = False.
    @param start_codons String: A comma-delimited list to represent 
    start codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no start codons within the sequence. Secondly, it can be used 
    to cap the start of a newly generated sequence.Default = 
    'TTG,CTG,ATG'.
    @param cap_start Boolean: Flag to determine whether to cap the 
    start of the randomly generated sequence with a start codon. 
    Default = False.
    @param stop_codons String: A comma-delimited list to represent 
    stop codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no stop codons within the sequence. Secondly, it can be used 
    to cap the stop of a newly generated sequence.Default = 
    'TAA,TAG,TGA'.
    @param cap_stop Boolean: Flag to determine whether to cap the 
    end of the randomly generated sequence with a stop codon. 
    Default = False.
    @param selection String: Definition of the atomic sequence(s), 
    defined as a semi-colon delimited definition of comma-delimited 
    definition, used for random sequence generation. Default = 
    'A,250;T,250;G,250;C,250' which means that the list of atomic 
    sequences is defined as 250 'A', 250 'T', 250 'G', and 250 'C'.
    @param source_seq String: A sequence to be used to generate the 
    atomic sequence(s), used for random sequence generation. The 
    atomic sequence(s) is/are essentially an single-character 
    alphanumeric of the sequence, where relative abundance of each 
    character is determined by the occurrence in the sequence.
    @param fasta Boolean: Flag to determine whether the output is to 
    be formatted as a FASTA file. Default = True.
    @param prefix String: Prefix to the running number for the 
    sequence name, if the output is to be a FASTA format. Default = 
    'Test'.
    '''
    o = _initial_RandomSequence(start_codons, stop_codons,
                                selection, source_seq)
    for i in range(int(n)):
        sequence = _generate_sequence(o, min_length, max_length, 
                                      allow_start, allow_stop, 
                                      cap_start, cap_stop)
        if fasta:
            title = '_'.join([str(prefix), str(i+1)])
            print('> %s' % title)
        print(sequence)
    
def shuffle(sequence):
     '''!
    Function to shuffle a given sequence.
    
    Usage:

        python randomseq.py shuffle --sequence=<sequence to shuffle>

    @param sequence String: Sequence to shuffle
    '''
    sequence = list(sequence)
    random.shuffle(sequence)
    return ''.join(sequence)


if __name__ == '__main__':
    exposed_functions = {'FLS': gFixedLength,
                         'VLS': gVariableLength,
                         'shuffle': shuffle}
    fire.Fire(exposed_functions)
