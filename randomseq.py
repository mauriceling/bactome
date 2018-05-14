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
        Initialize the object with a bag of items (atomic sequence) 
        from 

        @param ndict Dictionary: 
        '''
        for k in ndict:
            self.sseq = self.sseq + [k] * int(ndict[k])
            random.shuffle(self.sseq)
        for i in range(10):
            random.shuffle(self.sseq)

    def initiateStart(self, start_codons='TTG,CTG,ATG'):
        '''!
        '''
        start_codons = start_codons.strip()
        start_codons = [x.strip() for x in start_codons.split(',')]
        self.start_codons = start_codons

    def initiateStop(self, stop_codons='TAA,TAG,TGA'):
        '''!
        '''
        stop_codons = stop_codons.strip()
        stop_codons = [x.strip() for x in stop_codons.split(',')]
        self.stop_codons = stop_codons

    def _cleanStop(self, seq):
        '''!
        '''
        for stop in self.stop_codons:
            seq = seq.replace(stop, '')
        return seq

    def _cleanStart(self, seq):
        '''!
        '''
        for start in self.start_codons:
            seq = seq.replace(start, '')
        return seq

    def generateSequence(self, length, allow_start=False,
                         allow_stop=False):
        '''!
        '''
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

def gFixedLength(length, n, allow_stop=False, allow_start=False,
                 start_codons='TTG,CTG,ATG', cap_start=False,
                 stop_codons='TAA,TAG,TGA', cap_stop=False,
                 selection='A,250;T,250;G,250;C,250', 
                 source_seq='', fasta=True, prefix='Test'):
    '''!
    
    Usage:

        python randomseq.py FLS --length=100 --n=10 --allow_start=False --allow_stop=False --start_codons='TTG,CTG,ATG' --stop_codons='TAA,TAG,TGA' --cap_start=True --cap_stop=True --selection=A,250;T,250;G,250;C,250 --fasta=True, --prefix='Test'
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
                    allow_stop=False, allow_start=False,
                    start_codons='TTG,CTG,ATG', cap_start=False,
                    stop_codons='TAA,TAG,TGA', cap_stop=False,
                    selection='A,250;T,250;G,250;C,250', 
                    source_seq='', fasta=True, prefix='Test'):
    '''!
    
    Usage:

        python randomseq.py VLS --min_length=90 --max_length=110 --n=10 --allow_start=False --allow_stop=False --start_codons='TTG,CTG,ATG' --stop_codons='TAA,TAG,TGA' --cap_start=True --cap_stop=True --selection=A,250;T,250;G,250;C,250 --fasta=True, --prefix='Test'
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
    

if __name__ == '__main__':
    exposed_functions = {'FLS': gFixedLength,
                         'VLS': gVariableLength}
    fire.Fire(exposed_functions)
