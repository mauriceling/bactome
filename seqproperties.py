'''!
Nucleotide and Amino Acid Sequence Analysis Tool

Date created: 3rd May 2018

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
from Bio import SeqIO 

import fire


class CodonUsageBias(object):
    '''!
    Class to hold the core methods for codon usage bias analysis.
    '''

    def __init__(self):
        '''!
        Constructor method.
        '''
        self.seqNN = {}
        self.codonCount = {}
        self.aalist = ['A', 'C', 'D', 'E', 'F', 
                       'G', 'H', 'I', 'K', 'L', 
                       'M', 'N', 'P', 'Q', 'R', 
                       'S', 'T', 'V', 'W', 'Y', 
                       '*']

    def addSequencesFromFasta(self, fastafile):
        '''!
        Method to add sequences from FASTA file into data structure. 
        This allows for sequences to be added from multiple FASTA 
        files.

        @param fastafile String: Path to the FASTA file to add.
        '''
        for r in SeqIO.parse(fastafile, 'fasta'):
            self.seqNN[r.id] = [r.seq, r.description]

    def generateCodonCount(self, seq, genetic_code=1):
        '''!
        Method to generate codon counts from sequence. This method 
        will take the nucleotide sequence to translate into amino 
        acid sequence before codon count generation for each amino 
        acid, resulting in codon counts per amino acid.

        @param seq String: Nucleotide sequence string
        @param genetic_code Integer: Genetic code number to be used 
        for translation. Default = 1 (Standard Code). For more 
        information, see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>

        @return A nested dictionary of {'amino acid': {'codon': 
        count}}
        '''
        nnseq = str(seq)
        nnseq = [nnseq[i:3+i] for i in range(0, len(nnseq), 3)]
        aaseq = str(seq.translate(genetic_code))
        aaseq = [aaseq[i:1+i] for i in range(len(aaseq))]
        table = {'G': {}, 'P': {}, 'A': {}, 'V': {}, 'L': {},
                 'I': {}, 'M': {}, 'C': {}, 'F': {}, 'Y': {},
                 'W': {}, 'H': {}, 'K': {}, 'R': {}, 'Q': {}, 
                 'N': {}, 'E': {}, 'D': {}, 'S': {}, 'T': {},
                 '*': {}}
        for index in range(len(nnseq)):
            codon = nnseq[index]
            aa = aaseq[index]
            if codon in table[aa]:
                table[aa][codon] = table[aa][codon] + 1
            else:
                table[aa][codon] = 1
        return table

    def generateCodonCounts(self, genetic_code=1):
        '''!
        Method to iterate through all sequences to generate codon 
        counts.

        @param genetic_code Integer: Genetic code number to be used 
        for translation. Default = 1 (Standard Code). For more 
        information, see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
        '''
        for k in self.seqNN:
            codonCount = self.generateCodonCount(self.seqNN[k][0],
                                                 genetic_code)
            self.codonCount[k] = codonCount

def sequenceIDs(fastafile):
    '''!
    Function to print out the sequence IDs of all the FASTA records 
    in the FASTA file.

    Usage:

        python seqprop.py showIDs --fastafile=<FASTA file path>

    The output will be in the format of

        <count> : <sequence ID>

    where 
        - count is the numeric running order
        - sequence ID is the sequence ID of the FASTA record.

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    count = 1
    for k in o.seqNN:
        print('%i : %s' % (count, k))
        count = count + 1

def sequenceDescriptions(fastafile):
    '''!
    Function to print out the sequence IDs and descriptions of all the 
    FASTA records in the FASTA file.

    Usage:

        python seqprop.py showDesc --fastafile=<FASTA file path>

    The output will be in the format of

        <count> : <sequence ID> : <description>

    where 
        - count is the numeric running order
        - sequence ID is the sequence ID of the FASTA record
        - description is the description of the FASTA record

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    count = 1
    for k in o.seqNN:
        print('%i : %s : %s' % (count, k, o.seqNN[k][1]))
        count = count + 1

def translate(fastafile, genetic_code=1):
    '''!
    Function to translate all the FASTA records from nucleotide 
    sequence(s) to amino acid sequence(s).

    Usage:

        python seqprop.py translate --fastafile=<FASTA file path> 
        --genetic_code=<genetic code number>

    The output will be in the format of

        <sequence ID> : <amino acid sequence>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - amino acid sequence is the translated amino acid sequence 
        of the FASTA record

    @param fastafile String: Path to the FASTA file to be processed.
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        aaseq = o.seqNN[k][0].translate(genetic_code)
        print('%s : %s' % (k, str(aaseq)))

def aminoacidCount(fastafile, genetic_code=1):
    '''!
    Function to translate each nucleotide sequence (by FASTA record) 
    and generate a frequency table of the amino acids.

    Usage:

        python seqprop.py aacount --fastafile=<FASTA file path> 
        --genetic_code=<genetic code number>

    The output will be in the format of

        <sequence ID> : <A count> : <C count> : <D count> : <E count> 
        : <F count> : <G count> : <H count> : <I count> : <K count> : 
        <L count> : <M count> : <N count> : <P count> : <Q count> : 
        <R count> : <S count> : <T count> : <V count> : <W count> : 
        <Y count>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - the counts are the number of the respective amino acid; for 
        example, A count is the number of alanine in the peptide

    @param fastafile String: Path to the FASTA file to be processed.
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    header = ' : '.join(['SequenceID', 'A', 'C', 'D', 'E', 'F', 'G', 
        'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 
        'W', 'Y'])
    print(header)
    for k in o.seqNN:
        aaseq = o.seqNN[k][0].translate(genetic_code)
        aaseq = [aaseq[i:1+i] for i in range(len(aaseq))]
        data = [k]
        for aa in ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 
                   'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']:
            data.append(len([str(x) for x in aaseq if str(x) == aa]))
        data = ' : '.join([str(x) for x in data])
        print(data)

def peptideLength(fastafile, genetic_code=1):
    '''!
    Function to translate each nucleotide sequence (by FASTA record) 
    and count the number of amino acids (peptide length).

    Usage:

        python seqprop.py plength --fastafile=<FASTA file path> 
        --genetic_code=<genetic code number>

     The output will be in the format of

        <sequence ID> : <peptide length>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - peptide length is the number of amino acids in the peptide

    @param fastafile String: Path to the FASTA file to be processed.
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        aaseq = o.seqNN[k][0].translate(genetic_code)
        print('%s : %s' % (k, str(len(str(aaseq)))))

def nucleotideLength(fastafile):
    '''!
    Function to count the number of nucleotides (number of bases) in 
    each FASTA record.

    Usage:

        python seqprop.py nlength --fastafile=<FASTA file path> 

     The output will be in the format of

        <sequence ID> : <nucleotide length>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - nucleotide length is the number of bases in the peptide

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        nnseq = o.seqNN[k][0]
        print('%s : %s' % (k, str(len(str(nnseq)))))

def flattenCodonCount(CC):
    '''!
    Function to flatten the codon usage / frequency table by 
    removing the amino acid identifiers.

    @param CC: Codon usage table as a nested dictionary 
    of {'amino acid': {'codon': count}}.

    @return Dictionary of {'codon': count}.
    '''
    table = {'AAA': 0, 'AAT': 0, 'AAG': 0, 'AAC': 0,
             'ATA': 0, 'ATT': 0, 'ATG': 0, 'ATC': 0,
             'AGA': 0, 'AGT': 0, 'AGG': 0, 'AGC': 0,
             'ACA': 0, 'ACT': 0, 'ACG': 0, 'ACC': 0,
             'TAA': 0, 'TAT': 0, 'TAG': 0, 'TAC': 0,
             'TTA': 0, 'TTT': 0, 'TTG': 0, 'TTC': 0,
             'TGA': 0, 'TGT': 0, 'TGG': 0, 'TGC': 0,
             'TCA': 0, 'TCT': 0, 'TCG': 0, 'TCC': 0,
             'GAA': 0, 'GAT': 0, 'GAG': 0, 'GAC': 0,
             'GTA': 0, 'GTT': 0, 'GTG': 0, 'GTC': 0,
             'GGA': 0, 'GGT': 0, 'GGG': 0, 'GGC': 0,
             'GCA': 0, 'GCT': 0, 'GCG': 0, 'GCC': 0,
             'CAA': 0, 'CAT': 0, 'CAG': 0, 'CAC': 0,
             'CTA': 0, 'CTT': 0, 'CTG': 0, 'CTC': 0,
             'CGA': 0, 'CGT': 0, 'CGG': 0, 'CGC': 0,
             'CCA': 0, 'CCT': 0, 'CCG': 0, 'CCC': 0}
    for aa in CC:
        for codon in CC[aa]:
            table[codon] = CC[aa][codon]
    return table

def codonCount(fastafile, genetic_code=1):
    '''!
    Function to generate the codon usage frequency table by each 
    FASTA record.

    Usage:

        python seqprop.py codoncount --fastafile=<FASTA file path> 
        --genetic_code=<genetic code number>

    The output will be in the format of

        <sequence ID> : <AAA count> : <AAC count> : <AAG count> : 
        <AAT count> : <ACA count> : <ACC count> : <ACG count> : 
        <ACT count> : <AGA count> : <AGC count> : <AGG count> : 
        <AGT count> : <ATA count> : <ATC count> : <ATG count> : 
        <ATT count> : <CAA count> : <CAC count> : <CAG count> : 
        <CAT count> : <CCA count> : <CCC count> : <CCG count> : 
        <CCT count> : <CGA count> : <CGC count> : <CGG count> : 
        <CGT count> : <CTA count> : <CTC count> : <CTG count> : 
        <CTT count> : <GAA count> : <GAC count> : <GAG count> : 
        <GAT count> : <GCA count> : <GCC count> : <GCG count> : 
        <GCT count> : <GGA count> : <GGC count> : <GGG count> : 
        <GGT count> : <GTA count> : <GTC count> : <GTG count> : 
        <GTT count> : <TAA count> : <TAC count> : <TAG count> : 
        <TAT count> : <TCA count> : <TCC count> : <TCG count> : 
        <TCT count> : <TGA count> : <TGC count> : <TGG count> : 
        <TGT count> : <TTA count> : <TTC count> : <TTG count> : 
        <TTT count>

     where 
        - sequence ID is the sequence ID of the FASTA record
        - the counts are the number of each codon; for example, 
        AAA count is the number of AAA codon in the sequence

    @param fastafile String: Path to the FASTA file to be processed.
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    header = ' : '.join(['SequenceID', 'AAA', 'AAC', 'AAG', 'AAT', 
        'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA',
        'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 
        'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 
        'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 
        'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAA', 
        'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGA', 'TGC', 
        'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT'])
    print(header)
    for k in o.seqNN:
        try:
            codonCount = o.generateCodonCount(o.seqNN[k][0],
                                              genetic_code)
            fCCount = flattenCodonCount(codonCount)
            data = [k, fCCount['AAA'], fCCount['AAC'], fCCount['AAG'], 
                    fCCount['AAT'], fCCount['ACA'], fCCount['ACC'], 
                    fCCount['ACG'], fCCount['ACT'], fCCount['AGA'], 
                    fCCount['AGC'], fCCount['AGG'], fCCount['AGT'], 
                    fCCount['ATA'], fCCount['ATC'], fCCount['ATG'], 
                    fCCount['ATT'], fCCount['CAA'], fCCount['CAC'], 
                    fCCount['CAG'], fCCount['CAT'], fCCount['CCA'], 
                    fCCount['CCC'], fCCount['CCG'], fCCount['CCT'], 
                    fCCount['CGA'], fCCount['CGC'], fCCount['CGG'], 
                    fCCount['CGT'], fCCount['CTA'], fCCount['CTC'], 
                    fCCount['CTG'], fCCount['CTT'], fCCount['GAA'], 
                    fCCount['GAC'], fCCount['GAG'], fCCount['GAT'], 
                    fCCount['GCA'], fCCount['GCC'], fCCount['GCG'], 
                    fCCount['GCT'], fCCount['GGA'], fCCount['GGC'], 
                    fCCount['GGG'], fCCount['GGT'], fCCount['GTA'], 
                    fCCount['GTC'], fCCount['GTG'], fCCount['GTT'], 
                    fCCount['TAA'], fCCount['TAC'], fCCount['TAG'], 
                    fCCount['TAT'], fCCount['TCA'], fCCount['TCC'], 
                    fCCount['TCG'], fCCount['TCT'], fCCount['TGA'], 
                    fCCount['TGC'], fCCount['TGG'], fCCount['TGT'], 
                    fCCount['TTA'], fCCount['TTC'], fCCount['TTG'], 
                    fCCount['TTT']]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except: pass

def percentGC(fastafile):
    '''!
    Function to generate the %GC by each FASTA record.

    Usage:

        python seqprop.py gc --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%GC>

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        percent = sequence.count('G') + sequence.count('C') + \
                  sequence.count('g') + sequence.count('c')
        percent = percent / len(sequence)
        data = [k, percent]
        data = ' : '.join([str(x) for x in data])
        print(data)

def percentG(fastafile):
    '''!
    Function to generate the %G by each FASTA record.

    Usage:

        python seqprop.py g --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%G>

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        percent = sequence.count('G') + sequence.count('g')
        percent = percent / len(sequence)
        data = [k, percent]
        data = ' : '.join([str(x) for x in data])
        print(data)

def percentA(fastafile):
    '''!
    Function to generate the %A by each FASTA record.

    Usage:

        python seqprop.py a --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%A>

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        percent = sequence.count('A') + sequence.count('a')
        percent = percent / len(sequence)
        data = [k, percent]
        data = ' : '.join([str(x) for x in data])
        print(data)

def percentGCi(fastafile, i, j=3):
    '''!
    Function to generate the %GC of the i-th base in each codon 
    (of j length) by each FASTA record.

    Usage:

        python seqprop.py gci --i=1 --j=3 --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%GC of i-th base in codon of j bases>

    @param fastafile String: Path to the FASTA file to be processed.
    @param i Integer: Position of the base in the codon.
    @param j Integer: Size (length) of each codon. Default = 3.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    i = int(i)
    j = int(j)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        sequence = [sequence[i:i+j] 
                    for i in range(0, len(sequence), j)]
        temp = []
        for x in sequence:
            try: temp.append(x[i-1])
            except IndexError: pass
        sequence = temp
        percent = sequence.count('G') + sequence.count('C') + \
                  sequence.count('g') + sequence.count('c')
        percent = percent / len(sequence)
        data = ' : '.join([str(k), str(percent)])
        print(data)

def percentGi(fastafile, i, j=3):
    '''!
    Function to generate the %G of the i-th base in each codon 
    (of j length) by each FASTA record.

    Usage:

        python seqprop.py gi --i=1 --j=3 --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%G of i-th base in codon of j bases>

    @param fastafile String: Path to the FASTA file to be processed.
    @param i Integer: Position of the base in the codon.
    @param j Integer: Size (length) of each codon. Default = 3.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    i = int(i)
    j = int(j)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        sequence = [sequence[i:i+j] 
                    for i in range(0, len(sequence), j)]
        temp = []
        for x in sequence:
            try: temp.append(x[i-1])
            except IndexError: pass
        sequence = temp
        percent = sequence.count('G') + sequence.count('g')
        percent = percent / len(sequence)
        data = ' : '.join([str(k), str(percent)])
        print(data)

def percentAi(fastafile, i, j=3):
    '''!
    Function to generate the %A of the i-th base in each codon 
    (of j length) by each FASTA record.

    Usage:

        python seqprop.py ai --i=1 --j=3 --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <%A of i-th base in codon of j bases>

    @param fastafile String: Path to the FASTA file to be processed.
    @param i Integer: Position of the base in the codon.
    @param j Integer: Size (length) of each codon. Default = 3.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    i = int(i)
    j = int(j)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        sequence = [sequence[i:i+j] 
                    for i in range(0, len(sequence), j)]
        temp = []
        for x in sequence:
            try: temp.append(x[i-1])
            except IndexError: pass
        sequence = temp
        percent = sequence.count('A') + sequence.count('a')
        percent = percent / len(sequence)
        data = ' : '.join([str(k), str(percent)])
        print(data)


if __name__ == '__main__':
    exposed_functions = {'showIDs': sequenceIDs,
                         'showDesc': sequenceDescriptions,
                         'plength': peptideLength,
                         'nlength': nucleotideLength,
                         'translate': translate,
                         'codoncount': codonCount,
                         'aacount': aminoacidCount,
                         'gc': percentGC,
                         'g': percentG,
                         'a': percentA,
                         'gci': percentGCi,
                         'gi': percentGi,
                         'ai': percentAi}
    fire.Fire(exposed_functions)