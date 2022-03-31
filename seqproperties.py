'''!
Nucleotide and Amino Acid Sequence Analysis Tool

Date created: 3rd May 2018

License: GNU General Public License version 3 for academic or 
not-for-profit use only

Reference: Ling, MHT. 2020. SeqProperties: A Python Command-Line Tool 
for Basic Sequence Analysis. Acta Scientific Microbiology 3(6): 103-106.

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
import subprocess
import sys

try: 
    from Bio import Align
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'biopython',
                           '--trusted-host', 'pypi.org', 
                           '--trusted-host', 'files.pythonhosted.org'])
    from Bio import Align
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqUtils.ProtParam import ProteinAnalysis

import Bio
biopython_version = float(Bio.__version__)
if biopython_version < 1.78:
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna

try: 
    import fire
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'fire',
                           '--trusted-host', 'pypi.org', 
                           '--trusted-host', 'files.pythonhosted.org'])
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

        python seqproperties.py showIDs --fastafile=<FASTA file path>

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

        python seqproperties.py showDesc --fastafile=<FASTA file path>

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

        python seqproperties.py translate --fastafile=<FASTA file path> --genetic_code=<genetic code number>

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

def aminoacidCount(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to translate each nucleotide sequence (by FASTA record) 
    and generate a frequency table of the amino acids.

    Usage:

        python seqproperties.py aacount --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <A count> : <C count> : <D count> : <E count> : <F count> : <G count> : <H count> : <I count> : <K count> : <L count> : <M count> : <N count> : <P count> : <Q count> : <R count> : <S count> : <T count> : <V count> : <W count> : <Y count>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - the counts are the number of the respective amino acid; for 
        example, A count is the number of alanine in the peptide

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    header = ' : '.join(['SequenceID', 'A', 'C', 'D', 'E', 'F', 'G', 
        'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 
        'W', 'Y'])
    print(header)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        aacount = sequence.count_amino_acids()
        try:
            data = [k]
            for aa in ['A', 'C', 'D', 'E', 'F', 
                       'G', 'H', 'I', 'K', 'L', 
                       'M', 'N', 'P', 'Q', 'R', 
                       'S', 'T', 'V', 'W', 'Y']:
                data.append(aacount[aa])
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except IOError:
            data = ' : '.join([str(k), 'Error'])
            print(data)
        
def genericCount(fastafile):
    '''!
    Function to count the frequency of each character in each 
    nucleotide sequence (by FASTA record) and generate a frequency 
    table.

    Usage:

        python seqproperties.py count --fastafile=<FASTA file path>

    The output will be in the format of

        <sequence ID> : <length of sequence> : [list of counts delimited by " : "]

    where 
        - sequence ID is the sequence ID of the FASTA record
        - the counts are the number of the respective characters

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    char_set = set()
    for k in o.seqNN:
        sequence = [x for x in str(o.seqNN[k][0])]
        char_set.update(sequence)
    char_set = list(char_set)
    char_set.sort()
    header = ['SequenceID', 'Length'] + [c.upper() for c in char_set]
    header = ' : '.join(header)
    print(header)
    for k in o.seqNN:
        sequence = str(o.seqNN[k][0])
        data = [k, len(sequence)] + \
               [sequence.count(c) for c in char_set]
        data = ' : '.join([str(x) for x in data])
        print(data)

def peptideLength(fastafile):
    '''!
    Function to count the number of amino acids (peptide length) by 
    peptide FASTA record.

    Usage:

        python seqproperties.py plength --fastafile=<FASTA file path> 

     The output will be in the format of

        <sequence ID> : <peptide length>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - peptide length is the number of amino acids in the peptide

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        aaseq = o.seqNN[k][0]
        print('%s : %s' % (k, str(len(str(aaseq)))))

def nucleotideLength(fastafile):
    '''!
    Function to count the number of nucleotides (number of bases) in 
    each FASTA record.

    Usage:

        python seqproperties.py nlength --fastafile=<FASTA file path> 

     The output will be in the format of

        <sequence ID> : <nucleotide length>

    where 
        - sequence ID is the sequence ID of the FASTA record
        - nucleotide length is the number of bases in the sequence

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        nnseq = o.seqNN[k][0]
        print('%s : %s' % (k, str(len(str(nnseq)))))

def complement(fastafile):
    '''!
    Function to generate the complement sequence of each FASTA record. 
    This is done using reverse_complement function in Biopython - each 
    FASTA sequence is assumed to be in 5'-->3', this function will 
    generate the complementary sequence in 5'-->3' orientation rather 
    than 3'<--5' orientation.

    Usage:

        python seqproperties.py complement --fastafile=<FASTA file path> 

    The output will be in FASTA format.

    @param fastafile String: Path to the FASTA file to be processed.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        nnseq = o.seqNN[k][0]
        print("> %s" % k)
        print(str(nnseq.reverse_complement()))

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

        python seqproperties.py codoncount --fastafile=<FASTA file path> --genetic_code=<genetic code number>

    The output will be in the format of

        <sequence ID> : <AAA count> : <AAC count> : <AAG count> : <AAT count> : <ACA count> : <ACC count> : <ACG count> : <ACT count> : <AGA count> : <AGC count> : <AGG count> : <AGT count> : <ATA count> : <ATC count> : <ATG count> : <ATT count> : <CAA count> : <CAC count> : <CAG count> : <CAT count> : <CCA count> : <CCC count> : <CCG count> : <CCT count> : <CGA count> : <CGC count> : <CGG count> : <CGT count> : <CTA count> : <CTC count> : <CTG count> : <CTT count> : <GAA count> : <GAC count> : <GAG count> : <GAT count> : <GCA count> : <GCC count> : <GCG count> : <GCT count> : <GGA count> : <GGC count> : <GGG count> : <GGT count> : <GTA count> : <GTC count> : <GTG count> : <GTT count> : <TAA count> : <TAC count> : <TAG count> : <TAT count> : <TCA count> : <TCC count> : <TCG count> : <TCT count> : <TGA count> : <TGC count> : <TGG count> : <TGT count> : <TTA count> : <TTC count> : <TTG count> : <TTT count>

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

        python seqproperties.py gc --fastafile=<FASTA file path>

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

        python seqproperties.py g --fastafile=<FASTA file path>

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

        python seqproperties.py a --fastafile=<FASTA file path>

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

        python seqproperties.py gci --i=1 --j=3 --fastafile=<FASTA file path>

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

        python seqproperties.py gi --i=1 --j=3 --fastafile=<FASTA file path>

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

        python seqproperties.py ai --i=1 --j=3 --fastafile=<FASTA file path>

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

def _toPeptide(sequence, molecule, genetic_code=1, to_stop=True):
    '''
    Private function - Takes a sequence (DNA/RNA/amino acid) and 
    process it according to return a ProteinAnalysis object.

    @param sequence String: Nucleotide (DNA/RNA) or amino acid 
    sequence.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    @return: Bio.SeqUtils.ProtParam.ProteinAnalysis object
    '''
    if molecule.lower() == 'peptide':
        peptide = ProteinAnalysis(sequence)
    elif molecule.lower() == 'rna':
        rna = str(sequence)
        rna = Seq(rna, generic_rna)
        peptide = rna.translate(genetic_code, to_stop=to_stop)
        peptide = ProteinAnalysis(str(peptide))
    elif molecule.lower() == 'dna':
        dna = str(sequence)
        dna = Seq(dna, generic_dna)
        rna = dna.transcribe()
        peptide = rna.translate(genetic_code, to_stop=to_stop)
        peptide = ProteinAnalysis(str(peptide))
    return peptide

def molecularWeight(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the molecular weight, using Biopython, by 
    each FASTA record.

    Usage:

        python seqproperties.py mw --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <molecular weight>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        try:
            sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                                  genetic_code, to_stop)
            result = '%0.2f' % sequence.molecular_weight()
            data = [k, result]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def aromaticity(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the aromaticity index by each FASTA record. 
    This is calculated by BioPython library using method described in 
    Lobry JR, Gautier C. 1994. Hydrophobicity, expressivity and 
    aromaticity are the major trends of amino-acid usage in 999 
    Escherichia coli chromosome-encoded genes, Nucleic Acids Research 
    22(15):3174-3180. https://doi.org/10.1093/nar/22.15.3174

    Usage:

        python seqproperties.py aromaticity --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <aromaticity index>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            result = '%0.6f' % sequence.aromaticity()
            data = [k, result]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def instability(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the instability index by each FASTA record. 
    This is calculated by BioPython library using method described in 
    Guruprasad et al. 1990. Correlation between stability of a protein 
    and its dipeptide composition: a novel approach for predicting in 
    vivo stability of a protein from its primary sequence, Protein 
    Engineering, Design and Selection 4(2):155-161.

    Usage:

        python seqproperties.py instability --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <instability index>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            result = '%0.3f' % sequence.instability_index()
            data = [k, result]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def isoelectric(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the isoelectric point (pI) by each FASTA 
    record. This is calculated by BioPython library.

    Usage:

        python seqproperties.py isoelectric --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <isoelectric point>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            result = '%0.2f' % sequence.isoelectric_point()
            data = [k, result]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def secondaryStructure(fastafile, molecule, genetic_code=1, 
                       to_stop=True):
    '''!
    Function to calculate the secondary structure fractions by each 
    FASTA record. This is calculated by BioPython library.

    Usage:

        python seqproperties.py secstruct --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <helix fraction> : <turn fraction> : <sheet fraction>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        result = sequence.secondary_structure_fraction()
        helix = '%0.4f' % result[0]
        turn = '%0.4f' % result[1]
        sheet = '%0.4f' % result[2]
        data = [k, helix, turn, sheet]
        data = ' : '.join([str(x) for x in data])
        print(data)

def gravy(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the hydropathy, also known as GRAVY (Grand 
    Average of Hydropathy), by each FASTA record. This is calculated 
    by BioPython library using method described in Kyte J and 
    Doolittle RF. 1982. A simple method for displaying the hydropathic 
    character of a protein. J. Mol. Biol. 157, 105-132.

    Usage:

        python seqproperties.py gravy --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <GRAVY value>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            result = '%0.6f' % sequence.gravy()
            data = [k, result]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def flexibility(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to calculate the flexibility by each FASTA record. This 
    is calculated by BioPython library using method described in 
    Vihinen, et al. 1994. Accuracy of protein flexibility predictions. 
    Proteins, 19: 141-149.

    Usage:

        python seqproperties.py flexibility --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : {<flexibility value>}

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            data = [k] + sequence.flexibility()
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def extinction_coefficient(fastafile, molecule, genetic_code=1, 
                           to_stop=True):
    '''!
    Function to calculate the molar extinction coefficient by each 
    FASTA record. This is calculated by BioPython library.

    Usage:

        python seqproperties.py extinction --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

    Options for genetic_code and to_stop are needed if molecule type 
    is not peptide, as these options are needed for translation. The 
    output will be in the format of

        <sequence ID> : <extinction coefficient assuming reduced cysteine> : <extinction coefficient assuming non-reduced cysteine>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences (requires transcription and translation), and 'RNA' 
    for RNA sequence (requires translation).
    @param genetic_code Integer: Genetic code number to be used for 
    translation. Default = 1 (Standard Code). For more information, 
    see <https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi>
    @param to_stop Boolean: Flag to stop translation when first stop 
    codon is encountered. Default = True.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    for k in o.seqNN:
        sequence = _toPeptide(str(o.seqNN[k][0]), molecule, 
                              genetic_code, to_stop)
        try:
            result = sequence.molar_extinction_coefficient()
            data = [k, '%0.2f' % result[0], '%0.2f' % result[1]]
            data = ' : '.join([str(x) for x in data])
            print(data)
        except ZeroDivisionError:
            data = ' : '.join([str(k), 'undefined'])
            print(data)
        except KeyError:
            data = ' : '.join([str(k), 'KeyError'])
            print(data)
        except IndexError:
            data = ' : '.join([str(k), 'IndexError'])
            print(data)
        except:
            data = ' : '.join([str(k), 'Error'])
            print(data)

def _dictionaryGenerator(sequence, n, suffix=''):
    '''!
    Private method - Generates a dictionary of n-gram (as keys) 
    based on sequence.

    @param sequence List: Character sequence for n-gram generation.
    @param n Integer: Size of n-gram. If n=2, bigram will be generated.
    @return: Dictionary of n-grams where the keys will be the n-gram 
    and values will be zero.
    '''
    import itertools
    seqD = {}
    n = int(n) - len(suffix)
    for seq in itertools.product(sequence, repeat=n):
        seq = suffix.join(seq)
        seqD[seq] = 0
    return seqD

def nGram(fastafile, molecule, n):
    '''!
    Function to process n-grams by each FASTA record.

    Usage:

        python seqproperties.py ngram --fastafile=<FASTA file path> --molecule=<molecule type> --n=2

    The output will be in the format of:

        <sequence ID> : <list of n-gram counts> : <list of n-gram identities>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences, and 'RNA' for RNA sequence.
    @param n Integer: Size of n-gram. If n=2, bigram will be generated.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    if molecule == 'DNA':
        sequence = ['A', 'T', 'G', 'C', '*']
    elif molecule == 'RNA':
        sequence = ['A', 'U', 'G', 'C', '*']
    elif molecule == 'peptide':
        sequence = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                    'T', 'V', 'W', 'Y', '*']
    for k in o.seqNN:
        seq = o.seqNN[k][0]
        seqD = _dictionaryGenerator(sequence, n)
        for i in range(len(seq)-n):
            s = seq[i:i+n]
            seqD[s] = seqD[s] + 1
        table = list(seqD.keys())
        table.sort()
        result = ' : '.join([str(seqD[t]) for t in table])
        table = ' : '.join(table)
        print('%s : %s : %s' % (k, result, table))

def hasReverse(fastafile, molecule, min, max, suffix=''):
    '''!
    Function to process each FASTA record for the presence of 
    a sub-sequence and its reverse. For example, this is to see 
    if a DNA sequence has both GATCTA and ATCTAG in its sequence.

    Usage:

        python seqproperties.py reverse --fastafile=<FASTA file path> --molecule=<molecule type> --suffix=<substring to start sequence with> --min=3 --max=5

    The output will be in the format of:

        <sequence ID> : <length of reverse> : <sequence> : <reversed sequence>

    @param fastafile String: Path to the FASTA file to be processed.
    @param molecule String: Defines the type of molecule. Three 
    options are allowed: 'peptide' for amino acid sequences, 'DNA' for 
    DNA sequences, and 'RNA' for RNA sequence.
    @param min Integer: Minimum length of sub-sequence (including 
    suffix).
    @param max Integer: Maximum length of sub-sequence (including 
    suffix).
    @param suffix String: Defining the starting portion of the 
    sub-sequence - this is to enable the search for longer sub-
    sequence without running out of memory. Default = ''.
    '''
    def _generateReverseSequence(k, n, nonEmpty):
        for item in nonEmpty:
            rItem = ''.join(reversed(item))
            if rItem in nonEmpty:
                print('%s : %s : %s : %s' % (k, n, item, rItem))
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    if molecule == 'DNA':
        sequence = ['A', 'T', 'G', 'C', '*']
    elif molecule == 'RNA':
        sequence = ['A', 'U', 'G', 'C', '*']
    elif molecule == 'peptide':
        sequence = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 
                    'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 
                    'T', 'V', 'W', 'Y', '*']
    for k in o.seqNN:
        for n in range(min, max+1):
            seq = o.seqNN[k][0]
            seqD = _dictionaryGenerator(sequence, n, suffix)
            for i in range(len(seq)-n):
                s = seq[i:i+n]
                seqD[s] = seqD[s] + 1
            nonEmpty = [k1 for k1 in seqD if seqD[k1] > 0]
            _generateReverseSequence(k, n, nonEmpty)

def _readDiPFreq(datafile, separator=','):
    '''!
    Private function to read and process dipeptide frequency data 
    file to dictionary. The dipeptide frequency file should be in 
    the format of:

        <dipeptide>, <count>

    @param datafile String: Path to the dipeptide frequency file 
    (usually as comma-delimited file) to process.
    @param separator String: Separator in dipeptide frequency file. 
    Default=','.
    @return: Dictionary of dipeptides as keys and frequencies as 
    values.
    '''
    separator = str(separator)
    dataTable = {}
    data = open(datafile, 'r').readlines()
    if header:
        data = data[1:]
    data = [x[:-1].split(separator) for x in data]
    data = [[str(x[0].strip()), int(x[1].strip())] 
            for x in data]
    for freq in data:
        dataTable[freq[0]] = freq[1]
    return dataTable

def asymmetricFrequency(datafile, separator=',', header=False):
    '''!
    Function to process dipeptide frequency data file to asymmetric 
    frequency or C190 (Carugo O. 2013. Frequency of dipeptides and 
    antidipeptides. Computational and Structural Biotechnology 
    Journal 8, e201308001. doi:10.5936/csbj.201308001).

    Usage:

        python seqproperties.py asymfreq --datafile=<CSV file to process> --separator=<separator> --header=True

    The dipeptide frequency file should be in the format of:

        <dipeptide>, <count>

    The output will be in the format of:

        <dipeptide> : <antidipeptide> : <C190 score>

    @param datafile String: Path to the dipeptide frequency file 
    (usually as comma-delimited file) to process.
    @param separator String: Separator in dipeptide frequency file. 
    Default=','.
    @param header Boolean: Flag to indicate header row in dipeptide 
    frequency file. True, if the first row in the dipeptide 
    frequency file is header row. Default=False.
    '''
    dataTable = _readDiPFreq(datafile, separator)
    for seq in dataTable:
        dipeptideF = dataTable[seq]
        antiseq = seq[1] + seq[0]
        antidipeptideF = dataTable[antiseq]
        numerator = abs(dipeptideF - antidipeptideF)
        denominator = (dipeptideF + antidipeptideF) / 2
        result = numerator / denominator
        print('%s : %s : %s' % (seq, antiseq, result))

def propensity(datafile, separator=',', header=False):
    '''!
    Function to process dipeptide frequency data file to propensity 
    (Carugo O. 2013. Frequency of dipeptides and antidipeptides. 
    Computational and Structural Biotechnology Journal 8, 
    e201308001. doi:10.5936/csbj.201308001).

    Usage:

        python seqproperties.py propensity --datafile=<CSV file to process> --separator=<separator> --header=True

    The dipeptide frequency file should be in the format of:

        <dipeptide>, <count>

    The output will be in the format of:

        <dipeptide> : <propensity score>

    @param datafile String: Path to the dipeptide frequency file 
    (usually as comma-delimited file) to process.
    @param separator String: Separator in dipeptide frequency file. 
    Default=','.
    @param header Boolean: Flag to indicate header row in dipeptide 
    frequency file. True, if the first row in the dipeptide 
    frequency file is header row. Default=False.
    '''
    dataTable = _readDiPFreq(datafile, separator)
    nXX = 0
    for seq in dataTable:
        nXX = nXX + int(dataTable[seq])
    for seq in dataTable:
        nAB = dataTable[seq]
        nAX = 0
        for seq2 in dataTable:
            if seq[0] == seq2[0]:
                nAX = nAX + dataTable[seq2]
        nXB = 0
        for seq2 in dataTable:
            if seq[1] == seq2[1]:
                nXB = nXB + dataTable[seq2]
        result = (nAB / nXB) / (nAX / nXX)
        print('%s : %s' % (seq, result))

def pairwise_alignment(fastafile, algorithm='local'):
    '''!
    Function to take a FASTA file and calculate pairwise alignments 
    between all the sequences in the file.

    Usage:

        python seqproperties.py palign --fastafile=<FASTA file path> --algorithm=local

    The output will be in the format of

        <count> : <alignment score> : <sequence ID 1> : <sequence ID 2>

    where 
        - count is the numeric running order
        - alignment score is the calculated pairwise alignment 
        score
        - sequence ID 1 and 2 are the sequence IDs of the 2 FASTA 
        records used for pairwise alignment

    @param fastafile String: Path to the FASTA file to be processed.
    @param algorithm String: Type of pairwise alignment algorithm to 
    use. Allowable values are 'local' (Smith-Waterman algorithm) 
    and 'global' (Needleman-Wunsch algorithm). Default = local.
    '''
    o = CodonUsageBias()
    o.addSequencesFromFasta(fastafile)
    aligner = Align.PairwiseAligner()
    aligner.mode = str(algorithm)
    print(aligner)
    count = 1
    reduced_set = [k for k in o.seqNN]
    for k in o.seqNN:
        reduced_set = [k1 for k1 in reduced_set if k1 != k]
        for k1 in reduced_set:
            sequenceA = str(o.seqNN[k][0])
            sequenceB = str(o.seqNN[k1][0])
            score = aligner.score(sequenceA, sequenceB)
            print('%s : %s : %s : %s' % (str(count), str(score), 
                                         str(k), str(k1)))
            count = count + 1

def pairwise_alignment2(queryfile, dbfile, 
                        outfmt='summarize', algorithm='local'):
    '''!
    Function to take 2 FASTA files (a query FASTA file and a database 
    FASTA file) and calculate pairwise alignments in FASTA file to 
    database FASTA file.

    Usage:

        python seqproperties.py palign2 --queryfile=<FASTA file path> -dbfile=<FASTA file path> --algorithm=local --outfmt=summarize

    The output will be in the format of

        <count> : <alignment score> : <sequence ID from query> : <sequence ID from database>

    for full output format or the following for summarized output

        <count> : <minimum alignment score> : <average alignment score> : <standard deviation of alignment score> : <maximum alignment score> : <sequence ID from query>

    where 
        - count is the numeric running order
        - minimum, average, standard deviation, and maximum alignment 
        scores are calculated from pairwise alignment scores
        - sequence ID(s) is/are the sequence ID(s) of the FASTA 
        record(s) used for pairwise alignment

    @param queryfile String: Path to the FASTA file to be used as query.
    @param dbfile String: Path to the FASTA file to be used as database.
    @param outfmt String: Output format. Allowable values are 'full' 
    (one result line per comparison) and 'summarize' (one result line 
    per query FASTA record). Default = summarize.
    @param algorithm String: Type of pairwise alignment algorithm to 
    use. Allowable values are 'local' (Smith-Waterman algorithm) 
    and 'global' (Needleman-Wunsch algorithm). Default = local.
    '''
    q = CodonUsageBias()
    q.addSequencesFromFasta(queryfile)
    db = CodonUsageBias()
    db.addSequencesFromFasta(dbfile)
    aligner = Align.PairwiseAligner()
    aligner.mode = str(algorithm)
    print(aligner)
    count = 1
    if outfmt == 'full':
        print(' : '.join(['Count', 'Score', 'QuerySeqID', 
                          'DatabaseSeqID']))
    elif outfmt == 'summarize':
        print(' : '.join(['Count', 'Minimum Score', 'Average Score', 
                          'SD Score', 'Maximum Score', 'QuerySeqID', 
                          'DatabaseSeqID']))
    for qk in q.seqNN:
        querySeq = str(q.seqNN[qk][0])
        if outfmt == 'full':
            for dbk in db.seqNN:
                dbSeq = str(db.seqNN[dbk][0])
                score = aligner.score(querySeq, dbSeq)
                print('%s : %s : %s : %s' % (str(count), str(score), 
                                             str(qk), str(dbk)))
                count = count + 1
        elif outfmt == 'summarize':
            scores = [aligner.score(querySeq, str(db.seqNN[dbk][0])) 
                      for dbk in db.seqNN]
            min_score = min(scores)
            max_score = max(scores)
            avg_score = sum(scores) / len(scores)
            sd_score = [(s-avg_score) ** 2 for s in scores]
            sd_score = sum(sd_score) / len(scores)
            sd_score = sd_score ** 0.5
            print('%s : %s : %s : %s : %s : %s' % \
                  (str(count), str(min_score), str(avg_score),
                    str(sd_score), str(max_score), str(qk)))
            count = count + 1

def findORF(fastafile, min_length=33, max_length=105000, outfmt="CSV", 
            start_codons="TTG,CTG,ATG", stop_codons="TAA,TAG,TGA"):
    '''!
    Function to find open reading frames (ORF) for each FASTA record 
    in a given FASTA file. An ORF is basically computed as a stretch 
    of sequence flanked by a start and stop codon.

    Usage:

        python seqproperties.py orf --start_codons="TTG,CTG,ATG" --stop_codons="TAA,TAG,TGA" --fastafile=<fasta file path> --min_length=33 --max_length=105000 --outfmt=CSV-NS

    The CSV output will be in the format of:

        <count> : <sequence ID> : <start position> : <stop position> : <strand> : <length of ORF> : [<sequence of ORF>]

    The description line for FASTA output will be:

        <count>|<sequence ID>|<start position>|<stop position>|<strand>|<length of ORF>

    @param fastafile String: Path to the FASTA file to be processed.
    @param min_length Integer: Minimum length of ORF. Default = 33.
    @param max_length Integer: Maximum length of ORF. Default = 105000.
    @param output String: Defines type of output. Allowable types are 
    "CSV" (comma-separated file with sequence), "CSV-NS" (comma-separated 
    file without sequence) or "FASTA" (FASTA format). Default = "CSV".
    @param start_codons String: A comma-delimited list to represent 
    start codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no start codons within the sequence. Secondly, it can be used 
    to cap the start of a newly generated sequence.Default = 
    "TTG,CTG,ATG" (PMID 9278503).
    @param stop_codons String: A comma-delimited list to represent 
    stop codons. This is used for two purposes. Firstly, it can be 
    used to check the randomly generated sequence to ensure that there 
    is no stop codons within the sequence. Secondly, it can be used 
    to cap the stop of a newly generated sequence.Default = 
    "TAA,TAG,TGA".
    '''
    q = CodonUsageBias()
    q.addSequencesFromFasta(fastafile)
    if type(start_codons) is str:
        start_codons = start_codons.strip()
        start_codons = [x.strip() for x in start_codons.split(",")]
    elif type(start_codons) is tuple or type(start_codons) is list:
        start_codons = [x.strip() for x in start_codons]
    if type(stop_codons) is str:
        stop_codons = stop_codons.strip()
        stop_codons = [x.strip() for x in stop_codons.split(",")]
    elif type(stop_codons) is tuple or type(stop_codons) is list:
        stop_codons = [x.strip() for x in stop_codons]
    def find_all(string, substring):
        start = 0
        while True:
            start = string.find(substring, start)
            if start == -1: return
            yield start
            start = start + len(substring)
    def chunkstring(string, length):
        return (string[0+i:length+i] 
                for i in range(0, len(string), length))
    def find_coordinates(seq, start, stop_codons, 
                         min_length, max_length):
        s = list(chunkstring(seq[start:], 3))
        stop_location = [index for index in range(len(s)) 
                            if s[index] in stop_codons]
        if len(stop_location) > 0:
            stop_location = min(stop_location) * 3
            if stop_location >= min_length and \
                stop_location <= max_length:
                stop_location = stop_location + start + 3
                return (start, stop_location)
    count = 1
    if outfmt.upper() == 'CSV':
        print("Count : SequenceID : Start : Stop : Strand : Length : Sequence")
    elif outfmt.upper() == 'CSV-NS':
        print("Count : SequenceID : Start : Stop : Strand : Length")
    for k in q.seqNN:
        seq = str(q.seqNN[k][0])
        start_locations = [list(find_all(seq, start)) 
                           for start in start_codons]
        start_locations = [item for sublist in start_locations 
                            for item in sublist]
        coordinates = [find_coordinates(seq, start, stop_codons, 
                                        min_length, max_length)
                       for start in start_locations]
        coordinates = [coord for coord in coordinates 
                        if coord != None]
        for coord in coordinates:
            if outfmt.upper() == 'CSV':
                print("%s : %s : %s : %s : forward : %s : %s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0]), seq[coord[0]:coord[1]]))
            elif outfmt.upper() == 'CSV-NS':
                print("%s : %s : %s : %s : forward : %s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0])))
            elif outfmt.upper() == 'FASTA':
                print("> %s|%s|%s|%s|forward|%s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0])))
                print(seq[coord[0]:coord[1]])
            count = count + 1
        rev_seq = str(Seq(seq, generic_dna).reverse_complement())
        start_locations = [list(find_all(rev_seq, start)) 
                           for start in start_codons]
        start_locations = [item for sublist in start_locations 
                            for item in sublist]
        coordinates = [find_coordinates(rev_seq, start, stop_codons, 
                                        min_length, max_length)
                       for start in start_locations]
        coordinates = [coord for coord in coordinates 
                        if coord != None]
        for coord in coordinates:
            if outfmt.upper() == 'CSV':
                print("%s : %s : %s : %s : reverse : %s : %s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0]), 
                     rev_seq[coord[0]:coord[1]]))
            elif outfmt.upper() == 'CSV-NS':
                print("%s : %s : %s : %s : reverse : %s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0])))
            elif outfmt.upper() == 'FASTA':
                print("> %s|%s|%s|%s|reverse|%s" % \
                    (str(count), k, str(coord[0]), str(coord[1]), 
                     str(coord[1]-coord[0])))
                print(seq[coord[0]:coord[1]])
            count = count + 1

def random_selection(fastafile, n=250, with_replacement=True, 
                     outfmt='fasta'):
    '''!
    Function to select a random set of sequences from a given FASTA 
    file.

    Usage:

        python seqproperties.py rselect --fastafile=<fasta file path> --n=250 --with_replacement=True --outfmt=fasta

    The linear output format will be:

        <count> : <sequence ID> : <sequence>

    @param fastafile String: Path to the FASTA file to be processed.
    @param n Integer: Number of sequences to select. Default = 250.
    @param with_replacement String: Flag to indicate whether duplicated 
    selection is allowed. Allowable options are "True" (no duplicates 
    allowed) or "False" (duplicates allowed). Default = "True".
    @param outfmt String: Type of output. Allowable options are "linear" 
    (ID line and sequence in the same line) or "fasta" (FASTA format). 
    Default = 'fasta'
    '''
    q = CodonUsageBias()
    q.addSequencesFromFasta(fastafile)
    selection = []
    while len(selection) < int(n):
        s = random.sample(list(q.seqNN), k=1)[0]
        if str(with_replacement) == "True" and (s not in selection):
            selection.append((s, q.seqNN[s]))
        else:
            selection.append((s, q.seqNN[s]))
 
    count = 1
    for s in selection:
        if outfmt.lower() == 'linear':
            print("%s : %s : %s" % (str(count), s[0], s[1][0]))
        elif outfmt.lower() == 'fasta':
            print("> %s" % s[0])
            print(s[1][0])
        count = count + 1

def pointMutationOverGenerations(organisms=100, length=1000, bases="DNA", 
                                 mutations=10, mutation_rate=-1,
                                 algorithm="local", 
                                 generations=100, tests=100):
    """!
    Function to perform naive simulation of a population of sequences 
    over a number of generations and sample the sequence diversity at 
    each generation. 

    The simulation is based on the following: (a) an initial population 
    of identical sequences, (b) population size and sequence length do 
    not change, (c) number of bases/residues or probability of bases/
    residues to mutate does not change (probability, if more than zero, 
    will take precedence over number), (d) equal chance of mutation 
    per base/residue, (e) no selection or mating process, and (f) 
    random sampling of organisms to test for sequence diversity with 
    selfing and resampling possible.

    Usage:

        python seqproperties.py pmog --organisms=100 --length=1000 --bases=DNA --mutations=10 --mutation_rate=-1 --algorithm=local --generations=100 --tests=100

    The CSV output will be in the format of:

        <Generation>, {<Pairwise alignment score>}

    @param organisms Integer: Number of organisms. Default = 100
    @param length Integer: Sequence length per organism. Default = 1000
    @param bases String: Defines the type of bases, which can be 'peptide' 
    for amino acid sequences, 'DNA' for DNA sequences, and 'RNA' for RNA 
    sequence, or self-defined bases. Default = local
    @param mutations Integer: Number of bases/residues to mutate per 
    organism per generation. This will only be used if mutation_rate < 
    0. Default = 10
    @param mutation_rate float: Mutation rate where 0.01 or 1e-2 means 
    mutation rate is 1 in 100 or 1% chances of mutation per base/residue 
    per generation. If mutation_rate > 0, this will take precedence 
    over mutation parameter. Default = -1
    @param algorithm String: Type of pairwise alignment algorithm to 
    use. Allowable values are 'local' (Smith-Waterman algorithm) 
    and 'global' (Needleman-Wunsch algorithm). Default = local.
    @param generations Integer: Number of generations to simulate. 
    Default = 1000
    @param tests Integer: Number of pairwise alignments per generation. 
    Default = 100
    """
    if bases == "DNA":
        bases = [x for x in "ATGC"]
    elif bases == "RNA":
        bases = [x for x in "AUGC"]
    elif bases == "peptide":
        bases = [x for x in "ACDEFGHIKLMNPQRSTVWY"]
    else:
        bases = [x for x in bases]
    if algorithm == "local":
        aligner = Align.PairwiseAligner()
        aligner.mode = "local"
    elif algorithm == "global":
        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
    seq = [''.join([random.choice(bases) for i in range(int(length))])] * int(organisms)
    score_header = ','.join(["Score_" + str(i+1) for i in range(tests)])
    print("Generation,%s" % score_header)
    mutation_rate = float(mutation_rate)
    def mutate(s, mutation_rate, mutations):
        s = [base for base in s]
        if mutation_rate < 0:
            for m in range(int(mutations)):
                s[random.randint(0, len(s)-1)] = random.choice(bases)
        else:
            for position in range(len(s)):
                if random.random() <= mutation_rate:
                    s[position] = random.choice(bases)
        return ''.join(s)
    for gen in range(int(generations)+1):
        score = [str(aligner.score(random.choice(seq), random.choice(seq))) 
                 for test in range(int(tests))]
        print("%s,%s" % (str(gen), ",".join(score)))
        seq = [mutate(s, mutation_rate, mutations) for s in seq]

def extractFasta(fastafile, keyfile, outfile, match="start"):
    '''!
    Function to pull out / extract records from one FASTA file into 
    another.

    Usage:

        python seqproperties.py extractfasta --fastafile=<original FASTA file> --keyfile=<record list file> --outfile=<new output FASTA file> --match=<type of match>

    @param fastafile String: Path to the FASTA file to be processed - 
    the original FASTA file to extract data from. 
    @param keyfile String: Path to a file containing the list of 
    records - one description per line.
    @param outfile String: Path to the new FASTA file to be written.
    @param match String: Type of match. Allowable types are "start" 
    (matched as long as the start of the descriptor line is the same) 
    or "full" (full match between descriptor lines).
    '''
    infasta = CodonUsageBias()
    infasta.addSequencesFromFasta(fastafile)
    kfile = [x[:-1].strip() for x in open(keyfile, "r").readlines()]
    ofile = open(outfile, "w")
    ori_FASTA_count = str(len(infasta.seqNN))
    count = 0
    for key in kfile:
        for fastakey in infasta.seqNN:
            if match == "start":
                if fastakey.startswith(key):
                    print("Found: %s --> %s" % (str(key.strip()), 
                                                str(fastakey.strip())))
                    ofile.write(">" + str(fastakey.strip()) + '\n')
                    sequence = str(infasta.seqNN[fastakey][0])
                    ofile.write(str(sequence) + '\n')
                    count = count + 1
                    del infasta.seqNN[fastakey]
                    break
            elif match == "full":
                if fastakey.strip() == key.strip():
                    print("Found: %s --> %s" % (str(key.strip()), 
                                                str(fastakey.strip())))
                    ofile.write("> " + str(fastakey.strip()) + '\n')
                    sequence = str(infasta.seqNN[fastakey][0])
                    ofile.write(str(sequence) + '\n')
                    count = count + 1
                    del infasta.seqNN[fastakey]
                    break
    print("Number of records in original FASTA file: %s" % ori_FASTA_count)
    print("Number of items in Key File: %s" % str(len(kfile)))
    print("Number of records matched: %s" % str(count))
    ofile.close()

def differenceFasta(fastafileA, outfile, fastafileB=None, keyfile=None):
    '''!
    Function to pull out / extract records that are found in one FASTA 
    file (fastafileA) but not in the other FASTA file (fastafileB) or 
    descriptor listing (keyfile). The output FASTA file will consist of 
    (fastafileA - fastafileB) or (fastafileA - keyfile). The keyfile (if 
    given) will take precedence over fastafileB (which means fastafileB)

    Usage:

        python seqproperties.py difffasta --fastafileA=<original FASTA file> --fastafileB=<FASTA file to be subtracted> --keyfile=<record list file to be subtract> --outfile=<new output FASTA file>

    @param fastafileA String: Path to the FASTA file to be processed - 
    the original FASTA file to subtract data from. 
    @param fastafileB String: Path to the FASTA file containing record 
    to be subtracted.
    @param keyfile String: Path to a file containing the list of 
    records to be subtracted from fastafileA - one description per line.
    @param outfile String: Path to the new FASTA file to be written.
    '''
    fastaA = CodonUsageBias()
    fastaA.addSequencesFromFasta(fastafileA)
    fastaA_keys = list(fastaA.seqNN.keys())
    ofile = open(outfile, "w")
    if keyfile:
        kfile = [x[:-1].strip() 
                 for x in open(keyfile, "r").readlines()]
        subtracted_keys = [str(x) for x in fastaA_keys 
                            if x not in kfile]     
    else:
        fastaB = CodonUsageBias()
        fastaB.addSequencesFromFasta(fastafileB)
        fastaB_keys = list(fastaB.seqNN.keys())
        subtracted_keys = [str(x) for x in fastaA_keys 
                            if x not in fastaB_keys]
    for key in subtracted_keys:
        ofile.write(">" + key + '\n')
        sequence = str(fastaA.seqNN[key][0])
        ofile.write(str(sequence) + '\n')
    print("Number of records in FASTA file A: %s" % str(len(fastaA_keys)))
    if keyfile:
        print("Number of items in keyfile: %s" % str(len(kfile)))
    else:
        print("Number of records in FASTA file B: %s" % str(len(fastaB_keys)))
    print("Number of records in FASTA file A but not in FASTA file B: %s" % str(len(subtracted_keys)))
    ofile.close()

def cleanFasta(fastafile, outfile):
    '''!
    Function to clean out irrelevant data from a Fasta file - the resulting output Fasta file will only have Fasta records.

    Usage:

        python seqproperties.py cleanfasta --fastafile=<original FASTA file> --outfile=<new output FASTA file>

    @param fastafile String: Path to the FASTA file to be processed.
    @param outfile String: Path to the new FASTA file to be written.
    '''
    data = open(fastafile, "r").readlines()
    data = [x[:-1] for x in data]
    data = [x for x in data if x != ""]
    ofile = open(outfile, "w")
    fastaR = False
    for row in data:
        if row.startswith(">"):
            fastaR = True
            ofile.write(row + "\n")
        elif fastaR == True:
            ofile.write(row + "\n")
            fastaR = False
    ofile.close()

def fastaNGram(fastafile, n):
    '''!
    Function to generate starting N-grams of Fasta records.

    Usage:

        python seqproperties.py fastanname --fastafile=<FASTA file> --n=2

    @param fastafile String: Path to the FASTA file to be processed.
    @param n Integer: Size of n-gram. If n=2, bigram will be generated.
    '''
    infasta = CodonUsageBias()
    infasta.addSequencesFromFasta(fastafile)
    ngram = set([fastakey[:int(n)] for fastakey in infasta.seqNN])
    ngram = list(ngram)
    ngram.sort()
    print(", ".join(ngram))

def coexpression(expfile, method):
    '''!
    Function to generate gene co-expressions from expression data.

    Usage:

        python seqproperties.py coexp --expfile=<CSV file> --method=<coexpression method>

    @param expfile String: Path to the comma-separated value (CSV) file 
    containing gene expression data.
    @param method String: Co-expression measure. Allowable values are braycurtis (Bray and Curtis coefficient), cosine (Cosine coefficient) canberra (Canberra distance), euclidean (Euclidean distance), kendall (Kendall's tau), manhattan (Manhattan distance), pearson (Pearson's correlation), pointserial (Point biserial correlation), somer (Somer's D), spearman (Spearman's correlation), and tanimoto (Tanimoto coefficient).
    '''
    from scipy import stats
    from copads import objectdistance as d
    expData = {}
    for line in open(expfile, "r").readlines()[1:]:
        line = [x.strip() for x in line.split(',')]
        expData[line[0]] = [float(exp) for exp in line[1:]]
    idList1 = list(expData.keys())
    idList2 = list(expData.keys())
    count = 1
    for id1 in idList1:
        idList2 = [i for i in idList2 if i != id1]
        for id2 in idList2:
            if method == 'braycurtis': score = d.Bray_Curtis(expData[id1], expData[id2])
            if method == 'canberra': score = d.Canberra(expData[id1], expData[id2])
            if method == 'cosine': score = d.Cosine(expData[id1], expData[id2])
            if method == 'euclidean': score = d.Euclidean(expData[id1], expData[id2])
            if method == 'kendall': score = stats.kendalltau(expData[id1], expData[id2]).correlation
            if method == 'manhattan': score = d.Manhattan(expData[id1], expData[id2])
            if method == 'pearson': score = stats.pearsonr(expData[id1], expData[id2])[0]
            if method == 'pointbiserial': score = stats.pointbiserialr(expData[id1], expData[id2]).correlation
            if method == 'somer': score = stats.somersd(expData[id1], expData[id2]).statistic
            if method == 'spearman': score = stats.spearmanr(expData[id1], expData[id2]).correlation
            if method == 'tanimoto': score = d.Tanimoto(expData[id1], expData[id2])
            print('%s : %s : %s : %s' % (str(count), str(id1), 
                                         str(id2), str(score)))
            count = count + 1

def coexpression_randomization(expfile, method, n, replicate):
    '''!
    Function to generate randomized gene co-expressions from expression 
    data (for statistical testing).

    Usage:

        python seqproperties.py coexp_rand --expfile=<CSV file> --method=<coexpression method> --n=1000 --replicate=30

    @param expfile String: Path to the comma-separated value (CSV) file 
    containing gene expression data.
    @param method String: Co-expression measure. Allowable values are braycurtis (Bray and Curtis coefficient), cosine (Cosine coefficient) canberra (Canberra distance), euclidean (Euclidean distance), kendall (Kendall's tau), manhattan (Manhattan distance), pearson (Pearson's correlation), pointserial (Point biserial correlation), somer (Somer's D), spearman (Spearman's correlation), and tanimoto (Tanimoto coefficient).
    @param n Integer: Number of samples in each replicate.
    @param replicate Integer: Number of replicates.
    '''
    from scipy import stats
    from copads import objectdistance as d
    expData = {}
    for line in open(expfile, "r").readlines()[1:]:
        line = [x.strip() for x in line.split(',')]
        expData[line[0]] = [float(exp) for exp in line[1:]]
    idList = list(expData.keys())
    count = 1
    for rep in range(replicate):
        scores = []
        for i in range(n):
            d1 = expData[random.choice(idList)]
            d2 = expData[random.choice(idList)]
            if method == 'braycurtis': scores.append(d.Bray_Curtis(d1, d2))
            if method == 'canberra': scores.append(d.Canberra(d1, d2))
            if method == 'cosine': scores.append(d.Cosine(d1, d2))
            if method == 'euclidean': scores.append(d.Euclidean(d1, d2))
            if method == 'kendall': scores.append(stats.kendalltau(d1, d2).correlation)
            if method == 'manhattan': scores.append(d.Manhattan(d1, d2))
            if method == 'pearson': scores.append(stats.pearsonr(d1, d2)[0])
            if method == 'pointbiserial': scores.append(stats.pointbiserialr(d1, d2).correlation)
            if method == 'somer': scores.append(stats.somersd(d1, d2).statistic)
            if method == 'spearman': scores.append(stats.spearmanr(d1, d2).correlation)
            if method == 'tanimoto': scores.append(d.Tanimoto(d1, d2))
        mean_score = stats.describe(scores).mean
        print('%s : %s : %s' % (str(count), str(len(scores)), str(mean_score)))
        count = count + 1

def coexpression_filter(coexpfile, threshold=0, percentile=0, compare="above", absolute="yes", separator=":"):
    '''!
    Function to filter gene co-expressions file against a threshold.

    Usage:

        python seqproperties.py coexp_filter --compare=above --absolute=yes --separator=: --threshold=0.7 --coexpfile=<co-expression file>
        python seqproperties.py coexp_filter --compare=above --absolute=no --separator=: --percentile=0.95 --coexpfile=<co-expression file>

    @param coexpfile String: Path to file containing gene co-expressions.
    @param threshold Float: Threshold value, and will only be used if percentile is 0 or not given. Default = 0
    @param percentile Float: Percentile to filter, and takes precedence over threshold. Default = 0
    @param compare String: Type of comparison. Allowable types are "above" (filter co-expressions above the threshold) and "below" (filter co-expressions below the threshold). Default = above
    @param absolute String: Flag to take absolute values of co-expressions. Allowable types are "yes" (convert co-expression to absolute co-expression) and "no" (do not convert co-expression to absolute co-expression). Default = yes
    @param separator String: Separator in gene co-expression file. Default = :
    '''
    if percentile > 0:
        threshold = open(coexpfile, "r").readlines()
        threshold = [x[:-1].split(separator)[-1].strip() for x in threshold]
        if absolute == "yes": [abs(float(x)) for x in threshold]
        threshold.sort()
        index = int(float(percentile) * len(threshold))
        threshold = float(threshold[index])
    else:
        threshold = float(threshold)
    with open(coexpfile, "r") as f:
        for line in f:
            line = [x.strip() for x in line[:-1].split(separator)]
            line[-1] = float(line[-1])
            if absolute == "yes": line[-1] = abs(line[-1])
            if compare == "above":
                if line[-1] >= threshold:
                    line = [str(x) for x in line]
                    print(" : ".join(line))
            if compare == "below":
                if line[-1] <= threshold:
                    line = [str(x) for x in line]
                    print(" : ".join(line))

def coexpression_compare(coexpfile, truthfile, separator=":"):
    '''!
    Function to compare gene co-expressions file against model answer (truth file).

    Usage:

        python seqproperties.py coexp_compare --separator=: --coexpfile=<co-expression file> --truthfile=<actual protein-protein interaction file>

    Both gene co-expressions file and model answer (truth file) must have the following format: 

        <count> : <ID1> : <ID2> : ....

    where ":" is the separator; <ID1> and <ID2> pair represents significant gene co-expression or actual protein-protein interaction. This function only uses 3 columns.

    @param coexpfile String: Path to file containing gene co-expressions.
    @param truthfile String: Path to file containing protein-protein interactions for other truth file.
    @param separator String: Separator in gene co-expression file. Default = :
    '''
    truth = {}
    truthdata = [x[:-1] for x in open(truthfile).readlines()]
    truthdata = [[str(x.split(separator)[1].strip()), 
                  str(x.split(separator)[2].strip())] 
                for x in truthdata]
    for line in truthdata:
        truth[":".join(line)] = None
    true_positive = 0
    false_positive = 0
    false_negative = 0
    coexp = {}
    coexpdata = [x[:-1] for x in open(coexpfile).readlines()]
    coexpdata = [[str(x.split(separator)[1].strip()), 
                  str(x.split(separator)[2].strip())] 
                for x in coexpdata]
    for line in coexpdata:
        coexp[":".join(line)] = None
        if (":".join([line[0], line[1]]) in truth) or \
            (":".join([line[1], line[0]]) in truth): 
            true_positive = true_positive + 1
        elif (":".join([line[0], line[1]]) not in truth) and \
            (":".join([line[1], line[0]]) not in truth): 
            false_positive = false_positive + 1
    for x in truth:
        x = x.split(separator)
        if (":".join([x[0], x[1]]) not in coexp) and \
            (":".join([x[1], x[0]]) not in coexp):
            false_negative = false_negative + 1
    print("True positive = %s" % str(true_positive))
    print("False positive = %s" % str(false_positive))
    print("False negative = %s" % str(false_negative))
    precision = true_positive / (true_positive + false_positive)
    recall = true_positive / (true_positive + false_negative)
    F = (2 * precision * recall) / (precision + recall)
    CSI = true_positive / (true_positive + false_positive + false_negative)
    FM = (precision * recall) ** 0.5
    print("Precision = %.5f" % precision)
    print("Recall = %.5f" % recall)
    print("F1-Score = %.5f" % F)
    print("Critical Success Index = %.5f" % CSI)
    print("Fowlkes–Mallows Index = %.5f" % FM)

def checkCSV(file, separator=","):
    '''!
    Function to check for potentially problematic data rows in CSV file.

    Usage:

        python seqproperties.py checkcsv --file=<CSV file> --separator=,

    @param file String: Path to CSV file.
    @param separator String: Separator in CSV file. Default = ,
    '''
    import pandas as pd
    try:
        data = pd.read_csv(file, sep=separator, on_bad_lines=False)
    except:
        data = pd.read_csv(file, sep=separator, error_bad_lines=False, 
                           warn_bad_lines=True)

def sample_ClusterScan(file, filetype="excel", sheet_name=None, 
                       usecols=None, samplesbyrow=True, 
                       min_clusters=2, max_clusters=100):
    '''!
    Function to calculate Davies-Bouldin Score, Calinski and Harabasz 
    Score, and Silhouette Score; for the number of clusters in order 
    to determine the suitable number of clusters for K-means clustering 
    of the data.

    Usage:

        python seqproperties.py clusterscan --file=<data file> --filetype=<type of file> --sheet_name=<name of Excel sheet> --usecols=<columns to use> --samplesbyrow=<whether sample are by rows> --min_clusters=<minimum number of clusters> --max_clusters=<maximum number of clusters>

    For example,

        python seqproperties.py clusterscan --file=iAF692_fluxes.xlsx --filetype=excel --sheet_name=iAF692 --usecols=B:ZK --samplesbyrow=False

    @param file String: Path to data file.
    @param filetype String: Type of file. Allowable types are "csv" 
    (comma-separated file) and "excel" (Microsoft Excel workbook). 
    Default = excel
    @param sheet_name String: Name of sheet to analyse (only if filetype 
    = excel). Default = None
    @param usecols String: Columns to use. Default = None (all columns 
    will be used)
    @param samplesbyrow Boolean: Flag to determine if samples are by 
    rows. If False, it means samples are by columns. Default = True
    @param min_clusters Integer: Minimum number of clusters to scan. 
    Default = 2
    @param max_clusters Integer: Maximum number of clusters to scan. 
    Default = 100
    '''
    from sklearn.cluster import KMeans
    from sklearn.metrics import davies_bouldin_score
    from sklearn.metrics import calinski_harabasz_score
    from sklearn.metrics import silhouette_score
    import pandas as pd
    if filetype.lower() == "excel":
        df = pd.read_excel(file, sheet_name, usecols, engine="openpyxl")
    elif filetype.lower() == "csv":
        df = pd.read_csv(file, usecols)
    if samples.lower == "column":
        df = df.T     # columns are features, rows are samples
    print("Clusters, Davies-Bouldin Score, Calinski and Harabasz Score, Silhouette Score")
    for centre in range(int(min_clusters), int(max_clusters)+1):
        kmeans = KMeans(n_clusters=centre)
        model = kmeans.fit_predict(df)
        scores = [str(davies_bouldin_score(df, model)),
                  str(calinski_harabasz_score(df, model)),
                  str(silhouette_score(df, model))]
        print(centre, ",", ",".join(scores))

def sample_ClusterLabel(file, filetype="excel", sheet_name=None, 
                        usecols=None, samplesbyrow=True, 
                        clusters=10, label="cluster",
                        resultfile="result.csv"):
    '''!
    Function to cluster data from file by K-means clustering and append 
    the cluster label to the result file.

    Usage:

        python seqproperties.py clusterlabel --file=<data file> --filetype=<type of file> --sheet_name=<name of Excel sheet> --usecols=<columns to use> --samplesbyrow=<whether sample are by rows> --clusters=<number of clusters> --label=<cluster label> --resultfile=<result file name>

    For example,

        python seqproperties.py clusterlabel --file=iAF692_fluxes.xlsx --filetype=excel --sheet_name=iAF692 --usecols=B:ZK --samplesbyrow=False --clusters=10 --label=cluster --resultfile=iAF692_clustered.csv

    @param file String: Path to data file.
    @param filetype String: Type of file. Allowable types are "csv" 
    (comma-separated file) and "excel" (Microsoft Excel workbook). 
    Default = excel
    @param sheet_name String: Name of sheet to analyse (only if filetype 
    = excel). Default = None
    @param usecols String: Columns to use. Default = None (all columns 
    will be used)
    @param samplesbyrow Boolean: Flag to determine if samples are by 
    rows. If False, it means samples are by columns. Default = True
    @param clusters Integer: Number of clusters. Default = 10
    @param label String: Field name for cluster. Default = cluster
    @param resultfile String: Name of result file containing clusters. 
    Default = result.csv
    '''
    import pandas as pd
    from sklearn.cluster import KMeans
    if filetype.lower() == "excel":
        df = pd.read_excel(file, sheet_name, usecols, engine="openpyxl")
    elif filetype.lower() == "csv":
        df = pd.read_csv(file, usecols)
    if samples.lower == "column":
        df = df.T     # columns are features, rows are samples
    model = KMeans(n_clusters=int(clusters), random_state=0).fit(df)
    df[str(label)] = model.labels_
    df.to_csv(str(resultfile))

def overlap_statistics(file1, file2, separator, n_item, replicate=30):
    '''!
    Function to randomize 2 lists (given as files - file1 and file 2) 
    and generate the number of overlapping elements within the 2 randomized 
    list for a number of replicates. This can be used to test whether 
    the number of overlapping elements in 2 files is significantly 
    higher or lower than chance.

    Each file is in the format of <item 1><separator><item 2><separator>
    <item 3><separator> ... <separator><item N>

    Usage:

        python seqproperties.py overlap_stat --file1=overlap1.txt --file2=overlap2.txt --separator=: --n_item=2 --replicate=30

    @param file1 String: Path to data file 1.
    @param file2 String: Path to data file 2.
    @param separator String: Separator for items in the files. Default = :
    @param n_items Integer: Number of consecutive items to be considered as 
    an element.
    @param replicate Integer: Number of replicates.
    '''
    dataA = [x[:-1] for x in open(file1).readlines()]
    dataA = [[x.strip() for x in row.split(separator)] for row in dataA]
    dataA = [separator.join(row[:n_item]) for row in dataA]
    dataB = [x[:-1] for x in open(file2).readlines()]
    dataB = [[x.strip() for x in row.split(separator)] for row in dataB]
    dataB = [separator.join(row[:n_item]) for row in dataB]
    actual_overlap = sum([1 for x in dataA if x in dataB])
    print("Actual number of overlaps = " + str(actual_overlap))
    for i in range(replicate):
        combinelist = dataA + dataB
        random.shuffle(combinelist)
        listA = combinelist[:len(dataA)]
        listB = combinelist[len(dataA):]
        randomized_overlap = sum([1 for x in listA if x in listB])
        print("Randomized overlaps %s = %s" % (i+1, randomized_overlap))


if __name__ == '__main__':
    exposed_functions = {'a': percentA,
                         'aacount': aminoacidCount,
                         'ai': percentAi,
                         'aromaticity': aromaticity,
                         'asymfreq': asymmetricFrequency,
                         'checkcsv': checkCSV,
                         'cleanfasta': cleanFasta,
                         'clusterscan': sample_ClusterScan,
                         'clusterlabel': sample_ClusterLabel,
                         'codoncount': codonCount,
                         'coexp': coexpression,
                         'coexp_compare': coexpression_compare,
                         'coexp_filter': coexpression_filter,
                         'coexp_rand': coexpression_randomization,
                         'complement': complement,
                         'difffasta': differenceFasta,
                         'extinction': extinction_coefficient,
                         'extractfasta': extractFasta,
                         'fastanname': fastaNGram,
                         'count': genericCount,
                         'flexibility': flexibility,
                         'g': percentG,
                         'gc': percentGC,
                         'gci': percentGCi,
                         'gi': percentGi,
                         'gravy': gravy,
                         'instability': instability,
                         'isoelectric': isoelectric,
                         'mw': molecularWeight,
                         'ngram': nGram,
                         'nlength': nucleotideLength,
                         'orf': findORF,
                         'overlap_stat': overlap_statistics,
                         'palign': pairwise_alignment,
                         'palign2': pairwise_alignment2,
                         'plength': peptideLength,
                         'pmog': pointMutationOverGenerations,
                         'propensity': propensity,
                         'reverse': hasReverse,
                         'rselect': random_selection,
                         'secstruct': secondaryStructure,
                         'showDesc': sequenceDescriptions,
                         'showIDs': sequenceIDs,
                         'translate': translate}
    fire.Fire(exposed_functions)
