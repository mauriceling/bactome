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
import subprocess
import sys

try: 
    from Bio import Align
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna
    from Bio.Seq import Seq
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'biopython'])
    from Bio import Align
    from Bio import SeqIO
    from Bio.Alphabet import generic_dna
    from Bio.Alphabet import generic_rna
    from Bio.Seq import Seq
    from Bio.SeqUtils.ProtParam import ProteinAnalysis


try: 
    import fire
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'fire'])
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

        python seqproperties.py translate --fastafile=<FASTA file path> 
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

def aminoacidCount(fastafile, molecule, genetic_code=1, to_stop=True):
    '''!
    Function to translate each nucleotide sequence (by FASTA record) 
    and generate a frequency table of the amino acids.

    Usage:

        python seqproperties.py aacount --molecule=<molecule type> --genetic_code=<genetic code number> --to_stop=<Boolean flag> --fastafile=<FASTA file path>

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
    header = ['SequenceID'] + char_set
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

        python seqproperties.py codoncount --fastafile=<FASTA file path> 
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
    Function to calculate the molecular weight by each FASTA record.

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
    22(15):3174–3180. https://doi.org/10.1093/nar/22.15.3174

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
    Engineering, Design and Selection 4(2):155–161.

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

        <sequence ID> : <list of n-gram counts> : 
        <list of n-gram identities>

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

        <sequence ID> : <length of reverse> : <sequence> : 
        <reversed sequence>

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

        python seqproperties.py palign --queryfile=<FASTA file path> -dbfile=<FASTA file path> --algorithm=local --output=summarize

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
            sd_score = sum(numerator) / len(scores)
            sd_score = sd_score ** 0.5
            print('%s : %s : %s : %s : %s' % \
                  (str(count), str(min_score), str(avg_score),
                    str(sd_score), str(max_score), str(qk)))
            count = count + 1


if __name__ == '__main__':
    exposed_functions = {'a': percentA,
                         'aacount': aminoacidCount,
                         'ai': percentAi,
                         'aromaticity': aromaticity,
                         'asymfreq': asymmetricFrequency,
                         'codoncount': codonCount,
                         'count': genericCount,
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
                         'palign': pairwise_alignment,
                         'palign2': pairwise_alignment2,
                         'plength': peptideLength,
                         'propensity': propensity,
                         'reverse': hasReverse,
                         'secstruct': secondaryStructure,
                         'showDesc': sequenceDescriptions,
                         'showIDs': sequenceIDs,
                         'translate': translate}
    fire.Fire(exposed_functions)
