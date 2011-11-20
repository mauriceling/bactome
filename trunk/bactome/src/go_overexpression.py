"""
Extracting Gene List for Gene Ontology Over-representation.

This script only works with Python version < 3.0 due to print 
statements. 

Date created: 1st July 2011

Licence: Python Software Foundation License version 2

@see: Ling, MHT. 2011. Bactome II: Analyzing Gene List for Gene Ontology 
Over-Representation. The Python Papers Source Codes 3: 3.
"""

import os
import sys

def parse_GAF(gaf_file):
    """
    Reader for GO Annotation File - http://www.geneontology.org/
    GO.format.gaf-2_0.shtml
    
    @param gaf_file: file path for GO annotation file
    @returns: a tuple - (gaf, columns, gaf_ver, header) - 
        "gaf" is the annotation, 
        "columns" is the column header in the annotation, 
        "gaf_ver" is the GAF version (either 1.0 or 2.0), 
        "header" contains all the comments in the annotation file.
    
    @raise AttributeError: when there is no GAF version or if GAF 
    version is not 1.0 or 2.0
    """
    gaf = open(gaf_file, 'r').xreadlines()
    gaf = [x[:-1].strip() for x in gaf]
    header = [x for x in gaf if x.startswith('!')]
    gaf = [x for x in gaf if not x.startswith('!')]
    if '!gaf-version: 2.0' in header:
        gaf_ver = 2.0
    elif '!gaf-version: 1.0' in header:
        gaf_ver = 1.0
    else:
        raise AttributeError('File does not have GAF version')
    gaf = [x.split('\t') for x in gaf]
    gaf_16 = [x for x in gaf if len(x) == 16]
    gaf_17 = [x for x in gaf if len(x) == 17]
    gaf = [x for x in gaf if len(x) == 15]
    gaf = [(x[0].strip(), x[1].strip(), x[2].strip(), 
            [y.strip() for y in x[3].split('|')], x[4].strip(), 
            [y.strip() for y in x[5].split('|')], x[6].strip(), 
            [y.strip() for y in x[7].split('|')], x[8].strip(), 
            x[9].strip(), [y.strip() for y in x[10].split('|')], 
            x[11].strip(), [y.strip() for y in x[12].split('|')], 
            x[13].strip(), x[14].strip())
           for x in gaf]
    gaf_16 = [(x[0].strip(), x[1].strip(), x[2].strip(), 
               [y.strip() for y in x[3].split('|')], x[4].strip(), 
               [y.strip() for y in x[5].split('|')], x[6].strip(), 
               [y.strip() for y in x[7].split('|')], x[8].strip(), 
               x[9].strip(), [y.strip() for y in x[10].split('|')], 
               x[11].strip(), [y.strip() for y in x[12].split('|')], 
               x[13].strip(), x[14].strip(), 
               [y.strip() for y in x[15].split('|')])
              for x in gaf_16]
    gaf_17 = [(x[0].strip(), x[1].strip(), x[2].strip(), 
               [y.strip() for y in x[3].split('|')], x[4].strip(), 
               [y.strip() for y in x[5].split('|')], x[6].strip(), 
               [y.strip() for y in x[7].split('|')], x[8].strip(), 
               x[9].strip(), [y.strip() for y in x[10].split('|')], 
               x[11].strip(), [y.strip() for y in x[12].split('|')], 
               x[13].strip(), x[14].strip(), 
               [y.strip() for y in x[15].split('|')], 
               [y.strip() for y in x[16].split('|')])
              for x in gaf_17]
    gaf = gaf + gaf_16 + gaf_17
    if gaf_ver == 2.0:
        columns = ['DB', 'DB Object ID', 'DB Object Symbol', 
                   'Qualifier', 'GO ID', 'DB:Reference', 'Evidence Code', 
                   'With (or) From', 'Aspect', 'DB Object Name', 
                   'DB Object Synonym', 'DB Object Type', 'Taxon', 
                   'Date', 'Assigned By', 'Annotation Extension', 
                   'Gene Product Form ID']
    if gaf_ver == 1.0:        
        columns = ['DB', 'DB Object ID', 'DB Object Symbol', 
                   'Qualifier', 'GO ID', 'DB:Reference', 'Evidence Code', 
                   'With (or) From', 'Aspect', 'DB Object Name', 
                   'DB Object Synonym', 'DB Object Type', 'Taxon', 
                   'Date', 'Assigned By']
    return (gaf, columns, gaf_ver, header)

def process_association_file(assoc_file):
    """
    Extracts information from GO annotation file for over-representation
    analysis
    
    @param assoc_file: File path for GO annotation file
    
    @return: A tuple - (genome, gene_set, length of gene_set, 
        ontology_set, length of ontology_set) - 
        "genome" is a list of (<gene symbol>, <GO ID>) representing the 
        annotated genome
        "gene_set" is the set of genes in the genome
        "length of gene_set" is the number of genes in the genome
        "ontology_set" is the set of GO IDs used in the annotation
        "length of ontology_set" is the number of GO IDs used in the 
        annotation
    """
    genome = parse_GAF(assoc_file)[0]
    genome = [(x[2].lower(), x[4].lower())
              for x in genome
                  if len(x) > 5]
    gene_set = set([x[0] for x in genome])
    ontology_set = set([x[1] for x in genome])
    return (genome, gene_set, len(gene_set), ontology_set, 
            len(ontology_set))

def process_sample_file(samplefile, genome):
    """
    Processes the gene list file, the set of genes to be analyzed for 
    GO over-representation, and generate an annotated sample genome 
    (to represent the observed data in Chi Square test)
    
    @param samplefile: File path for the gene list for over-representation
        analysis - one gene per line
    @param genome: Genome to base expectation on. The genome is a list 
        of (<gene symbol>, <GO ID>)
    """
    sample_list = open(samplefile, 'r').readlines()
    sample_list = [x[:-1].lower() for x in sample_list]
    sample_genome = [x for x in genome
                     if x[0].lower() in sample_list]
    gene_set = set([x[0] for x in sample_genome])
    sample_ontology = set([x[1] for x in sample_genome])
    return(sample_list, sample_genome, gene_set, len(gene_set), 
           sample_ontology)

def find_ontology_by_gene(gene, g):
    """
    Find the set of GO IDs for a given gene
    
    @param gene: Gene symbol
    @param g: Genome to base expectation on. The genome is a list of 
        (<gene symbol>, <GO ID>)
        
    @return: A set of GO IDs
    """
    return set([x[1] for x in g if x[0] == gene.lower()])

def find_gene_by_ontology(ontology, g):
    """
    Find the set of genes by a given ontology ID
    
    @param ontology: Gene ontology ID
    @param g: Genome to base expectation on. The genome is a list of 
        (<gene symbol>, <GO ID>)
        
    @return: A set of gene symbols
    """
    return set([x[0] for x in g if x[1] == ontology.lower()])

def calculate_expectation(ontology, genome, genome_count):
    """
    Calculated the expectation (number of expected in Chi Square test) 
    for each ontology ID and the given genome
    
    @param ontology: Gene ontology ID to calculate expectation
    @param genome: Genome to base expectation on. The genome is a list
        of (<gene symbol>, <GO ID>)
    @param genome_count: Number of genes in the genome
    
    @return: Expectation as a float number.
    """
    gene_w_ontology = len(find_gene_by_ontology(ontology, genome))
    return float(gene_w_ontology) / float(genome_count)

def _write_headers(f, assoc_file, samplefile, outputfile,
                  genome_count, ontology_count, sample_count,
                  mapped_count, genes_not_mapped, genes_mapped,
                  test_count, display=True):
    """
    Writes the information header for the analysis - not to be used 
    externally
    """
    f.write(' '.join(['Gene Ontology association file:',
                      str(assoc_file), '\n']))
    f.write(' '.join(['Sample file:', str(samplefile), '\n']))
    f.write(' '.join(['Output file:', str(outputfile), '\n']))
    f.write(' '.join(['Number of genes in association file:',
                      str(genome_count), '\n']))
    f.write(' '.join(['Number of ontological terms in association file:',
                      str(ontology_count), '\n']))
    f.write(' '.join(['Number of genes in sample file:',
                      str(sample_count), '\n']))
    f.write(' '.join(['Number of genes in sample file mapped:',
                      str(mapped_count), '\n']))
    f.write(' '.join(['Genes in sample file mapped:',
                      str(genes_mapped), '\n']))                                
    if genes_not_mapped == '':
        genes_not_mapped = '(None - All genes are mapped)'
    f.write(' '.join(['Genes in sample file not mapped:',
                      str(genes_not_mapped), '\n']))
    f.write(' '.join(['Number of statistical tests:',
                      str(test_count), '\n']))
    if display:
        print ' '.join(['Gene Ontology association file:',
                      str(assoc_file)])
        print ' '.join(['Sample file:', str(samplefile)])
        print ' '.join(['Output file:', str(outputfile)])
        print ' '.join(['Number of genes in association file:',
                      str(genome_count)])
        print ' '.join(['Number of ontological terms in association file:',
                      str(ontology_count)])
        print ' '.join(['Number of genes in sample file:',
                      str(sample_count)])
        print ' '.join(['Number of genes in sample file mapped:',
                      str(mapped_count)])
        print ' '.join(['Genes in sample file mapped:',
                      str(genes_mapped)])
        print ' '.join(['Genes in sample file not mapped:',
                      str(genes_not_mapped)])
        print ' '.join(['Number of statistical tests:',
                      str(test_count)])

def main(assoc_file, samplefile, outputfile, display=True):
    """
    Runner for the GO over-representation analysis
    
    @param assoc_file: File path for GO annotation file
    @param samplefile: File path for the gene list for over-representation
        analysis - one gene per line
    @param outputfile: File path to write the analysis output
    @param display: Flag to display the results on console (default = True)
    """
    (genome, gene_set, genome_count, ontology_set, ontology_count) = \
             process_association_file(assoc_file)
    (sample_list, sample_genome, gene_set, sample_count, sample_ontology) = \
            process_sample_file(samplefile, genome)
    output_f = open(outputfile, 'w')
    genes_not_mapped = ', '.join([x for x in sample_list 
                                    if x not in gene_set])
    genes_mapped = ', '.join([x for x in sample_list 
                                    if x in gene_set])
    _write_headers(output_f, assoc_file, samplefile, outputfile,
                  genome_count, ontology_count, len(sample_list),
                  len(gene_set), genes_not_mapped, genes_mapped,
                  len(sample_ontology))
    output_f.write(','.join(['Term', 'Expected Count', '(Not) Expected',
                              'Actual Count', '(Not) Actual', '\n']))
    if display:
        print '\t'.join(['Term', 'Expected Count', '(Not) Expected',
                         'Actual Count', '(Not) Actual'])
    for ontology in sample_ontology:
        expectation = calculate_expectation(ontology, genome, genome_count)
        expected_count = expectation * sample_count
        expected_uncount = sample_count - expected_count
        actual_count = len(find_gene_by_ontology(ontology, sample_genome))
        actual_uncount = sample_count - actual_count
        output_f.write(','.join([ontology,
                                 str(expected_count),
                                 str(expected_uncount),
                                 str(actual_count),
                                 str(actual_uncount), '\n']))
        if display:
            print '\t'.join([ontology,
                             str(expected_count),
                             str(expected_uncount),
                             str(actual_count),
                             str(actual_uncount)])
    output_f.close()
        
def test():
    """
    Example execution using Escherichia coli GAF and Pseudomonas 
    aeruginosa GAF
    """
    main('gene_association.ecocyc', 'e_coli_gene.txt', 'e_coli_output.csv')
    main('gene_association.pseudocap', 'p_aeruginosa.txt', 
         'p_aeruginosa_output.csv')

if __name__ == '__main__':
    if len(sys.argv) != 4 and sys.argv[1] != 'test':
        print """Usage: python go_overexpression.py <association file> 
        <gene list file> <output file> \n
        where \n
        <association file> is the Gene Ontology annotation file for the 
        organism downloaded from http://www.geneontology.org/
        GO.downloads.annotations.shtml \n
        <gene list file> is the text file of gene list - one gene per line \n
        <output file> is the name of CSV file (comma-separated file) where 
        output is stored \n
        For example: python go_overexpression.py gene_association.ecocyc 
        e_coli_gene.txt e_coli_output.csv'"""
    elif sys.argv[1] == 'test':
        test()
    else:
        main(str(sys.argv[1]), str(sys.argv[2]), str(sys.argv[3]))