"""
Test script for src/go_overexpression.py

This script only works with Python version < 3.0 due to print 
statements. 

Date created: 1st July 2011

Licence: Python Software Foundation License version 2

@see: Ling, MHT. 2011. Bactome II: Analyzing Gene List for Gene Ontology 
Over-Representation. The Python Papers Source Codes 3: 3.
"""

import sys
import os
import unittest

sys.path.append(os.path.join(os.path.dirname(os.getcwd()), 'src'))

import go_overexpression as g

assocfile_eco = 'gene_association.ecocyc'
assocfile_pse = 'gene_association.pseudocap'
gene_eco = 'e_coli_gene.txt'
gene_pse = 'p_aeruginosa.txt'

class testGO(unittest.TestCase):
    def testVersion(self):
        self.assertTrue(g.parse_GAF(assocfile_eco)[2] == 2.0)
        self.assertTrue(g.parse_GAF(assocfile_pse)[2] == 2.0)
        
    def testAssoc(self):
        assoc1 = ('EcoCyc', '3-OXOACYL-ACP-SYNTHII-MONOMER', 'FabF',
                  [''], 'GO:0006633', ['GOA:interpro', 'GO_REF:0000002'],
                  'IEA', ['InterPro:IPR017568'], 'P', 'FabF', ['fabF',
                  'vtr', 'cvc', 'fabJ', 'vtrB', 'b1095', 'ECK1081'],
                  'protein', ['taxon:511145'], '20110110', 'InterPro')
        assoc2 = ('EcoCyc', 'AMINEOXID-MONOMER', 'TynA', [''], 'GO:0008131',
                  ['GOA:interpro', 'GO_REF:0000002'], 'IEA',
                  ['InterPro:IPR000269'], 'F', 'TynA', ['tynA', 'feaA',
                  'maoA', 'b1386', 'ECK1383'], 'protein', ['taxon:511145'],
                  '20110110', 'InterPro')
        assoc3 = ('PseudoCAP', 'PA4967', 'parE', [''], 'GO:0006259',
                  ['PMID:99169751'], 'IDA', [''], 'P',
                  'topoisomerase IV subunit B', [''], 'protein',
                  ['taxon:208964'], '20060206', 'PseudoCAP')
        assoc4 = ('PseudoCAP', 'PA4932', 'rplI', [''], 'GO:0044267',
                  ['PMID:96374488'], 'ISS', [''], 'P',
                  '50S ribosomal protein L9', [''], 'protein',
                  ['taxon:208964'], '20060206', 'PseudoCAP')
        self.assertTrue(assoc1 in g.parse_GAF(assocfile_eco)[0])
        self.assertTrue(assoc2 in g.parse_GAF(assocfile_eco)[0])
        self.assertTrue(assoc3 in g.parse_GAF(assocfile_pse)[0])
        self.assertTrue(assoc4 in g.parse_GAF(assocfile_pse)[0])

    def testOntologyByGene(self):
        eco_genome = g.process_association_file(assocfile_eco)[0]
        pse_genome = g.process_association_file(assocfile_pse)[0]
        onto1 = ['go:0006633', 'go:0005515', 'go:0005829', 'go:0008415',
                 'go:0016740', 'go:0003824', 'go:0009058', 'go:0016747',
                 'go:0008152', 'go:0008610']
        onto2 = ['go:0005737', 'go:0008652', 'go:0050661', 'go:0006561',
                 'go:0004735', 'go:0008152', 'go:0055114']
        onto3 = ['go:0009276', 'go:0009103']
        onto4 = ['go:0006091', 'go:0006810']
        self.assertEqual(list(g.find_ontology_by_gene('fabf', eco_genome)),
                         onto1)
        self.assertEqual(list(g.find_ontology_by_gene('proc', eco_genome)),
                         onto2)
        self.assertEqual(list(g.find_ontology_by_gene('rmd', pse_genome)),
                         onto3)
        self.assertEqual(list(g.find_ontology_by_gene('nosf', pse_genome)),
                         onto4)

    def testGeneByOntology(self):
        eco_genome = g.process_association_file(assocfile_eco)[0]
        pse_genome = g.process_association_file(assocfile_pse)[0]
        gene1 = ['acpp', 'ychm', 'fabr', 'accd', 'fabb', 'acca', 'fabd',
                 'accc', 'accb', 'acps', 'fabg', 'fabf', 'faba', 'fabz',
                 'plsx', 'lipoyl-acp', 'fabi', 'acph', 'fabh']
        gene2 = ['waac', 'wzm', 'glmu', 'lpxk', 'wzz', 'waaf', 'wbpm',
                 'wzy', 'wzt', 'wbpg', 'rmla', 'rmlb', 'wbpy', 'wbpx',
                 'rmd', 'glmm', 'wbpw']
        self.assertEqual(list(g.find_gene_by_ontology('go:0006633',
                                                      eco_genome)), gene1)
        self.assertEqual(list(g.find_gene_by_ontology('go:0009103',
                                                      pse_genome)), gene2)

    def testExpectation(self):
        eco_assoc = g.process_association_file('gene_association.ecocyc')
        pse_assoc = g.process_association_file('gene_association.pseudocap')
        self.assertAlmostEqual(g.calculate_expectation('go:0006633',
                                                       eco_assoc[0],
                                                       eco_assoc[2]),
                               0.00492994)
        self.assertAlmostEqual(g.calculate_expectation('go:0005515',
                                                       eco_assoc[0],
                                                       eco_assoc[2]),
                               0.23040996)
        self.assertAlmostEqual(g.calculate_expectation('go:0009276',
                                                       pse_assoc[0],
                                                       pse_assoc[2]),
                               0.07378129)
    def testReadSample(self):
        eco_assoc = g.process_association_file('gene_association.ecocyc')
        pse_assoc = g.process_association_file('gene_association.pseudocap')
        eco_list = ['mdog', 'dapa', 'crp', 'hslv', 'mrdb', 'fucu', 'yjgp',
                    'yigc', 'sun', 'gor', 'hflb', 'yqib', 'murg', 'yrbg',
                    'yejk', 'yfga', 'hflx', 'spot', 'holc', 'xerd', 'tolb',
                    'yhes', 'ntpa', 'yabb', 'lola', 'yggd', 'pnp', 'yrbb',
                    'rnc', 'xerc', 'rfaf', 'yigp', 'gyrb', 'nagc', 'nrdr',
                    'hemd', 'phet', 'frr', 'cls']
        pse_list = ['ppk', 'kdsa', 'dnaq', 'pcm', 'rpsk', 'thrh', 'thrb',
                    'zipa', 'rpsl', 'pcrv', 'rsal', 'msud', 'msue', 'flio',
                    'flip', 'fliq', 'flir', 'flhb', 'flha', 'aruc', 'exab',
                    'arnt', 'parc', 'pare', 'pa2018', 'pa2019', 'ftse', 'ftsx',
                    'leus', 'coae', 'yrfi', 'grx', 'pa4827', 'frr', 'sbcd',
                    'wzx', 'phaf', 'bdha', 'metf', 'biob', 'ftsa', 'ftsq',
                    'upps', 'glcf', 'glce', 'glcd', 'glcc', 'typa', 'rplo',
                    'prpl', 'rplc', 'gyrb', 'pepa', 'kdpa', 'kdpb', 'kdpc',
                    'rnt', 'ppka', 'moea2', 'lola', 'ribe', 'pchf', 'pche',
                    'phhc', 'dnab', 'tata', 'tatb', 'tatc', 'exaa', 'biof',
                    'trmu', 'rlud', 'arsr', 'arsb', 'arsc', 'prlc', 'waap',
                    'waag', 'hflc', 'hflk', 'pa4133', 'aotj', 'aotq', 'aotm',
                    'aotp', 'argr', 'phne', 'kdpf', 'sure', 'nrda', 'hemn',
                    'wbpw', 'rmd']
        self.assertEqual(g.process_sample_file(gene_eco, eco_assoc[0])[0],
                         eco_list)
        self.assertEqual(g.process_sample_file(gene_pse, pse_assoc[0])[0],
                         pse_list)

        def testSampleAnnotation(self):
            eco_assoc = g.process_association_file('gene_association.ecocyc')
            pse_assoc = g.process_association_file('gene_association.pseudocap')
            self.assertTrue(['go:0016051', 'go:0009376', 'go:0042150'] in \
                            g.process_sample_file(gene_eco, eco_assoc[0])[4])
            self.assertTrue(['go:0016787', 'go:0006950', 'go:0009636'] in \
                            g.process_sample_file(gene_pse, pse_assoc[0])[4])
        
if __name__ == '__main__':
    unittest.main()
