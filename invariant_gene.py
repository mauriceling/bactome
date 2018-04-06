"""
Methods for identifying invariant genes (suitable reference/
normalization genes) from expression data 

@see: Chan, OYW, Keng, BMH, Ling, MHT. 2014. Correlation and Variation 
Based Method for Reference Genes Identification from Large Datasets. 
Electronic Physician 6(1): 719-727.

Date created: 29th March 2012

License: GNU General Public License version 3

Bactome package is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import random    
from copads.samplestatistics import SingleSample, TwoSample

def selfed_correlation(data, randomsize):
    """
    The average absolute pairwise correlation between dataset X 
    and sample non-X data (n = randomsize) where small value 
    represents higher stability
    
    Pseudocode:
    1. tested_gene <-- list of expression for gene to be tested
    2. dataset_X <-- list of expression for any other gene
    3. q_set <-- {dataset_X(i) | 1 < i < n}
    4. correlation <-- {|r(dataset_X, q_set)|}, repeat steps 2 to 4
    for all genes that are not tested_gene
    5. average_correlation(tested_gene) <-- average of correlation
    6. repeat steps 1 to 5 to test for all genes
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @param randomsize: number of random samples to take
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    test_list = data.keys()
    if len(test_list) > randomsize:
        test_list = random.sample(test_list, randomsize)
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        pcc = []
        for tester in test_list:
            tester = [float(x) for x in data[tester]]
            sample = TwoSample(gene, 'genename', tester, 'tester')
            try: pcc.append(abs(sample.pearson()))
            except: pass
        results[genename] = float(sum(pcc)) / len(pcc)
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results
    
def selfed_ratio_correlation(data, randomsize):
    """
    The average absolute pairwise correlation between dataset X 
    and the quotient of X and sample of non-X data (n = randomsize) 
    where small value represents higher stability
    
    Pseudocode:
    1. tested_gene <-- list of expression for gene to be tested
    2. dataset_X <-- list of expression for any other gene
    3. q_set <-- {dataset_X(i) / tested_gene(i) | 1 < i < n}
    4. correlation <-- {|r(dataset_X, q_set)|}, repeat steps 2 to 4
    for all genes that are not tested_gene
    5. average_correlation(tested_gene) <-- average of correlation
    6. repeat steps 1 to 5 to test for all genes
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @param randomsize: number of random samples to take
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    test_list = data.keys()
    if len(test_list) > randomsize:
        test_list = random.sample(test_list, randomsize)
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        pcc = []
        for tester in test_list:
            tester = [float(x) for x in data[tester]]
            tester = [tester[i] / gene[i]
                      for i in range(len(gene))]
            sample = TwoSample(gene, 'genename', tester, 'tester')
            try: pcc.append(abs(sample.pearson()))
            except: pass
        results[genename] = float(sum(pcc)) / len(pcc)
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results

def selfed_product_correlation(data, randomsize):
    """
    The average absolute pairwuse correlation between dataset X 
    and the product of X and sample of non-X data (n = randomsize) 
    where large value represents higher stability
    
    Pseudocode:
    1. tested_gene <-- list of expression for gene to be tested
    2. dataset_X <-- list of expression for any other gene
    3. q_set <-- {dataset_X(i) * tested_gene(i) | 1 < i < n}
    4. correlation <-- {|r(dataset_X, q_set)|}, repeat steps 2 to 4
    for all genes that are not tested_gene
    5. average_correlation(tested_gene) <-- average of correlation
    6. repeat steps 1 to 5 to test for all genes
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @param randomsize: number of random samples to take
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    test_list = data.keys()
    if len(test_list) > randomsize:
        test_list = random.sample(test_list, randomsize)
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        pcc = []
        for tester in test_list:
            tester = [float(x) for x in data[tester]]
            tester = [tester[i] * gene[i]
                      for i in range(len(gene))]
            sample = TwoSample(gene, 'genename', tester, 'tester')
            try: pcc.append(sample.pearson())
            except: pass
        results[genename] = float(sum(pcc)) / len(pcc)
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results

def regression_ratio(data):
    """
    The quotient of coefficient of determination and gradient where
    large value represents higher stability.

    @see Lee et al. 2007. Identification of novel universal
    housekeeping genes by statistical analysis of microarray data.
    Journal of Biochemistry and Molecular Biology 40(2):226-231.
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        tester = range(len(gene))
        sample = TwoSample(gene, 'genename', tester, 'nominal')
        (gradient, intercept) = sample.linear_regression()
        if gradient == 0.0: gradient = 0.001
        pcc = sample.pearson()
        results[genename] = (pcc * pcc) / gradient
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results

def average_stdev(data):
    """
    The product of arithmetic mean and standard deviation where
    ?small? value represents higher stability.

    @see Lee et al. 2007. Identification of novel universal
    housekeeping genes by statistical analysis of microarray data.
    Journal of Biochemistry and Molecular Biology 40(2):226-231.
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        gene = SingleSample(gene)
        results[genename] = (gene.variance() ** 0.5) * \
                            gene.arithmeticMean()
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results

def cv(data):
    """
    The coefficient of variation where small value represents
    higher stability.
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        gene = SingleSample(gene)
        results[genename] = (gene.variance() ** 0.5) / \
                            gene.arithmeticMean()
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results

def gradient(data):
    """
    The linear regression gradient where small value represents
    higher stability.
    
    @param data: input data as a dictionary (processed by datafile function 
    in oliver.py) where key is ProbeName and value is a list of 
    SampleValues
    @return: dictionary where key is ProbeName and value is result from 
    current reference gene identification method
    """
    results = {}
    count = 0
    for genename in data.keys():
        gene = [float(x) for x in data[genename]]
        tester = range(len(gene))
        sample = TwoSample(gene, 'genename', tester, 'nominal')
        (gradient, intercept) = sample.linear_regression()
        results[genename] = abs(gradient)
        count = count + 1
        if count % 500 == 0: print count, 'gene/probe processed'
    return results
