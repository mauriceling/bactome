"""
Methods for identifying invariant genes (suitable reference/
normalization genes) from expression data 

Date created: 16th April 2013

License: GNU General Public License version 3 for academic or 
not-for-profit use only
"""

import sys
import time
import math
import sets

from invariant_gene import selfed_correlation
from invariant_gene import selfed_ratio_correlation
from invariant_gene import selfed_product_correlation
from invariant_gene import regression_ratio
from invariant_gene import average_stdev
from invariant_gene import cv
from invariant_gene import gradient

def datafile(filename):
    print 'Reading data from', filename
    t = time.time()
    f = open(filename, 'r').readlines()
    if len(f) == 1:
        f = f[0].split('\r')
    f = [x[:-1] for x in f]
    f = [x.split(',') for x in f]
    data = {}
    for x in f:
        try:
            x = [y for y in x if y != '']
            if len(x) > 2: 
                name = str(x[0])
                values = [float(y) for y in x[1:]]
                data[name] = values
        except: 'Error formatting:', x
    print '     Time taken:', time.time()-t, 'seconds'
    return data

def p_gradient(data, options, methods, results):
    print 'Calculating gradients ...'
    t = time.time()
    t_results = gradient(data)
    methods.append('gradient')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['gradient'] = t_results[genename]
        else:
            results[genename] = {'gradient': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)

def p_cv(data, options, methods, results):
    print 'Calculating coefficient of variations ...'
    t = time.time()
    t_results = cv(data)
    methods.append('cv')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['cv'] = t_results[genename]
        else:
            results[genename] = {'cv': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)

def p_regression_ratio(data, options, methods, results):
    print 'Calculating regression ratios (R^2/slope) ...'
    t = time.time()
    t_results = regression_ratio(data)
    methods.append('R^2/slope')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['R^2/slope'] = t_results[genename]
        else:
            results[genename] = {'R^2/slope': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)
    
def p_average_stdev(data, options, methods, results):
    print 'Calculating average x stdev ...'
    t = time.time()
    t_results = average_stdev(data)
    methods.append('avg*SD')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['avg*SD'] = t_results[genename]
        else:
            results[genename] = {'avg*SD': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)
    
def p_selfed_product_correlation(data, options, methods, results):
    print 'Calculating selfed product correlation ...'
    t = time.time()
    if not options.has_key('rss'): options['rss'] = 30
    t_results = selfed_product_correlation(data, options['rss'])
    methods.append('prod_corr')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['prod_corr'] = t_results[genename]
        else:
            results[genename] = {'prod_corr': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)
    
def p_selfed_ratio_correlation(data, options, methods, results):
    print 'Calculating selfed ratio correlation ...'
    t = time.time()
    if not options.has_key('rss'): options['rss'] = 30
    t_results = selfed_ratio_correlation(data, options['rss'])
    methods.append('ratio_corr')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['ratio_corr'] = t_results[genename]
        else:
            results[genename] = {'ratio_corr': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)
    
def p_selfed_correlation(data, options, methods, results):
    print 'Calculating selfed correlation ...'
    t = time.time()
    if not options.has_key('rss'): options['rss'] = 30
    t_results = selfed_correlation(data, options['rss'])
    methods.append('self_corr')
    for genename in t_results:
        if results.has_key(genename):
            results[genename]['self_corr'] = t_results[genename]
        else:
            results[genename] = {'self_corr': t_results[genename]}
    print '     Time taken:', time.time()-t, 'seconds'
    return (methods, results)

def p_geomean_expratio_cv(data, options, methods, results):
    print 'Calculating geometric mean of ratio correlation and CV ...'
    (methods, results) = p_selfed_ratio_correlation(data, options, 
                                                    methods, results)
    (methods, results) = p_cv(data, options, methods, results)
    methods.append('geomean_expratio_cv')
    for genename in results.keys():
        eratio = math.e ** results[genename]['ratio_corr']
        cv = results[genename]['cv']
        geomean = math.sqrt(eratio * cv)
        if results.has_key(genename): 
            results[genename]['geomean_expratio_cv'] = geomean
        else:
            results[genename] = {'geomean_expratio_cv': geomean}
    return (methods, results)

def p_avgexpratio_avgcv(data, options, methods, results):
    print 'Calculating sum of averaged ratio correlation and averaged CV ...'
    (methods, results) = p_selfed_ratio_correlation(data, options, 
                                                    methods, results)
    (methods, results) = p_cv(data, options, methods, results)
    methods.append('avgexpratio_avgcv')
    avgexpratio = [math.e ** results[genename]['ratio_corr'] 
                   for genename in results.keys()]
    avgexpratio = sum(avgexpratio) / float(len(avgexpratio))
    avgcv = [results[genename]['cv'] 
             for genename in results.keys()]
    avgcv = sum(avgcv) / float(len(avgcv))
    for genename in results.keys():
        eratio = math.e ** results[genename]['ratio_corr']
        cv = results[genename]['cv']
        r = (float(eratio) / avgexpratio) + (float(cv) / avgcv)
        if results.has_key(genename): 
            results[genename]['avgexpratio_avgcv'] = r
        else:
            results[genename] = {'avgexpratio_avgcv': r}
    return (methods, results)
    
def p_xxx(data, options, methods, results):
    return (methods, results)
    
def processor(data, options):
    results = {}
    methods = []
    if options.has_key('all_methods') or options.has_key('gradient'):
        # Method: gradient
        (methods, results) = p_gradient(data, options, methods, results)
    if options.has_key('all_methods') or options.has_key('cv'):
        # Method: CV
        (methods, results) = p_cv(data, options, methods, results)
    if options.has_key('all_methods') or options.has_key('regression'):
        # Method: regression ratios (R^2/slope)
        (methods, results) = p_regression_ratio(data, options, 
                                                methods, results)
    if options.has_key('all_methods') or options.has_key('avgstd'):
        # Method: average x stdev
        (methods, results) = p_average_stdev(data, options, 
                                             methods, results)
    if options.has_key('all_methods') or options.has_key('pcorr'):
        # Method: selfed product correlation
        (methods, results) = p_selfed_product_correlation(data, options, 
                                                          methods, results)
    if options.has_key('all_methods') or options.has_key('rcorr'):
        # Method: selfed ratio correlation
        (methods, results) = p_selfed_ratio_correlation(data, options, 
                                                        methods, results)
    if options.has_key('all_methods') or options.has_key('scorr'):
        # Method: selfed correlation
        (methods, results) = p_selfed_correlation(data, options, 
                                                  methods, results)
    # Method: geomean(e^ratio, CV)
    (methods, results) = p_geomean_expratio_cv(data, options, 
                                               methods, results)
    # Method: (e^ratio/avg(e^ratio)) + (cv/avg(cv))
    (methods, results) = p_avgexpratio_avgcv(data, options, 
                                             methods, results)
    return (methods, results)

def results_writer(filename, methods, results):
    out = open(filename, 'w')
    out.write('ResultFile,' + ','.join(methods) + '\n')
    for genename in results:
        values = ','.join([str(results[genename][m])
                           for m in methods])
        out.write(','.join([genename, values]) + '\n')
    out.close()

def option_processor(argv):
    options = {}
    options['application'] = argv[0]
    argv = [x[1:] for x in argv[1:] if x.startswith('-')]
    for arg in argv:
        arg = arg.split(':')
        if len(arg) == 2:
            options[arg[0]] = arg[1]
        else:
            options[arg[0]] = ''
    return options

def print_usage():
    print '''
    OLIgonucleotide Variable Expression Ranker (OLIVER) 1.0:
    A tool for identifying suitable reference (invariant) genes from 
    large transcriptomme datasets.
    
    Date created: 16th April 2013
    License: GNU General Public License version 3 for academic or not-
    for-profit use only

    Usage: python oliver.py <expression filename> <results filename> [options]
    where
    
    <expression filename> is a comma-delimited file in the
    format of
    ProbeName1,Sample1Value,Sample2value,Sample3Value, ...
    ProbeName2,Sample1Value,Sample2value,Sample3Value, ...
    ProbeName3,Sample1Value,Sample2value,Sample3Value, ...
    ...

    [options] in the format of -<option key>:<option value>
    For example, -outfmt:5 (option with value) and -fmt
    (option without value)

    Allowed options (please refer to user manual for detailed description):
    
    -all_methods    Calculate for all available methods of reference genes 
                    identification
    -avgstd         Calculates average x standard deviation
    -cv             Calculates coefficient of variation
    -gradient       Calculates gradient
    -help           Prints / displays usage message
    -regression     Calculates regression ratio (R^2/slope)
    -pcorr          Calculates selfed product correlation
    -rcorr          Calculates selfed ratio correlation
    -rss            Define the number of random samples for pair-wise 
                    correlations. If not given, the default is 30.
    -scorr          Calculates selfed correlation
    '''
    
if __name__=='__main__':
    if len(sys.argv) < 3:
        print print_usage()
    else:
        print '''
        OLIgonucleotide Variable Expression Ranker (OLIVER) 1.0:
        A tool for identifying suitable reference (invariant) genes from 
        large transcriptomme datasets.
    
        Date created: 16th April 2013
        License: GNU General Public License version 3 for academic or 
        not-for-profit use only
        '''
        options = option_processor(sys.argv)
        if options.has_key('help'):
            print print_usage()
            sys.exit(0)
        data = datafile(sys.argv[1])
        (methods, results) = processor(data, options)
        methods = list(sets.Set(methods))
        results_writer(sys.argv[2], methods, results)
        print
        print 'Summary ......'
        print 'Input data file:', sys.argv[1]
        print 'Number of genes/probes in input data file:', len(data)
        print 'Output results file:', sys.argv[2]
        print 'Number of methods used:', len(methods)
        print 'Methods:', ', '.join(methods)
        print '===== end of summary ====='
        
