"""
Functions for computer-intensive methods in statistics

Copyright 2009 Maurice Ling <mauriceling@acm.org>

Licence: GNU General Public License version 3

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
    
def randomization(statfunction, fulldata, testdata, iteration=1000):
    """
    Test harness for statistical randomization test, also known as
    permutation test or bootstrapping.

    Pitman, EJG. 1938. Significance tests which may be applied to samples
    from any population. Part III. The analysis of variance test.
    Biometrika 29: 322-335. 
    Manly, Bryan FJ. 2007. Randomization, Bootstrap and Monte Carlo
    Methods in Biology (3rd edition). Chapman & Hall/CRC.

    Parameters:
        statfunction = Python function to calculate the test statistic.
        fulldata = list of entire data set (for randomization)
        testdata = list of sample data of interest (for testing)
        iteration = number of times random samples are drawn.
            Default = 1000

    Returns:
        Proportion where statistic from randomized sample is greater
        or equal to the statistic from testdata. This can be taken as
        p-value for the test
    """
    from copads.Operations import sample_wr
    test_statistic = statfunction(testdata)
    randomized_statistic = [statfunction(sample_wr(fulldata, len(testdata))) 
                            for i in range(iteration)]
    #print randomized_statistic, test_statistic
    more = [1 for x in randomized_statistic
                if x >= test_statistic]
    return sum(more) / float(len(randomized_statistic))

def sample_bootstrap(data, iteration=1000):
    """
    Calculates bootstrapped means and standard deviation.
    
    Efron, Bradley. 1979. Bootstrap methods: Another look at the jackknife.
    The Annals of Statistics 7:1-26. 
    Manly, Bryan FJ. 2007. Randomization, Bootstrap and Monte Carlo
    Methods in Biology (3rd edition). Chapman & Hall/CRC.

    Parameters:
        data = list of data to generate bootstrapped samples from.
        iterations = number of times bootstrapped samples are drawn.
            Default = 1000

    Returns:
        (bootstrapped means, bootstrapped standard deviation)
    """
    from copads.Operations import sample_wr
    len_data = len(data)
    pool = [sum(sample_wr(data, len_data)) / float(len_data)
            for x in range(iteration)]
    poolmeans = float(sum(pool)) / float(len(pool))
    div = [(x - poolmeans)**2 for x in pool]
    stdev = math.sqrt(sum(div)) / (len(pool) - 1)
    return (poolmeans, stdev)
	
def sample_jackknife(data):
	"""
	Estimates mean and standard deviation from a single sample
	using Jackknife method.
	
	Wu, CFJ. 1986. Jackknife, Bootstrap and other resampling methods 
	in regression analysis. The Annals of Statistics 14:1261–1295.
	
	Parameters:
        data = list of data to generate bootstrapped samples from.
		
	Returns:
        (estimated means, estimated standard deviation)
	"""
	len_data = len(data)
	pool = [sum(data.pop(i)) / float(len_data -1) 
			for i in range(len_data)]
	poolmeans = float(sum(pool)) / float(len(pool))
    div = [(x - poolmeans)**2 for x in pool]
    stdev = math.sqrt(sum(div)) / (len(pool) - 1)
    return (poolmeans, stdev)

def electrophoresis_mw(marker, samples):
    """
    Parameters:
        marker = data of the marker lane in the format of
            [list of migration distance, list of molecular weights]
        samples = 
    """
    from copads.SampleStatistics import TwoSample
    if len(marker[0]) != len(marker[1]):
        raise AttributeError, 'Number of molecular weights not equal to the \
        number of migration distances'
    marker[1] = [math.log10(x) for x in marker[1]]
    (gradient, intercept) = TwoSample(marker[0], 'dist', 
                                      marker[1], 'mw').linear_regression()
    return [[10 ** (gradient * band) + 10 ** intercept 
            for band in sample] 
                for sample in samples]
    

