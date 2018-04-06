'''
Calculates generation time of suspension cells.

Copyright 2009 Maurice Ling <mauriceling@acm.org> on behalf of
Maurice Ling, Chin-How Lee, Jack Si-Hao Oon, Kun-Cheng Lee.

@see: Lee, CH, Lee, KC, Oon, JSH, Ling, MHT. 2010. Bactome, I: 
Python in DNA Fingerprinting. In: Peer-Reviewed Articles from 
PyCon Asia-Pacific 2010. The Python Papers 5(3): 6.

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
'''

def generation(A2, A1, min2, min1,
               conv=50000000/0.3,
               limit=0.301,
               gradient=52137400,
               intercept=118718650):
    """
    Calculates the number of generations and generation time using
    turbidity measurement of cell culture (usually OD600 readings)
    between 2 time points (min1 and min2).
    
    Parameters:
        A2 = turbidity measurement at min2.
        A1 = turbidity measurement at min1.
        conv = conversion ratio between OD readings to cell count.
            Default is 0.3 OD = 5e7 (50 million) cells, standard for
            Escherichia coli.
        limit = OD reading at which limit of proportionality of the
            conversion ratio no longer holds, such as cell size decrease
            or increase with higher OD readings.
            Default is 0.301, for Escherichia coli.
        gradient and intercept = the gradient and intercept for cell
            number correction when OD is above limit of proportionality.
            Correction equation:
                cell count = gradient * ln(OD) + intercept
            Default value of gradient is 52137400 and intercept is
            118718650.

        Limit of proportionality, correction gradient and intercept
        are reverse engineered from the graph of Sezonov, Guennadi,
        Joseleau-Petit, Daniele and D'Ari, Richard. 2007. Escherichia
        coli Physiology in Luria-Bertani Broth. Journal of Bacteriology
        189(23):8746-8749
    
    Returns:
        (number of generations between min1 and min2,
        length of time (in minutes) for each generation (generation time),
        number of cells per ml at min2,
        number of cells per ml at min1)
    """
    if A1 < limit:
        y = A1 * conv
    else:
        y = gradient * math.log(A1, math.e) + intercept
    if A2 < limit:
        x = A2 * conv
    else:
        x = gradient * math.log(A2, math.e) + intercept
    gen = (math.log10(x) - math.log10(y))/math.log10(2)
    time_diff = min2 - min1 
    return (gen, time_diff/gen, x, y)

def totalgen(ODs, p='yes', dilution=100,
             conv=50000000/0.3, limit=0.301,
             gradient=52137400, intercept=118718650):
    """
    Parameters:
        ODs = list of OD readings at the start of next subculture.
        p = flag to determine if the data is to be printed. Default = yes
        dilution = dilution factor between each passage. Default = 100
        conv = conversion ratio between OD readings to cell count.
            Default is 0.3 OD = 5e7 (50 million) cells, standard for
            Escherichia coli.
        limit = OD reading at which limit of proportionality of the
            conversion ratio no longer holds, such as cell size decrease
            or increase with higher OD readings.
            Default is 0.301, for Escherichia coli.
        gradient and intercept = the gradient and intercept for cell
            number correction when OD is above limit of proportionality.
            Correction equation:
                cell count = gradient * ln(OD) + intercept
            Default value of gradient is 52137400 and intercept is
            118718650.

        Limit of proportionality, correction gradient and intercept
        are reverse engineered from the graph of Sezonov, Guennadi,
        Joseleau-Petit, Daniele and D'Ari, Richard. 2007. Escherichia
        coli Physiology in Luria-Bertani Broth. Journal of Bacteriology
        189(23):8746-8749
    """
    ODs = [float(x) for x in ODs if x != '']
    gen_list = [0] * (len(ODs) - 1)
    for passage in range(len(ODs)-1):
        if ODs[passage] < limit:
            initial_count = ODs[passage] * conv
        else:
            initial_count = gradient * math.log(ODs[passage], math.e) + \
                            intercept
        if ODs[passage+1] < limit:
            final_count = ODs[passage+1] * limit
        else:
            final_count = gradient * math.log(ODs[passage+1], math.e) + \
                          intercept
        initial_count = initial_count / float(dilution)
        gen = (math.log10(final_count) - \
               math.log10(initial_count)) / math.log10(2)
        if gen < 0: gen = 0.0
        gen_list[passage] = gen
        if p == 'yes':
            print 'Number of cells inoculate from Passage ' + \
                  str(passage) + ' = ' +str(initial_count)
            print 'Number of cells at subculture from Passage ' + \
                  str(passage+1) + ' = ' +str(final_count)
            print 'Between Passage ' + str(passage) + ' and Passage ' + \
                  str(passage+1) + ', Number of generations = ' + \
                  str(gen_list[passage])
            print 'Total number of generations since Passage 0 = ' + \
                  str(sum(gen_list))
            print
    return gen_list