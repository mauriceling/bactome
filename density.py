'''!
Density Generator

Date created: 15th April 2019

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
import math

import fire

def cumulative_density(inputfile, column=1, increment=1, header=1):
    '''!
    Function to generate cumulative density (commonly known as 
    cumulative frequency) from a data set.

    Usage:

        python density.py CDF --inputfile=<input file path> --column=<column to use> --increment=<incremental value> --header=True

    @param inputfile String: Path of file to process.
    @param column Integer: Positional value of column to use. Default 
    = 1 (first column).
    @param increment Float: Incremental value for bin generation. 
    Default = 1.
    @param header Integer: Denotes the number of header rows to be 
    removed. Default = 1
    '''
    column = int(column) - 1
    data = open(inputfile, 'r').readlines()
    data = data[int(header):]
    data = [x.split(',')[column] for x in data]
    data = [float(x) for x in data]
    threshold = float(math.floor(min(data)))
    max_value = math.floor(max(data))
    num_of_data = len(data)
    print(' : '.join(['Threshold', 'Density']))
    while threshold <= max_value:
        temp = [x for x in data if x <= threshold]
        density = len(temp) / num_of_data
        print('%s : %s' % (str(threshold), str(density)))
        threshold = threshold + float(increment)

def probability_density(inputfile, column=1, increment=1, header=1):
    '''!
    Function to generate probability density (commonly known as 
    frequency) from a data set.

    Usage:

        python density.py PDF --inputfile=<input file path> --column=<column to use> --increment=<incremental value> --header=True

    @param inputfile String: Path of file to process.
    @param column Integer: Positional value of column to use. Default 
    = 1 (first column).
    @param increment Float: Incremental value for bin generation. 
    Default = 1.
    @param header Integer: Denotes the number of header rows to be 
    removed. Default = 1
    '''
    column = int(column) - 1
    data = open(inputfile, 'r').readlines()
    data = data[int(header):]
    data = [x.split(',')[column] for x in data]
    data = [float(x) for x in data]
    threshold = float(math.floor(min(data)))
    max_value = math.floor(max(data))
    num_of_data = len(data)
    print(' : '.join(['Threshold', 'Density']))
    while threshold <= max_value:
        temp = [x for x in data 
                if x <= threshold and \
                x > threshold - float(increment)]
        density = len(temp) / num_of_data
        print('%s : %s' % (str(threshold), str(density)))
        threshold = threshold + float(increment)


if __name__ == '__main__':
    exposed_functions = {'CDF': cumulative_density,
                         'PDF': probability_density}
    fire.Fire(exposed_functions)
