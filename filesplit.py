'''!
File Splitter

Date created: 8th May 2018

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
from itertools import islice

import fire

def splitByLine(inputfile, size=5000):
    '''!
    python filesplit.py line --inputfile=test.txt --size=5000
    '''
    size = int(size)
    count = 1
    with open(inputfile, 'r') as f_in:
        while True:
            block = list(islice(f_in, size))
            if not block:
                break
            outputfile = inputfile + '.' + str(count) + '.txt'
            f_out = open(outputfile, 'w')
            for line in block:
                f_out.write(line)
            f_out.close()
            count = count + 1


if __name__ == '__main__':
    exposed_functions = {'line': splitByLine}
    fire.Fire(exposed_functions)
