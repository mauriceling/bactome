'''
Helper utilities / functions for Bactome

Date created: 22nd April 2018

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

def sequenceSelector(sequence, start=0, end=-1):
    if int(start) == 0 and int(end) == -1:
        pass
    elif int(start) > 0 and int(end) == -1:
        sequence = sequence[int(start):]
    elif int(start) == 0 and int(end) > -1:
        sequence = sequence[:int(end)]
    elif int(start) > 0 and int(end) > -1:
        sequence = sequence[int(start):int(end)]
    return sequence

def outputWriter(output, header, data):
    f = open(output, 'w')
    f.write(str(header) + '\n')
    for line in data:
        f.write(str(line) + '\n')
    f.close()