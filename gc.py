'''
Percent GC Calculator

Date created: 20th April 2018

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

import fire

from bactome_utils import sequenceSelector
from genbank import recordSelector


def sliderGC(gbfile='DSM6083.gb', RecordIndex=0, RecordID=None,
             start=0, end=-1, window=10000, interval=1):
    '''
    python gc.py slidegc --gbfile=test/DSM6083.gb --RecordIndex=0 --start=0 --end=-1 --window=10000 --interval=10 
    '''
    sequence = recordSelector(gbfile, RecordIndex, RecordID, 
                              'sequence')
    sequence = sequenceSelector(sequence, start, end)
    window = int(window)
    interval = int(interval)
    pointer = 1
    pdata = []
    while (pointer - 1 + window) < len(sequence):
        s = sequence[pointer-1:pointer-1+window]
        gc = float(s.count('G') + s.count('C')) * 100/window
        print('%i : %i : %.4f' % (pointer+start, 
                                  pointer-1+start+window,
                                  gc))
        pointer = pointer + interval

def blockGC(gbfile='DSM6083.gb', RecordIndex=0, RecordID=None,
            start=0, end=-1, window=10000):
    '''
    python gc.py blockgc --gbfile=test/DSM6083.gb --RecordIndex=0 --start=0 --end=-1 --window=10000
    '''
    sequence = recordSelector(gbfile, RecordIndex, RecordID,
                              'sequence')
    sequence = sequenceSelector(sequence, start, end)
    window = int(window)
    pointer = 1
    pdata = []
    while (pointer - 1 + window) < len(sequence):
        s = sequence[pointer-1:pointer-1+window]
        gc = float(s.count('G') + s.count('C')) * 100/window
        print('%i : %i : %.4f' % (pointer+start, 
                                  pointer-1+start+window,
                                  gc))
        pointer = pointer + window

def randomize(gbfile='DSM6083.gb', RecordIndex=0, RecordID=None,
              start=0, end=-1, window=10000, n=5000):
    '''
    python gc.py random --gbfile=test/DSM6083.gb --RecordIndex=0 --start=0 --end=-1 --window=10000 --n=200
    '''
    sequence = recordSelector(gbfile, RecordIndex, RecordID,
                              'sequence')
    sequence = sequenceSelector(sequence, start, end)
    choices = list(range(len(sequence)))
    window = int(window)
    count = 0
    data = []
    pdata = []
    while count < int(n):
        index = [random.choice(choices) for i in range(window)]
        s = [sequence[i] for i in index]
        s = ''.join(s)
        gc = float(s.count('G') + s.count('C')) * 100/window
        data.append(gc)
        print('%i : %.4f' % (count+1, gc))
        count = count + 1
    print()
    mean_gc = sum(data) / len(data)
    stdev_gc = (sum([(i-mean_gc)**2 for i in data]) / \
        (len(data) - 1)) ** 0.5
    print('Average Percent GC: %.4f' % mean_gc)
    print('Standard Deviation Percent GC: %.5f' % stdev_gc)


if __name__ == '__main__':
    exposed_functions = {'slidegc': sliderGC,
                         'blockgc': blockGC,
                         'random': randomize}
    fire.Fire(exposed_functions)
