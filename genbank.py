'''!
Genbank file reader and processor

Date created: 7th May 2017

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
import fire
from Bio import SeqIO


class GenBankFile(object):
    '''!
    Class to read a Genbank file, parse it using BioPython library, 
    into a Python dictionary. If the Genbank file consists of more 
    than one Genbank records, the respective Genbank record name will 
    be used as key in the dictionary.

    A Genbank record consists of the following keys (as described in 
    https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html):
        - 'name' is the Genbank record name (usually accession number)
        - 'description'
        - 'accessions'
        - 'comment'
        - 'data_file_division'
        - 'date'
        - 'features' is a list of feature dictionary
        - 'keywords'
        - 'organism'
        - 'sequence'
        - 'sequence_version'
        - 'source'
        - 'taxonomy'
        - 'topology'

    where feature dictionary consists of the following keys:
        - 'id'
        - 'start'
        - 'end'
        - 'type'
        - 'strand'
    '''

    def __init__(self):
        '''!
        Constructor method.
        '''
        pass

    def _read_gb_features(self, gb_record):
        '''!
        Private method - used to read the features of a Genbank 
        record.
        '''
        features = [0] * len(gb_record.features)
        for i in range(len(gb_record.features)):
            f = gb_record.features[i]
            record = {'type': f.type,
                      'start': f.location.start.position,
                      'end': f.location.end.position}
            try: record['id'] = f.id
            except: pass
            try: record['strand'] = f.strand
            except: pass
            for key in f.qualifiers.keys():
                record[key] = f.qualifiers[key]
            features[i] = record
        return features

    def _read_gb_annotations(self, record, gb_record):
        '''!
        Private method - used to read the annotations of a Genbank 
        record.
        '''
        annotations = gb_record.annotations
        try:
            record['accessions'] = annotations['accessions']
        except KeyError: pass
        try:
            comment = annotations['comment'].split('\n')
            record['comment'] = '|'.join([x.strip() for x in comment])
        except KeyError: pass
        try:
            record['data_file_division'] = annotations['data_file_division']
        except KeyError: pass
        try:
            record['date'] = annotations['date']
        except KeyError: pass
        try:
            record['keywords'] = annotations['keywords']
        except KeyError: pass
        try:
            record['organism'] = annotations['organism']
        except KeyError: pass
        try:
            record['sequence_version'] = annotations['sequence_version']
        except KeyError: pass
        try:
            record['source'] = annotations['source']
        except KeyError: pass
        try:
            record['taxonomy'] = annotations['taxonomy']
        except KeyError: pass
        try:
            record['topology'] = annotations['topology']
        except KeyError: pass
        return record
    
    def readGB(self, filepath):
        '''!
        Method to read a GenBank file and parse it into a dictionary. 

        @param filepath string: file path of the Genbank file
        '''
        self.filepath = filepath
        gb = open(self.filepath, 'r')
        self.records = {}
        for gb_record in SeqIO.parse(gb, 'genbank'):
            record = {'name': gb_record.name,
                      'description': gb_record.description}
            record['features'] = self._read_gb_features(gb_record)
            record['sequence'] = str(gb_record.seq)
            record = self._read_gb_annotations(record, gb_record)
            self.records[gb_record.name] = record

    def getIDs(self):
        '''!
        Method to return all the identifying IDs for the Genbank 
        records.

        @return A list of names (identifying IDs) for all the read 
        Genbank records
        '''
        return list(self.records.keys())

    def getRecord(self, ID):
        '''!
        Method to return the data for a specific Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record

        @return Dictionary representing the Genbank record 
        '''
        return self.records[ID]

    def _getItem(self, ID, key):
        '''!
        Private method - used to get the value for a field for a 
        specific Genbank record.
        '''
        if key in self.records[ID]:
            return self.records[ID][key]
        else:
            return ''

    def getName(self, ID):
        '''!
        Method to get the name of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'name')

    def getDescription(self, ID):
        '''!
        Method to get the description of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'description')

    def getSequence(self, ID):
        '''!
        Method to get the sequence of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'sequence')

    def getOrganism(self, ID):
        '''!
        Method to get the organism name of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'organism')

    def getTaxonomy(self, ID):
        '''!
        Method to get the taxonomy of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'taxonomy')

    def getAccessions(self, ID):
        '''!
        Method to get the Accession number of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'accessions')

    def getComment(self, ID):
        '''!
        Method to get the comment of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'comment')

    def getDate(self, ID):
        '''!
        Method to get the date of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'date')

    def getKeywords(self, ID):
        '''!
        Method to get the keyword(s) of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'keywords')

    def getSequenceVersion(self, ID):
        '''!
        Method to get the sequence version of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'sequence_version')

    def getDataFileDivision(self, ID):
        '''!
        Method to get the data file division of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'data_file_division')

    def getSource(self, ID):
        '''!
        Method to get the source of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'source')

    def getTopology(self, ID):
        '''!
        Method to get the topology of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'topology')

    def getFeatures(self, ID):
        '''!
        Method to get the feature(s) of the Genbank record.

        @param ID String: Unique name (identifying IDs) for the 
        specific Genbank record
        '''
        return self._getItem(ID, 'features')


def recordSelector(gbfile, RecordIndex=0, RecordID=None,
                   returntype='sequence'):
    '''!
    Function to select a Genbank record in a multi-record Genbank
    file. 
    
    @param gbfile string: file path of the Genbank file
    @param RecordIndex Integer: Relative Genbank record number in 
    the given Genbank file where 0 is the first record and 1 is 
    the secord record, and so on
    @param RecordID String: Unique name (identifying IDs) for the 
    specific Genbank record
    @param returntype String: Defines the data to return; currently, 
    only 'sequence' and 'features' are acceptable
    '''
    gb = GenBankFile()
    gb.readGB(gbfile)
    IDs = list(gb.getIDs())
    if RecordIndex == 0 or RecordID == None:
        ID = IDs[0]
    elif RecordIndex != 0 and RecordID == None:
        ID = int(RecordIndex)
    elif RecordID != None:
        ID = RecordIndex
    if returntype == 'sequence':
        return gb.getSequence(ID)
    elif returntype == 'features':
        return gb.getFeatures(ID)

def readGenBankFile(gbfile='DSM6083.gb'):
    '''!
    Function to read a Genbank file and prints out descriptions 
    for each Genbank record in the file.

    Usage:
        
        python genbank.py readgb --gbfile=<Genbank file path>

    @param gbfile string: file path of the Genbank file
    '''
    gb = GenBankFile()
    gb.readGB(gbfile)
    print('Number of Records: ' + str(len(gb.getIDs())))
    print()
    for ID in list(gb.getIDs()):
        print('Name: ' + gb.getName(ID))
        print('Accession: ' + '|'.join(gb.getAccessions(ID)))
        print('Description: ' + gb.getDescription(ID))
        print('Source: ' + gb.getSource(ID))
        print('Topology: ' + gb.getTopology(ID))
        print('Organism: ' + gb.getOrganism(ID))
        print('Taxonomy: ' + '|'.join(gb.getTaxonomy(ID)))
        print('Data File Division: ' + gb.getDataFileDivision(ID))
        print('Comment: ' + gb.getComment(ID))
        print('Date: ' + gb.getDate(ID))
        print('Keywords: ' + '|'.join(gb.getKeywords(ID)))
        print('Sequence Version: ' + str(gb.getSequenceVersion(ID)))
        print('Length of Sequence: ' + str(len(gb.getSequence(ID))))
        features = gb.getFeatures(ID)
        print('Number of Features: ' + str(len(features)))

def displayFeatureTypes(gbfile='DSM6083.gb'):
    '''!
    Function to read a Genbank file and prints out the number of each 
    feature for each Genbank record in the file.

    Usage:
    
        python genbank.py displayfeaturetypes --gbfile=<Genbank 
        file path>

    The format of the output is

        <feature type> : <count>

    where count is the number of features for each feature type.

    @param gbfile string: file path of the Genbank file
    '''
    gb = GenBankFile()
    gb.readGB(gbfile)
    print('Number of Records: ' + str(len(gb.getIDs())))
    print()
    for ID in list(gb.getIDs()):
        features = gb.getFeatures(ID)
        count = {}
        for f in features:
            if f['type'] not in count:
                count[f['type']] = 1
            else:
                count[f['type']] = count[f['type']] + 1
        print('Name: ' + gb.getName(ID))
        for k in count:
            print('%s : %i' % (k, count[k]))
        print()

def featureMap(feature='CDS', RecordIndex=0, RecordID=None,
               gbfile='DSM6083.gb', expand=False, resolution=1):
    '''!
    Function to print out the base position(s) or base range for 
    a specific feature type in a specific Genbank record within 
    the given Genbank file.

    Usage 1: Prints the range of bases for a required feature
    
        python genbank.py featuring --feature=<feature type> 
        --RecordIndex=<Record number> --RecordID=<Genbank record name> 
        --gbfile=<Genbank file path> 

    Usage 2: Prints out individual base position for a required 
    feature
    
        python genbank.py featuring --feature=<feature type> 
        --RecordIndex=<record number> --RecordID=<Genbank record name> 
        --gbfile=<Genbank file path> 
        --expand=False --resolution=<resolution>

    where
        - record number is the relative Genbank record number in 
        the given Genbank file where 0 is the first record and 
        1 is the secord record, and so on
        - Genbank record name is the identifying ID for the Genbank 
        record
        - resolution is to control the resolution of printouts. 
        For example, resolution of 1 will print out every base 
        position for a given feature type; which will generate 
        a long printout. Increasing the resolution will reduce 
        the volume of printout but may not accurately pin-point 
        the end base position of the feature 

    In both usage, it requires either RecordIndex or RecordID.

    Usage 1 will generate the following output:

        <feature type> : <starting base position> : <ending 
        base position>

    Usage 2 will generate the following output:

        <feature type> : <base position>

    @param feature String: Type of feature
    @param RecordIndex Integer: Relative Genbank record number in 
    the given Genbank file where 0 is the first record and 1 is 
    the secord record, and so on
    @param RecordID String: Unique name (identifying IDs) for the 
    specific Genbank record
    @param gbfile string: file path of the Genbank file
    @param expand Boolean: Flag to define whether to expand from 
    start/end base positions to indivisual base positions
    @param resolution Integer: Control the resolution of printouts. 
    For example, resolution of 1 will print out every base 
    position for a given feature type; which will generate a long 
    printout. Increasing the resolution will reduce the volume of 
    printout but may not accurately pin-point the end base position 
    of the feature 
    '''
    features = recordSelector(gbfile, RecordIndex, RecordID,
                              'features')
    required_feature = str(feature)
    expand = str(expand).lower()
    for f in features:
        if f['type'] == required_feature:
            if expand == 'false':
                print('%s : %s : %s' % \
                      (feature, str(f['start']), str(f['end'])))
            if expand == 'true':
                resolution = int(resolution)
                for base in range(f['start'], f['end']+1,
                                  resolution):
                    print('%s : %s' % (feature, base))


if __name__ == '__main__':
    exposed_functions = {'readgb': readGenBankFile,
                         'displayfeaturetypes': displayFeatureTypes,
                         'featuring': featureMap}
    fire.Fire(exposed_functions)
