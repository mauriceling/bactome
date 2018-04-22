'''
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
    '''
    Class to read a Genbank file, parse it using BioPython library, into a Python dictionary. If the Genbank file consists of more than one Genbank records, the respective Genbank record name will be used as key in the dictionary.

    A Genbank record consists of the following keys (as described in 
    https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html):
        - 'name' is the Genbank record name (usually accession number)
        - 'description'
        - 'accessions'
        - 'comment'
        - 'data_file_division'
        - 'date'
        - 'features' is a list of <feature dictionary>
        - 'keywords'
        - 'organism'
        - 'sequence'
        - 'sequence_version'
        - 'source'
        - 'taxonomy'
        - 'topology'

    where <feature dictionary> consists of the following keys:
        - 'id'
        - 'start'
        - 'end'
        - 'type'
        - 'strand'
    '''

    def __init__(self):
        pass

    def _read_gb_features(self, gb_record):
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
        return list(self.records.keys())

    def getRecord(self, ID):
        return self.records[ID]

    def _getItem(self, ID, key):
        if key in self.records[ID]:
            return self.records[ID][key]
        else:
            return ''

    def getName(self, ID):
        return self._getItem(ID, 'name')

    def getDescription(self, ID):
        return self._getItem(ID, 'description')

    def getSequence(self, ID):
        return self._getItem(ID, 'sequence')

    def getOrganism(self, ID):
        return self._getItem(ID, 'organism')

    def getTaxonomy(self, ID):
        return self._getItem(ID, 'taxonomy')

    def getAccessions(self, ID):
        return self._getItem(ID, 'accessions')

    def getComment(self, ID):
        return self._getItem(ID, 'comment')

    def getDate(self, ID):
        return self._getItem(ID, 'date')

    def getKeywords(self, ID):
        return self._getItem(ID, 'keywords')

    def getSequenceVersion(self, ID):
        return self._getItem(ID, 'sequence_version')

    def getDataFileDivision(self, ID):
        return self._getItem(ID, 'data_file_division')

    def getSource(self, ID):
        return self._getItem(ID, 'source')

    def getTopology(self, ID):
        return self._getItem(ID, 'topology')

    def getFeatures(self, ID):
        return self._getItem(ID, 'features')


def recordSelector(gbfile='DSM6083.gb', RecordIndex=0, RecordID=None,
                   returntype='sequence'):
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
    '''
    python genbank.py readgb --gbfile=--gbfile=test/DSM6083.gb
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
    '''
    python genbank.py displayfeaturetypes --gbfile=test/DSM6083.gb
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
            print('%s: %i' % (k, count[k]))
        print()

def featureMap(feature='CDS', RecordIndex=0, RecordID=None,
               gbfile='DSM6083.gb', output='feature_result.txt'):
    '''
    python genbank.py featuring --feature=rRNA --RecordIndex=0 --gbfile=test/DSM6083.gb --output=feature_result.txt
    '''
    features = recordSelector(gbfile, RecordIndex, RecordID,
                              'features')
    required_feature = str(feature)
    data = []
    for f in features:
        if f['type'] == required_feature:
            data.append((f['start'], f['end']))
    f = open(output, 'w')
    for region in data:
        for base in range(int(region[0]), (region[1])+1):
            f.write(str(base) + '\n')
    f.close()


if __name__ == '__main__':
    exposed_functions = {'readgb': readGenBankFile,
                         'displayfeaturetypes': displayFeatureTypes,
                         'featuring': featureMap}
    fire.Fire(exposed_functions)
