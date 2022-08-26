'''!
Brainopy

Date created: 18th August 2022

License: GNU General Public License version 3 for academic or 
not-for-profit use only.

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
import sqlite3
import uuid

class brainopy(object):
    def __init__(self, brainDB=None):
        if brainDB == None:
            self.con = None
            self.cur = None 
        else:
            self.connectBrain(brainDB)

    def connectBrain(self, brainDB):
        """!
        Connects to the brain database specified by the brainDB, which is a SQLite database.
        @param brainDB: Path to Brain database file
        @return: (Connection object or None, Cursor object or None)
        """
        self.con = sqlite3.connect(brainDB)
        self.cur = self.con.cursor()
        self.cur.execute("CREATE TABLE IF NOT EXISTS ID_table (ID text primary key, table_name text)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS neurotransmitter (neurotransmitter text primary key, description text)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS neuron_state (ID text, neurotransmitter text, value real)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS dendrite_state (ID text, neurotransmitter text, value real)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS axon_state (ID text, neurotransmitter text, value real)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS synapse_state (ID text, neurotransmitter text, value real)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS neuron_body (ID text, neuron_state_ID text, axon_state_ID text)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS neuron_dendrite (ID text, dendrite_state_ID text)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS axon_synapse_link (axon_state_ID text, synapse_state_ID text)")
        self.cur.execute("CREATE TABLE IF NOT EXISTS synapse_dendrite_link (synapse_state_ID text, dendrite_state_ID text)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS ID_table_index ON ID_table (ID, table_name)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS ID_table_table ON ID_table (table_name)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_state_ID ON neuron_state (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS dendrite_state_ID ON dendrite_state (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS axon_state_ID ON axon_state (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS synapse_state_ID ON synapse_state (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_body_ID ON neuron_body (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_dendrite_ID ON neuron_dendrite (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS axon_synapse_link_index ON axon_synapse_link (axon_state_ID, synapse_state_ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS synapse_dendrite_link_index ON synapse_dendrite_link (synapse_state_ID, dendrite_state_ID)")
        self.con.commit()

    def disconnectBrain(self):
        self.con.commit()
        self.con.close()

    def addNeurotransmitters(self, neurotransmitters):
        for key in neurotransmitters:
            try: 
                self.cur.execute("INSERT INTO neurotransmitter VALUES ('%s', '%s')" % (key, neurotransmitters[key]))
            except: 
                pass
        self.con.commit()

    def _getUniqueID(self):
        exist = True
        while exist:
            ID = str(uuid.uuid4())
            self.cur.execute("SELECT 1 FROM ID_table WHERE ID = '%s'" % ID)
            if self.cur.fetchone() == None: 
                exist = False
        return ID

    def getIDs(self, table):
        self.cur.execute("SELECT DISTINCT ID from ID_table WHERE table_name = '%s'" % table)
        return [x[0] for x in self.cur.fetchall()]

    def getNeurotransmitters(self):
        self.cur.execute("SELECT neurotransmitter from neurotransmitter")
        return [x[0] for x in self.cur.fetchall()]

    def _addState(self, table):
        neurotransmitters = self.getNeurotransmitters()
        ID = self._getUniqueID()
        self.cur.execute("INSERT INTO ID_table VALUES ('%s', '%s')" % (ID, table))
        for ntrans in neurotransmitters:
            self.cur.execute("INSERT INTO %s VALUES ('%s', '%s', %f)" % (table, ID, ntrans, 0.0))
        self.con.commit()
        return ID

    def addNeuron(self, n=1):
        IDList = []
        for i in range(int(n)):
            neurotransmitters = self.getNeurotransmitters()
            neuron_ID = self._getUniqueID()
            dendrite_state_ID = self._addState("dendrite_state")
            neuron_state_ID = self._addState("neuron_state")
            axon_state_ID = self._addState("axon_state")
            self.cur.execute("INSERT INTO ID_table VALUES ('%s', '%s')" % (neuron_ID, 'neuron_body'))
            self.cur.execute("INSERT INTO neuron_body VALUES ('%s', '%s', '%s')" % (neuron_ID, neuron_state_ID, axon_state_ID))
            self.cur.execute("INSERT INTO neuron_dendrite VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
            self.con.commit()
            IDs = (neuron_ID, dendrite_state_ID, neuron_state_ID, axon_state_ID)
            IDList.append(IDs)
        return IDList

    def addSynapse(self, n=1):
        synapse_state_IDs = []
        for i in range(int(n)):
            synapse_state_ID = self._addState("synapse_state")
            synapse_state_IDs.append(synapse_state_ID)
        return synapse_state_IDs

    def addDendrite(self, neuron_ID):
        neurotransmitters = self.getNeurotransmitters()
        dendrite_state_ID = self._addState("dendrite_state")
        self.cur.execute("INSERT INTO neuron_dendrite VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
        self.con.commit()
        return dendrite_state_ID

    def linkAxonSynapse(self, axon_state_ID, synapse_state_ID):
        self.cur.execute("INSERT INTO axon_synapse_link VALUES ('%s', '%s')" % (axon_state_ID, synapse_state_ID))
        self.con.commit()
        return [(axon_state_ID, synapse_state_ID)]

    def linkRandomAxonSynapse(self, n=1):
        axon_state_IDs = self.getIDs("axon_state")
        synapse_state_IDs = self.getIDs("synapse_state")
        linkages = []
        for i in range(int(n)):
            axon_state_ID = random.choice(axon_state_IDs)
            synapse_state_ID = random.choice(synapse_state_IDs)
            link = self.linkAxonSynapse(axon_state_ID, synapse_state_ID)
            linkages.append(link[0])
        return linkages

    def linkSynapseDendrite(self, synapse_state_ID, dendrite_state_ID):
        self.cur.execute("INSERT INTO synapse_dendrite_link VALUES ('%s', '%s')" % (synapse_state_ID, dendrite_state_ID))
        self.con.commit()
        return [(synapse_state_ID, dendrite_state_ID)]

    def linkRandomSynapseDendrite(self, n=1):
        dendrite_state_IDs = self.getIDs("dendrite_state")
        synapse_state_IDs = self.getIDs("synapse_state")
        linkages = []
        for i in range(int(n)):
            dendrite_state_ID = random.choice(dendrite_state_IDs)
            synapse_state_ID = random.choice(synapse_state_IDs)
            link = self.linkSynapseDendrite(synapse_state_ID, dendrite_state_ID)
            linkages.append(link[0])
        return linkages

    def tfSynapseDendrite(self, neuron_ID):
        pass

    def mfDendrite(self, neuron_ID):
        pass

    def tfDendriteNeuron(self, neuron_ID):
        pass

    def mfNeuron(self, neuron_ID):
        pass

    def tfNeuronAxon(self, neuron_ID):
        pass

    def mfAxon(self, neuron_ID):
        pass

    def tfAxonSynapse(self, neuron_ID):
        pass

    def mfSynapse(self, synapse_state_IDs):
        pass 
    
    def tfSynapseAxon(self, neuron_ID):
        pass

    def inputSignal(self, synapse_state_ID, signal_state):
        pass