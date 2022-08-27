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
        self.logging = False
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
        # CREATE TABLE statements
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
        self.cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS axon_synapse_link_index ON axon_synapse_link (axon_state_ID, synapse_state_ID)")
        self.cur.execute("CREATE UNIQUE INDEX IF NOT EXISTS synapse_dendrite_link_index ON synapse_dendrite_link (synapse_state_ID, dendrite_state_ID)")
        # CREATE VIEW statements
        self.cur.execute("CREATE VIEW IF NOT EXISTS synapse_dendrite (neuron_ID, dendrite_state_ID, synapse_state_ID) AS SELECT nd.ID, sdl.dendrite_state_ID, sdl.synapse_state_ID FROM neuron_dendrite nd INNER JOIN synapse_dendrite_link sdl WHERE nd.dendrite_state_ID = sdl.dendrite_state_ID")
        self.cur.execute("CREATE VIEW IF NOT EXISTS neuron (neuron_ID, dendrite_state_ID, neuron_state_ID, axon_state_ID) AS SELECT neuron_dendrite.ID, neuron_dendrite.dendrite_state_ID, neuron_body.neuron_state_ID, neuron_body.axon_state_ID FROM neuron_dendrite INNER JOIN neuron_body WHERE neuron_dendrite.ID = neuron_body.ID")
        self.con.commit()

    def disconnectBrain(self):
        self.con.commit()
        self.con.close()

    def logger(self, function, message):
        self.cur.execute("INSERT INTO log (function , message) VALUES ('%s', '%s')" % (function , message))
        self.con.commit()

    def addNeurotransmitters(self, neurotransmitters):
        for key in neurotransmitters:
            try: 
                self.cur.execute("INSERT INTO neurotransmitter (neurotransmitter, description) VALUES ('%s', '%s')" % (key, neurotransmitters[key]))
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
        self.cur.execute("INSERT INTO ID_table (ID, table_name) VALUES ('%s', '%s')" % (ID, table))
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
            self.cur.execute("INSERT INTO ID_table (ID, table_name) VALUES ('%s', '%s')" % (neuron_ID, 'neuron_body'))
            self.cur.execute("INSERT INTO neuron_body (ID, neuron_state_ID, axon_state_ID) VALUES ('%s', '%s', '%s')" % (neuron_ID, neuron_state_ID, axon_state_ID))
            self.cur.execute("INSERT INTO neuron_dendrite (ID, dendrite_state_ID) VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
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
        self.cur.execute("INSERT INTO neuron_dendrite (ID, dendrite_state_ID) VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
        self.con.commit()
        return dendrite_state_ID

    def linkAxonSynapse(self, axon_state_ID, synapse_state_ID):
        try:
            self.cur.execute("INSERT INTO axon_synapse_link (axon_state_ID, synapse_state_ID) VALUES ('%s', '%s')" % (axon_state_ID, synapse_state_ID))
            self.con.commit()
        except sqlite3.IntegrityError: pass
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
        try:
            self.cur.execute("INSERT INTO synapse_dendrite_link (synapse_state_ID, dendrite_state_ID) VALUES ('%s', '%s')" % (synapse_state_ID, dendrite_state_ID))
            self.con.commit()
        except sqlite3.IntegrityError: pass
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
        """
        Synapse to Dendrite Transfer Function
        """
        neurotransmitters = self.getNeurotransmitters()
        self.cur.execute("SELECT DISTINCT s.dendrite_state_ID, s.synapse_state_ID FROM synapse_dendrite s WHERE s.neuron_ID = '%s'" % neuron_ID)
        synapse_dendrite_List = [(x[0], x[1]) for x in self.cur.fetchall()]
        print("synapse_dendrite_List: " + str(synapse_dendrite_List))
        dendriteList = list(set([x[0] for x in synapse_dendrite_List]))
        for dendrite in dendriteList:
            print("process dendrite " + dendrite)
            dendrite_neuro = {}
            for n in neurotransmitters:
                dendrite_neuro[n] = []
            synapseList = list(set([x[1] for x in synapse_dendrite_List if x[0] == dendrite]))
            for synapse in synapseList:
                print("process dendritic synapse " + synapse)
                self.cur.execute("SELECT neurotransmitter, value FROM synapse_state WHERE ID = '%s'" % synapse)
                for state in [(x[0], x[1]) for x in self.cur.fetchall()]:
                    dendrite_neuro[state[0]] = dendrite_neuro[state[0]] + [float(state[1])]
            for n in neurotransmitters:
                dendrite_neuro[n] = str(sum(dendrite_neuro[n]) / len(dendrite_neuro[n]))
                self.cur.execute("UPDATE dendrite_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (dendrite_neuro[n], dendrite, n))
        self.con.commit()

    def mfDendrite(self, neuron_ID):
        """
        Dendrite Modulator
        """
        neurotransmitters = self.getNeurotransmitters()

    def tfDendriteNeuron(self, neuron_ID):
        """
        Dendrite to Neuron Transfer Function 
        """
        neurotransmitters = self.getNeurotransmitters()
        self.cur.execute("SELECT DISTINCT dendrite_state_ID, neuron_state_ID FROM neuron WHERE neuron_ID = '%s'" % neuron_ID)

    def mfNeuron(self, neuron_ID):
        """
        Neuron Modulator
        """
        neurotransmitters = self.getNeurotransmitters()

    def tfNeuronAxon(self, neuron_ID):
        """
        Neuron to Axon Transfer Function 
        """
        neurotransmitters = self.getNeurotransmitters()

    def mfAxon(self, neuron_ID):
        """
        Axon Modulator
        """
        neurotransmitters = self.getNeurotransmitters()

    def tfAxonSynapse(self, neuron_ID):
        """
        Axon to Synapse Transfer Function 
        """
        neurotransmitters = self.getNeurotransmitters()

    def mfSynapse(self, synapse_state_IDs):
        """
        Synapse Modulator
        """
        neurotransmitters = self.getNeurotransmitters()
    
    def tfSynapseAxon(self, neuron_ID):
        """
        Synapse to Axon Transfer Function 
        """
        neurotransmitters = self.getNeurotransmitters()

    def neuronFunction(self, neuron_ID):
        self.tfSynapseDendrite(neuron_ID)
        self.mfDendrite(neuron_ID)
        self.tfDendriteNeuron(neuron_ID)
        self.mfNeuron(neuron_ID)
        self.tfNeuronAxon(neuron_ID)
        self.mfAxon(neuron_ID)
        self.tfAxonSynapse(neuron_ID)
        
    def mtNeuronGrowth(self):
        """
        Neuronal Growth Function (NGF)
        """
        pass

    def mtSynapseGrowth(self):
        """
        Synaptic Growth Function (SGF)
        """
        pass

    def mtNeuronPrune(self):
        """
        Neuronal Prune Function (NPF)
        """
        pass

    def mtSynapsePrune(self):
        """
        Synaptic Prune Function (SGF)
        """
        pass

    def mtGlobal(self):
        """
        Global Maintenance Function (GMF)
        """
        pass

    def maintenanceFunction(self):
        self.mtNeuronGrowth()
        self.mtNeuronPrune()
        self.mtSynapseGrowth()
        self.mtSynapsePrune()

    def inputSignal(self, synapse_state_ID, signal_state):
        for neurotransmitter in signal_state:
            value = float(signal_state[neurotransmitter])
            self.cur.execute("UPDATE synapse_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (value, synapse_state_ID, neurotransmitter))
        self.con.commit()
