'''!
Brainopy: A SQLite-Based Neural Network (Brain) Library

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
        self.cur.execute("CREATE TABLE IF NOT EXISTS log (ID integer primary key autoincrement, function text, message text)")
        # CREATE INDEX statements
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
        if self.logging: self.logger("connectBrain", "connectBrain")

    def disconnectBrain(self):
        self.con.commit()
        if self.logging: self.logger("disconnectBrain", "disconnectBrain")
        self.con.close()

    def logger(self, function, message):
        try: 
            self.cur.execute("INSERT INTO log (function , message) VALUES ('%s', '%s')" % (function , message))
        except sqlite3.OperationalError:
            print("INSERT INTO log (function , message) VALUES ('%s', '%s')" % (function , message))
        self.con.commit()

    def addNeurotransmitters(self, neurotransmitters):
        for key in neurotransmitters:
            try: 
                self.cur.execute("INSERT INTO neurotransmitter (neurotransmitter, description) VALUES ('%s', '%s')" % (key, neurotransmitters[key]))
                if self.logging: self.logger("addNeurotransmitters", "neurotransmitter=" + str(key) + "/value=" + str(value))
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
            if self.logging: self.logger("addNeuron", "1/neuron_ID=" + str(neuron_ID))
            dendrite_state_ID = self._addState("dendrite_state")
            if self.logging: self.logger("addNeuron", "2/new_dendrite_state/dendrite_state_ID=" + str(dendrite_state_ID))
            neuron_state_ID = self._addState("neuron_state")
            if self.logging: self.logger("addNeuron", "3/new_neuron_state/neuron_state_ID=" + str(neuron_state_ID))
            axon_state_ID = self._addState("axon_state")
            if self.logging: self.logger("addNeuron", "4/new_axon_state/axon_state_ID=" + str(axon_state_ID))
            self.cur.execute("INSERT INTO ID_table (ID, table_name) VALUES ('%s', '%s')" % (neuron_ID, 'neuron_body'))
            self.cur.execute("INSERT INTO neuron_body (ID, neuron_state_ID, axon_state_ID) VALUES ('%s', '%s', '%s')" % (neuron_ID, neuron_state_ID, axon_state_ID))
            self.cur.execute("INSERT INTO neuron_dendrite (ID, dendrite_state_ID) VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
            self.con.commit()
            if self.logging: self.logger("addNeuron", "5/insert_tables")
            IDs = (neuron_ID, dendrite_state_ID, neuron_state_ID, axon_state_ID)
            IDList.append(IDs)
        return IDList

    def addSynapse(self, n=1):
        synapse_state_IDs = []
        for i in range(int(n)):
            synapse_state_ID = self._addState("synapse_state")
            if self.logging: self.logger("addSynapse", "new_synapse_state/synapse_state_ID=" + str(synapse_state_ID))
            synapse_state_IDs.append(synapse_state_ID)
        return synapse_state_IDs

    def addDendrite(self, neuron_ID):
        neurotransmitters = self.getNeurotransmitters()
        dendrite_state_ID = self._addState("dendrite_state")
        if self.logging: self.logger("addDendrite", "1/new_dendrite_state/dendrite_state_ID=" + str(dendrite_state_ID))
        self.cur.execute("INSERT INTO neuron_dendrite (ID, dendrite_state_ID) VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
        self.con.commit()
        if self.logging: self.logger("addDendrite", "2/insert_tables")
        return dendrite_state_ID

    def linkAxonSynapse(self, axon_state_ID, synapse_state_ID):
        try:
            self.cur.execute("INSERT INTO axon_synapse_link (axon_state_ID, synapse_state_ID) VALUES ('%s', '%s')" % (axon_state_ID, synapse_state_ID))
            self.con.commit()
            if self.logging: self.logger("linkAxonSynapse", "new_axon_synapse/axon_state_ID=" + str(axon_state_ID) + "/synapse_state_ID=" + str(synapse_state_ID))
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
            if self.logging: self.logger("linkSynapseDendrite", "new_synapse_dendrite/dendrite_state_ID=" + str(dendrite_state_ID) + "/synapse_state_ID=" + str(synapse_state_ID))
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
        if self.logging: self.logger("tfSynapseDendrite", "1/get_links")
        dendriteList = list(set([x[0] for x in synapse_dendrite_List]))
        for dendrite in dendriteList:
            if self.logging: self.logger("tfSynapseDendrite", "2/process_dendrite/dendrite_state_ID=" + str(dendrite))
            dendrite_neuro = {}
            for n in neurotransmitters:
                dendrite_neuro[n] = []
            synapseList = list(set([x[1] for x in synapse_dendrite_List if x[0] == dendrite]))
            for synapse in synapseList:
                if self.logging: self.logger("tfSynapseDendrite", "3/process_dendritic_synapse/dendrite_state_ID=" + str(dendrite) + "/synapse_state_ID=" + str(synapse))
                self.cur.execute("SELECT neurotransmitter, value FROM synapse_state WHERE ID = '%s'" % synapse)
                for state in [(x[0], x[1]) for x in self.cur.fetchall()]:
                    dendrite_neuro[state[0]] = dendrite_neuro[state[0]] + [float(state[1])]
            for n in neurotransmitters:
                dendrite_neuro[n] = str(sum(dendrite_neuro[n]) / len(dendrite_neuro[n]))
                self.cur.execute("UPDATE dendrite_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (dendrite_neuro[n], dendrite, n))
                if self.logging: self.logger("tfSynapseDendrite", "4/update_dendrite_state/dendrite_state_ID=" + str(dendrite) + "/neurotransmitter=" + str(n) + "/value=" + str(dendrite_neuro[n]))
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
        dendrite_neuron_List = [(x[0], x[1]) for x in self.cur.fetchall()]
        if self.logging: self.logger("tfDendriteNeuron", "1/get_links")
        neuronList = list(set([x[1] for x in dendrite_neuron_List]))
        for neuron in neuronList:
            if self.logging: self.logger("tfDendriteNeuron", "2/process_neuron/neuron_state_ID=" + str(neuron))
            neuron_neuro = {}
            for n in neurotransmitters:
                neuron_neuro[n] = []
            dendriteList = list(set([x[0] for x in dendrite_neuron_List if x[1] == neuron]))
            for dendrite in dendriteList:
                if self.logging: self.logger("tfDendriteNeuron", "3/process_dendrite/neuron_state_ID=" + str(neuron) + "/dendrite_state_ID=" + str(dendrite))
                self.cur.execute("SELECT neurotransmitter, value FROM dendrite_state WHERE ID = '%s'" % dendrite)
                for state in [(x[0], x[1]) for x in self.cur.fetchall()]:
                    neuron_neuro[state[0]] = neuron_neuro[state[0]] + [float(state[1])]
            for n in neurotransmitters:
                neuron_neuro[n] = str(sum(neuron_neuro[n]) / len(neuron_neuro[n]))
                self.cur.execute("UPDATE neuron_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (neuron_neuro[n], neuron, n))
                if self.logging: self.logger("tfDendriteNeuron", "4/update_neuron_state/dendrite_state_ID=" + str(neuron) + "/neurotransmitter=" + str(n) + "/value=" + str(neuron_neuro[n]))
        self.con.commit()

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
        self.cur.execute("SELECT DISTINCT neuron_state_ID, axon_state_ID FROM neuron_body where ID = '%s'" % neuron_ID)
        neuron_axon = [(x[0], x[1]) for x in self.cur.fetchall()]
        if self.logging: self.logger("tfNeuronAxon", "1/get_link")
        if len(neuron_axon) > 1:
            print("BRAIN CORRUPTED!!! - MORE THAN ONE NEURON-AXON PAIR!")
            print(neuron_axon)
        else:
            neuron_state_ID = neuron_axon[0][0]
            axon_state_ID = neuron_axon[0][1]
        self.cur.execute("SELECT neurotransmitter, value FROM neuron_state WHERE ID = '%s'" % neuron_state_ID)
        stateList = [(x[0], x[1]) for x in self.cur.fetchall()]
        if self.logging: self.logger("tfNeuronAxon", "2/process_axon/neuron_state_ID=" + str(neuron_state_ID) + "/daxon_state_ID=" + str(axon_state_ID))
        for state in stateList:
            self.cur.execute("UPDATE axon_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (state[1], axon_state_ID, state[0]))
            if self.logging: self.logger("tfNeuronAxon", "3/update_axon_state/axon_state_ID=" + str(axon_state_ID) + "/neurotransmitter=" + str(state[0]) + "/value=" + str(state[1]))
        self.con.commit()

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
        if self.logging: self.logger("inputSignal", "1/input_signal/synapse_state_ID=" + str(synapse_state_ID))
        for neurotransmitter in signal_state:
            value = float(signal_state[neurotransmitter])
            self.cur.execute("UPDATE synapse_state SET value = '%s' WHERE ID = '%s' AND neurotransmitter = '%s'" % (value, synapse_state_ID, neurotransmitter))
            if self.logging: self.logger("inputSignal", "2/update_synapse_state/synapse_state_ID=" + str(synapse_state_ID) + "/neurotransmitter=" + str(neurotransmitter) + "/value=" + str(value))
        self.con.commit()
