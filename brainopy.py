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
        """ create a database connection to the SQLite database
            specified by the db_file
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
        self.cur.execute("CREATE INDEX IF NOT EXISTS ID_table_index (ID, table_name)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS ID_table_table (table_name)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_state_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS dendrite_state_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS axon_state_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS synapse_state_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_body_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS neuron_dendrite_ID (ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS axon_synapse_link_index (axon_state_ID, synapse_state_ID)")
        self.cur.execute("CREATE INDEX IF NOT EXISTS synapse_dendrite_link_index (synapse_state_ID, dendrite_state_ID)")
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

    def addNeuron(self):
        neurotransmitters = self.getNeurotransmitters()
        neuron_ID = self._getUniqueID()
        dendrite_state_ID = self._addState("dendrite_state")
        neuron_state_ID = self._addState("neuron_state")
        axon_state_ID = self._addState("axon_state")
        self.cur.execute("INSERT INTO ID_table VALUES ('%s', '%s')" % (ID, 'neuron_body'))
        self.cur.execute("INSERT INTO neuron_body VALUES ('%s', '%s', '%s')" % (neuron_ID, neuron_state_ID, axon_state_ID))
        self.cur.execute("INSERT INTO neuron_dendrite VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
        self.con.commit()
        return (neuron_ID, dendrite_state_ID, neuron_state_ID, axon_state_ID)

    def addSynapse(self):
        synapse_state_ID = self._addState("synapse_state")
        return synapse_state_ID

    def addDendrite(self, neuron_ID):
        neurotransmitters = self.getNeurotransmitters()
        dendrite_state_ID = self._addState("dendrite_state")
        self.cur.execute("INSERT INTO neuron_dendrite VALUES ('%s', '%s')" % (neuron_ID, dendrite_state_ID))
        self.con.commit()
        return dendrite_state_ID

    def linkAxonSynapse(self, axon_state_ID, synapse_state_ID):
        self.cur.execute("INSERT INTO axon_synapse_link VALUES ('%s', '%s')" % (axon_state_ID, synapse_state_ID))
        self.con.commit()

    def linkSynapseDendrite(self, synapse_state_ID, dendrite_state_ID):
        self.cur.execute("INSERT INTO synapse_dendrite_link VALUES ('%s', '%s')" % (synapse_state_ID, dendrite_state_ID))
        self.con.commit()
