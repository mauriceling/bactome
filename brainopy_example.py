from brainopy import brainopy

b = brainopy("brain.db")
b.logging = True

neurotransmitters = {"Ach": "acetylcholine",
                     "DA": "dopamine", 
                     "GLU": "glutamate",
                     "NE": "norepinephrine",
                     "5HT": "serotonin",
                     "GABA": "gamma-Aminobutyric acid"}
b.addNeurotransmitters(neurotransmitters)

neuronList = b.addNeuron(10)
# Neuron List = [(neuron_ID, dendrite_state_ID, neuron_state_ID, axon_state_ID)]
print("Neuron List: " + str(neuronList))

synapseList = b.addSynapse(10)
# Synapse List = [synapse_state_IDs]
print("Synapse List: " + str(synapseList))

linkages = b.linkRandomAxonSynapse(30)
# Axon-Synapse Linkages = [(axon_state_ID, synapse_state_ID)]
print("Axon-Synapse Linkages: " + str(linkages))

linkages = b.linkRandomSynapseDendrite(30)
# Synapse-Dendrite Linkages = [(synapse_state_ID, dendrite_state_ID)]
print("Synapse-Dendrite Linkages: " + str(linkages))

inputSignal = {"Ach": 0.11,
               "DA": 0.15, 
               "GLU": 0.21,
               "NE": 0.25,
               "5HT": 0.31,
               "GABA": 0.35}
for synapse_state_ID in synapseList:
    b.inputSignal(synapse_state_ID, inputSignal)
    print("Inserted input signal into synapse " + synapse_state_ID)

neuronList = b.getIDs("neuron_body")
# Neuron List = [neuron_ID]
print("Neurons to run brain: " + str(neuronList))

b.nameID(synapseList[3], "synapse4", "synapse number 4")
print("Label " + synapseList[3] + " as synapse4")

synapse4 = b.readNeurotransmitters("synapse4")
print("Original state of synapse4 = " + str(synapse4))

for cycle in range(1, 11):
    b.runBrain(neuronList)
    print("State of synapse4 after cycle " + str(cycle))
    print("By name: " + str(b.readNeurotransmitters("synapse4")))
    print("By ID: " + str(b.readNeurotransmitters(synapseList[3], "ID")))
    print("")

b.disconnectBrain()
