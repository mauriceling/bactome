from brainopy import brainopy

b = brainopy("brain.db")

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

linkages = b.linkRandomAxonSynapse(10)
# Axon-Synapse Linkages = [(axon_state_ID, synapse_state_ID)]
print("Axon-Synapse Linkages: " + str(linkages))

linkages = b.linkRandomSynapseDendrite(10)
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

dendriteList = b.getIDs("dendrite_state")
for dendrite in dendriteList:
    b.tfSynapseDendrite(dendrite)
    print("Dendrite state updated - " + dendrite)

b.disconnectBrain()
