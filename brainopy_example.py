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
for neuron_ID in neuronList:
    b.tfSynapseDendrite(neuron_ID)
    print("Dendrite state updated for neuron " + neuron_ID)

b.disconnectBrain()
