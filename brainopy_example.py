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
print("Neuron List: " + str(neuronList))

synapseList = b.addSynapse(10)
print("Synapse List: " + str(synapseList))

linkages = b.linkRandomAxonSynapse(10)
print("Axon-Synapse Linkages: " + str(linkages))

linkages = b.linkRandomSynapseDendrite(10)
print("Synapse-Dendrite Linkages: " + str(linkages))

b.disconnectBrain()
