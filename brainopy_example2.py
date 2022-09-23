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


neuron_names = ["N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10"]
neuronList = [b.addNamedNeuron(name) for name in neuron_names]
print("Neuron List: " + str(neuronList))

synapse_names = ["S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8", "S9", "S10"]
synapseList = [b.addNamedSynapse(name) for name in synapse_names]
print("Synapse List: " + str(synapseList))