neuron_count = 1000
logging = True

import time
start_time = time.time()

from brainopy import brainopy

b = brainopy("brain.db")
b.logging = logging

neurotransmitters = {"Ach": "acetylcholine", "DA": "dopamine", "GLU": "glutamate","NE": "norepinephrine", "5HT": "serotonin", "GABA": "gamma-Aminobutyric acid"}
b.addNeurotransmitters(neurotransmitters)

_ = b.addNeuron(neuron_count)
neuronList = b.getIDs("neuron_body")
for neuron in neuronList: b.addDendrite(neuron)
for neuron in neuronList: b.addDendrite(neuron)

synapseList = b.addSynapse(2*neuron_count)
_ = b.linkRandomAxonSynapse(0.5*neuron_count)
_ = b.linkRandomSynapseDendrite(1.5*neuron_count)

inputSignal = {"Ach": 0.11, "DA": 0.15, "GLU": 0.21, "NE": 0.25, "5HT": 0.31, "GABA": 0.35}
for synapse_state_ID in synapseList: b.inputSignal(synapse_state_ID, inputSignal)

b.nameID(synapseList[0], "TestSynapse", "Synapse to measure")
print("Label " + synapseList[0] + " as TestSynapse")
synapse4 = b.readNeurotransmitters("TestSynapse")
print("Original state of TestSynapse = " + str(synapse4))

end_setup = time.time()
print("Setup time (s): " + str(int(end_setup - start_time)))

for cycle in range(1, 11):
    cycle_start = time.time()
    b.runBrain(neuronList)
    print("State of TestSynapse after cycle " + str(cycle))
    print(b.readNeurotransmitters("TestSynapse"))
    cycle_end = time.time()
    print("Cycle time (s): " + str(int(cycle_end - cycle_start)))
    print("")

b.disconnectBrain()
