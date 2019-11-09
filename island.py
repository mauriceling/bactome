'''!
Island: A Simple Forward Simulation Tool for Population Genetics

Date created: 8th November 2019

License: GNU General Public License version 3 for academic or 
not-for-profit use only


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
import os
import random
import subprocess
import sys

try: 
    import fire
except ImportError:
    subprocess.check_call([sys.executable, '-m', 'pip', 
                           'install', 'fire'])
    import fire

def read_parameter_file(parameterfile, cmdline=True):
    """!
    Function to read a simulation parameter file containing the 
    allelic frequencies of all genes and print out the data for 
    checking.

    The simulation parameter file is a comma-delimited file, consisting 
    of a header row followed by row(s) of allelic frequencies. Each 
    allelic frequency is defined as <Gene Name>, [<Allelic Frequency 1>, 
    <Allelic Frequency 2>, ...]. Each Gene can have a different number 
    of alleles but the sum of all allelic frequencies for a specific 
    gene must add up to 1. A sample of simulation parameter file is 
    given as island_parameter.csv.

    Usage:
        
        python island.py readpf --parameterfile=island_parameter.csv

    @param parameterfile String: Relative or absolute path to the 
    simulation parameter file.
    @param cmdline Boolean: Flag to use this function as a command-line 
    function, which means results are not returned. Default = True 
    (results are not returned).
    """
    parameterfile = os.path.abspath(parameterfile)
    pfile = open(parameterfile, 'r').readlines()
    pfile = pfile[1:]
    pfile = [row[:-1] for row in pfile]
    pfile = [row.split(',') for row in pfile]
    pfile = [[x.strip() for x in row if x != ""] for row in pfile]
    gene_sequence = []
    pop_param = {}
    for gene in pfile:
        cumulativeAF = [float(af) for af in gene[1:]]
        caf = 0
        t_caf = []
        for af in cumulativeAF:
            caf = caf + af
            t_caf.append(caf)
        pop_param[str(gene[0])] = t_caf
        gene_sequence.append(str(gene[0]))
    if cmdline:
        print("Gene Name --> Cumulative Allelic Frequencies")
        for gene in pop_param:
            allelic_frequency = [str("%.5f" % af) 
                                 for af in pop_param[gene]]
            allelic_frequency = " | ".join(allelic_frequency)
            print("%s --> %s" % (gene, allelic_frequency))
    else:
        return (pop_param, gene_sequence)

def _generate_organism(pop_param, gene_sequence, ploidy):
    organism = []
    for i in range(int(ploidy)):
        ploid = []
        for gene in gene_sequence:
            freq = random.random()
            allele_bin = 0
            while freq >= pop_param[gene][allele_bin]:
                allele_bin = allele_bin + 1
            ploid.append(allele_bin + 1)
        organism.append(ploid)
    return organism

def generate_population(parameterfile, populationfile,
                        population_size=10, ploidy=2, 
                        generation_count=0):
    """!
    Function to generate the population file, which will be used as 
    input for simulation, from the simulation parameter file.

    Usage:
    
        python island.py gpop --populationfile=test_pop.pop --ploidy=2 --generation_count=0 --population_size=10 --parameterfile=island_parameter.csv

    @param parameterfile String: Relative or absolute path to the 
    simulation parameter file.
    @param populationfile String: Relative or absolute path for writing 
    out the population file for simulation.
    @param population_size Integer: Population size / number of 
    organisms to generate. Default = 10.
    @param ploidy Integer: Number of chromosome sets. Default = 2 
    (diploid).
    @param generation_count Integer: The base generation count. 
    Default = 0.
    """
    populationfile = os.path.abspath(populationfile)
    parameterfile = os.path.abspath(parameterfile)
    pop_file = open(populationfile, "w")
    (pop_param, gene_sequence) = read_parameter_file(parameterfile, False)
    allele_count = [str(len(pop_param[gene])) for gene in gene_sequence]
    stdout = '%s>%s' % ("|".join(gene_sequence), "|".join(allele_count))
    pop_file.write(stdout + "\n")
    print(stdout)
    for gene in gene_sequence:
        allelic_frequency = ["%.5f" % pop_param[gene][0]] + \
            ["%.5f" % (pop_param[gene][i] - pop_param[gene][i-1]) 
                for i in range(1, len(pop_param[gene]))]
        stdout = "A>%s>" % gene
        stdout = stdout + "|".join([str(af) for af in allelic_frequency])
        pop_file.write(stdout + "\n")
        print(stdout)
    for organism_count in range(int(population_size)):
        organism = _generate_organism(pop_param, gene_sequence, ploidy)
        stdout = "O>%s|%s|0|0>%s" % \
            (str(organism_count), str(int(generation_count)),
             ";".join(["|".join([str(allele) for allele in ploid]) 
                        for ploid in organism]))
        pop_file.write(stdout + "\n")
        print(stdout)
    pop_file.close()

def _simulation_writeout(filename, organisms, headerData):
    """!
    Private function called by simulate_<simulation type>() functions 
    to write out the corresponding population files for each generation.

    @param filename String: Relative or absolute path of the new
    population file.
    @param organisms Dictionary: Dictionary of organisms from 
    simulate_<simulation type>() functions to be written out.
    @param headerData List: Header data of population file (consisting 
    of gene list and allelic frequencies) from simulate_<simulation 
    type>() functions, to enable write out of simulated populations.
    """
    outputfile = open(filename, "w")
    for header in headerData:
        outputfile.write(header + "\n")
    for organism in organisms:
        genome = ["|".join(organisms[organism]['genome'][i]) 
                  for i in range(organisms[organism]['polyploid'])]
        stdout = "O>%s|%s|0|0>%s" % (str(organisms[organism]['organism']), 
                                     str(organisms[organism]['generation']), 
                                     str(";".join(genome)))
        outputfile.write(stdout + "\n")
    outputfile.close()

def simulate_simple(populationfile, generations, organisms, headerData):
    """!
    Function to perform simple simulation (simulation type = simple). 
    The features of this simulation types are:

        - assumes diploid (only the first 2 sets of chromosomes are used)
        - one random crossover per chromosome pair
        - no mutations
        - random mating with possibility of self-mating
        - mating only within generation

    @param populationfile String: Relative or absolute path of the 
    population file for simulation.
    @param generations Integer: Number of generations to simulate. 
    Note that generation count is not incremental from population, 
    which means that despite the generation in population file may be 
    50, the generation count in the results file will begin with 1. 
    @param organisms Dictionary: Dictionary of organisms from 
    simulate_population() function.
    @param headerData List: Header data of population file (consisting 
    of gene list and allelic frequencies) from simulate_population() 
    function, to enable write out of simulated populations.
    """
    for gen_count in range(int(generations)):
        gen_count = gen_count + 1
        outputfile = '.'.join([populationfile, str(gen_count)]) + '.pop'
        # Generating crossovers
        for organism in organisms:
            position = random.randint(0, len(organisms[organism]['genome'][0])-1)
            print("Generation count = %s; Organism = %s; Crossover position = %s" % \
                (str(gen_count), str(organism), str(position)))
            chromosomeA = organisms[organism]['genome'][0][0:position] + organisms[organism]['genome'][1][position:]
            chromosomeB = organisms[organism]['genome'][1][0:position] + organisms[organism]['genome'][0][position:]
            print("Orginal Polyploid A: %s" % "|".join(organisms[organism]['genome'][0]))
            print("Orginal Polyploid B: %s" % "|".join(organisms[organism]['genome'][1]))
            organisms[organism]['genome'][0] = chromosomeA
            organisms[organism]['genome'][1] = chromosomeB
            print("Crossed Polyploid A: %s" % "|".join(organisms[organism]['genome'][0]))
            print("Crossed Polyploid B: %s" % "|".join(organisms[organism]['genome'][1]))
        # Generating next generation
        organismList = list(organisms.keys())
        new_organisms = {}
        for i in range(len(organisms)):
            parentA = random.choice(organismList)
            parentB = random.choice(organismList)
            chromosomeA = organisms[parentA]['genome'][random.choice([0, 1])]
            chromosomeB = organisms[parentB]['genome'][random.choice([0, 1])]
            org = {'organism': str(i),
                   'generation': str(gen_count),
                   'parentA': str(parentA),
                   'parentB': str(parentB),
                   'polyploid': 2,
                   'genome': [chromosomeA, chromosomeB]}
            new_organisms[str(i)] = org
        _simulation_writeout(outputfile, new_organisms, headerData)
        organisms = new_organisms

def simulate_population(populationfile, simulation_type='simple',
                        generations=10):
    """!
    Function to simulate population over generations, given a 
    population. Allowable simulation types are simple. For more 
    description of the simulation types, please read documentation 
    in simulate_<simulation type>() function.

    Usage:
    
        python island.py simulate --populationfile=test_pop.pop --simulation_type=simple --generations=10

    @param populationfile String: Relative or absolute path of the 
    population file for simulation.
    @param simulation_type String: Type of simulation to run. Default 
    = simple.
    @param generations Integer: Number of generations to simulate. 
    Note that generation count is not incremental from population, 
    which means that despite the generation in population file may be 
    50, the generation count in the results file will begin with 1. 
    Default = 10.
    """
    populationfile = os.path.abspath(populationfile)
    inputfile = open(populationfile, "r").readlines()
    inputfile = [x[:-1] for x in inputfile]
    headerData = [x for x in inputfile if not x.startswith("O")]
    organismData = [x for x in inputfile if x.startswith("O")]
    organismData = [[x.split(">")[1].split("|"), 
                     x.split(">")[2].split(";")] 
                    for x in organismData]
    organismData = [[x[0], [ploid.split("|") for ploid in x[1]]] 
                    for x in organismData]
    organisms = {}
    for organism in organismData:
        org = {'organism': str(organism[0][0]),
               'generation': str(organism[0][1]),
               'parentA': str(organism[0][2]),
               'parentB': str(organism[0][3]),
               'genome': organism[1]}
        organisms[org['organism']] = org
    if simulation_type.upper() == 'SIMPLE':
        simulate_simple(populationfile, generations, 
                        organisms, headerData)


if __name__ == '__main__':
    exposed_functions = {
        'gpop': generate_population,
        'readpf': read_parameter_file,
        'simulate': simulate_population
        }
    fire.Fire(exposed_functions)
