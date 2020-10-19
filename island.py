'''!
Island: A Simple Forward Simulation Tool for Population Genetics

Date created: 8th November 2019

License: GNU General Public License version 3 for academic or 
not-for-profit use only

Reference: Ling, MHT. 2019. Island: A Simple Forward Simulation Tool 
for Population Genetics. Acta Scientific Computer Sciences 1(2): 20-22.


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

######################################################################
# Section 1: Read file operations
######################################################################
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

def read_population_file(populationfile, cmdline=True):
    """!
    Function to read/prepare population file for simulation.

    Usage:

        python island.py readpop --populationfile=test_pop.pop

    @param populationfile String: Relative or absolute path of the 
    population file.
    @param cmdline Boolean: Flag to use this function as a command-line 
    function, which means results are not returned. Default = True 
    (results are not returned).
    """
    populationfile = os.path.abspath(populationfile)
    inputfile = open(populationfile, "r").readlines()
    inputfile = [x[:-1] for x in inputfile]
    geneData = inputfile[0]
    alleleData = [x for x in inputfile if x.startswith("A")]
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
    if cmdline:
        print("Gene Names : %s" % geneData.split(">")[0])
        print("Number of Alleles : %s" % geneData.split(">")[1])
        print("")
        print("Allelic Frequencies ...")
        for allele in alleleData:
            print("%s : %s" % (allele.split(">")[1], 
                               allele.split(">")[2]))
        print("")
        print("Organism Data ...")
        for org in organisms:
            print("Organism = %s; Generation = %s; ParentA = %s; ParentB = %s" % \
                (str(organisms[org]['organism']),
                 str(organisms[org]['generation']),
                 str(organisms[org]['parentA']),
                 str(organisms[org]['parentB'])))
            for i in range(len(organisms[org]['genome'])):
                print("Polyloid %i = %s" % \
                    (i, "|".join(organisms[org]['genome'][i])))
    else:
        return (geneData, alleleData, organismData, organisms)

######################################################################
# Section 2: Generate population from parameters operations
######################################################################
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
    
        python island.py gpop --populationfile=test_pop --ploidy=2 --generation_count=0 --population_size=10 --parameterfile=island_parameter.csv

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
    # print(stdout)
    print("Number of Genes = %i" % len(gene_sequence))
    for gene in gene_sequence:
        allelic_frequency = ["%.5f" % pop_param[gene][0]] + \
            ["%.5f" % (pop_param[gene][i] - pop_param[gene][i-1]) 
                for i in range(1, len(pop_param[gene]))]
        stdout = "A>%s>" % gene
        stdout = stdout + "|".join([str(af) for af in allelic_frequency])
        pop_file.write(stdout + "\n")
        # print(stdout)
    for organism_count in range(int(population_size)):
        organism = _generate_organism(pop_param, gene_sequence, ploidy)
        stdout = "O>%s|%s|0|0>%s" % \
            (str(organism_count), str(int(generation_count)),
             ";".join(["|".join([str(allele) for allele in ploid]) 
                        for ploid in organism]))
        pop_file.write(stdout + "\n")
        # print(stdout)
        if (organism_count + 1) % 100 == 0:
            print("%s organisms generated" % str(organism_count + 1))
    print("Total %s organisms generated" % str(organism_count + 1))
    pop_file.close()

######################################################################
# Section 3: Simulation operations, given population
######################################################################
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
        stdout = "O>%s|%s|%s|%s>%s" % \
            (str(organisms[organism]['organism']), 
             str(organisms[organism]['generation']), 
             str(organisms[organism]['parentA']),
             str(organisms[organism]['parentB']),
             str(";".join(genome)))
        outputfile.write(stdout + "\n")
    outputfile.close()

def simulate_simple(populationfile, generations, organisms, 
                    population_size, headerData):
    """!
    Function to perform simple simulation (simulation type = simple). 
    The features of this simulation types are:

        - assumes diploid (only the first 2 sets of chromosomes are used)
        - one random crossover per chromosome pair
        - crossover is generated prior to mating to simulate random haploid
        - no mutations
        - random mating without possibility of self-mating
        - mating only within generation
        - population size may be changed

    Required options are:

        - populationfile
        - simulation_type
        - population_size
        - generations

    Usage:

        python island.py simulate --populationfile=test_pop --simulation_type=simple --population_size=10 --generations=10

    @param populationfile String: Relative or absolute path of the 
    population file for simulation.
    @param generations Integer: Number of generations to simulate. 
    Note that generation count is not incremental from population, 
    which means that despite the generation in population file may be 
    50, the generation count in the results file will begin with 1. 
    @param organisms Dictionary: Dictionary of organisms from 
    simulate_population() function.
    @param population_size Integer: Population size. If the population 
    size in the population file is lesser than the required population 
    size, this function will expand the population size to the required 
    population size at the first generation.
    @param headerData List: Header data of population file (consisting 
    of gene list and allelic frequencies) from simulate_population() 
    function, to enable write out of simulated populations.
    """
    for gen_count in range(int(generations)):
        gen_count = gen_count + 1
        outputfile = '.'.join([populationfile, str(gen_count)])
        # Generating next generation
        organismList = list(organisms.keys())
        new_organisms = {}
        i = 0
        while i < int(population_size):
            parentA = 0
            parentB = 0
            while parentA == parentB:
                parentA = random.choice(organismList)
                parentB = random.choice(organismList)
            chromosomeA1 = organisms[parentA]['genome'][0]
            chromosomeA2 = organisms[parentA]['genome'][1]
            chromosomeB1 = organisms[parentB]['genome'][0]
            chromosomeB2 = organisms[parentB]['genome'][1]
            position = random.randint(0, len(chromosomeA1))
            genomeA = [chromosomeA1[0:position] + chromosomeA2[position:],
                       chromosomeA2[0:position] + chromosomeA1[position:]]
            genomeB = [chromosomeB1[0:position] + chromosomeB2[position:],
                       chromosomeB2[0:position] + chromosomeB1[position:]]
            org = {'organism': str(i),
                   'generation': str(gen_count),
                   'parentA': str(parentA),
                   'parentB': str(parentB),
                   'polyploid': 2,
                   'genome': [genomeA[random.randint(0, 1)], 
                              genomeB[random.randint(0, 1)]]}
            new_organisms[str(i)] = org
            if i+1 % 100 == 0:
                print("%s organisms produced" % str(i))
            i = i + 1
        print("Total %s organisms produced" % str(i))
        _simulation_writeout(outputfile, new_organisms, headerData)
        organisms = new_organisms

def simulate_population(populationfile, 
                        population_size, 
                        simulation_type='simple',
                        generations=10):
    """!
    Function to simulate population over generations, given a 
    population. Allowable simulation types are simple. For more 
    description of the simulation types and the corresponding options 
    required for each simulation (as not all simulation types will 
    require the same set of options), please read documentation in 
    simulate_<simulation type>() function.

    Usage:
    
        python island.py simulate --populationfile=test_pop --simulation_type=simple --population_size=10 --generations=10

    @param populationfile String: Relative or absolute path of the 
    population file for simulation.
    @param population_size Integer: Population size / number of 
    organisms to generate. Default = 10.
    @param simulation_type String: Type of simulation to run. Default 
    = simple.
    @param generations Integer: Number of generations to simulate. 
    Note that generation count is not incremental from population, 
    which means that despite the generation in population file may be 
    50, the generation count in the results file will begin with 1. 
    Default = 10.
    """
    (geneData, alleleData, organismData, organisms) = \
        read_population_file(populationfile, False)
    headerData = [geneData] + alleleData
    if simulation_type.upper() == 'SIMPLE':
        simulate_simple(populationfile, generations, organisms, 
                        population_size, headerData)

######################################################################
# Section 4: Analyze population operations
######################################################################
def _generate_expected_allelic_counts(alleleData, pop_size, ploidy):
    """!
    Private function to generate expected allelic counts from allelic 
    frequencies.

    @param alleleData List: List of allele data from population file.
    @param pop_size Integer: Population size.
    @param ploidy Integer: Number of chromosome sets.
    """
    alleleData = [[x.split(">")[1].strip(), x.split(">")[2].strip()] 
                  for x in alleleData]
    allele_counts = {}
    for allele in alleleData:
        counts = [float(x) * pop_size * ploidy 
                  for x in allele[1].split("|")]
        allele_counts[allele[0]] = counts
    # print(allele_counts)
    return allele_counts

def _generate_observed_allelic_counts(organisms, geneData):
    """!
    Private function to tabulate observed allelic counts from population.

    @param organisms Dictionary: Dictionary of organisms.
    @param geneData String: Gene Names and the corresponding number of 
    alleles, which is the first row of population file.
    """
    geneData = geneData.split(">")
    geneNames = [x.strip() for x in geneData[0].split("|")]
    geneAlleles = [x.strip() for x in geneData[1].split("|")]
    geneData = [x for x in zip(geneNames, geneAlleles)]
    allele_counts = {}
    for pair in geneData:
        allele_counts[pair[0]] = [0] * int(pair[1])
    # print(allele_counts)
    for position in range(len(geneNames)):
        alleles = [int(organisms[org]['genome'][ploid][position])
                   for org in organisms
                        for ploid in range(len(organisms[org]['genome']))]
        # print(alleles)
        geneName = geneNames[position]
        # print(geneName, allele_counts[geneName])
        for allele in alleles:
            allele_counts[geneName][allele-1] = allele_counts[geneName][allele-1] + 1
    # print(allele_counts)
    return allele_counts

def tabulate_allelic_counts(populationfile, ploidy=2):
    """!
    Function to tabulate the expected and actual allelic counts from 
    the population and reports Chi-Square statistic.

    Usage:
    
        python island.py tabulateCount --ploidy=2 --populationfile=test_pop

    @param populationfile String: Relative or absolute path of the 
    population file to tabulate.
    @param ploidy Integer: Number of chromosome sets. Default = 2 
    (diploid).
    """
    (geneData, alleleData, organismData, organisms) = \
        read_population_file(populationfile, False)
    pop_size = len(organisms)
    exp_allele_counts = _generate_expected_allelic_counts(alleleData, 
                                                          pop_size, 
                                                          ploidy)
    obs_allele_counts = _generate_observed_allelic_counts(organisms, 
                                                          geneData)
    chiSq = 0
    df = -1
    print("Gene Name : Allele : Expected Count : Observed Count")
    for gene in exp_allele_counts:
        for allele in range(len(exp_allele_counts[gene])):
            print("%s : %s : %s : %s" % \
                (str(gene), str(allele),
                 str(exp_allele_counts[gene][allele]),
                 str(obs_allele_counts[gene][allele])))
            df = df + 1
            chiSq = chiSq + \
                (((obs_allele_counts[gene][allele] - \
                    exp_allele_counts[gene][allele]) ** 2)  / \
                exp_allele_counts[gene][allele])
    print("Chi Square Statistic : %s" % str(chiSq))
    print("Degrees of Freedom : %i" % df)
    print("Normalized Chi Square Statistic : %s" % str(chiSq / (df + 1)))

######################################################################
# Section 5: Utility operations
######################################################################
def combine_populations(populationfile1, populationfile2, outputfile):
    """!
    Function to combine the organisms in 2 population files into a 
    single population file. This function assumes that the gene names 
    and allelic frequency table are identical in both population files.

    Usage:

        python island.py combinepop --populationfile1=test_pop.2 --populationfile2=test_pop.3 --outputfile=test_pop.comb

    @param populationfile1 String: Relative or absolute path of the 
    first population file to combine.
    @param populationfile2 String: Relative or absolute path of the 
    second population file to combine.
    @param putputfile String: Relative or absolute path for writing 
    out the combined population file.
    """
    populationfile1 = os.path.abspath(populationfile1)
    populationfile2 = os.path.abspath(populationfile2)
    outputmapfile = outputfile + ".map"
    outputfile = os.path.abspath(outputfile)
    outputfile = open(outputfile, "w")
    outputmapfile = os.path.abspath(outputmapfile)
    outputmapfile = open(outputmapfile, "w")
    print("Read population file %s" % str(populationfile1))
    (geneData1, alleleData1, organismData1, _) = \
        read_population_file(populationfile1, False)
    print("Read population file %s" % str(populationfile2))
    (_, _, organismData2, _) = \
        read_population_file(populationfile2, False)
    print("Combine populations")
    outputfile.write(geneData1 + "\n")
    # print(geneData1)
    for allele in alleleData1:
        outputfile.write(allele + "\n")
    # print(alleleData1)
    count = 0
    for org in organismData1:
        genome = ["|".join(org[1][i]) for i in range(len(org[1]))]
        genome = ";".join(genome)
        mapdata = "%s:%s>%s" % (populationfile1, str(org[0][0]), str(count))
        print(mapdata)
        outputmapfile.write(mapdata + "\n")
        org[0][0] = str(count)
        organism = "O>" + "|".join(org[0])
        organismData = ">".join([organism, genome])
        outputfile.write(organismData + "\n")
        # print(organismData)
        count = count + 1
    for org in organismData2:
        genome = ["|".join(org[1][i]) for i in range(len(org[1]))]
        genome = ";".join(genome)
        mapdata = "%s:%s>%s" % (populationfile2, str(org[0][0]), str(count))
        print(mapdata)
        outputmapfile.write(mapdata + "\n")
        org[0][0] = str(count)
        organism = "O>" + "|".join(org[0])
        organismData = ">".join([organism, genome])
        outputfile.write(organismData + "\n")
        # print(organismData)
        count = count + 1
    print("Total number of combined organisms = %s" % \
        str(count))
    outputmapfile.close()
    outputfile.close()

def randomly_select_population(populationfile, outputfile, n):
    """!
    Function to randomly select n organisms from a population file into 
    another population file.

    Usage:

        python island.py random --populationfile=test_pop --outputfile=test_pop.rand --n=30

    @param populationfile String: Relative or absolute path of the 
    population file to randomly select population from.
    @param putputfile String: Relative or absolute path for writing 
    out the randomly selected population file.
    @param n Integer: Number of organisms to select.
    """
    populationfile = os.path.abspath(populationfile)
    outputfile = os.path.abspath(outputfile)
    print("Sample %s organisms from %s into %s" % \
        (str(int(n)), populationfile, outputfile))
    outputfile = open(outputfile, "w")
    (geneData, alleleData, organismData, _) = \
        read_population_file(populationfile, False)
    outputfile.write(geneData + "\n")
    # print(geneData)
    for allele in alleleData:
        outputfile.write(allele + "\n")
    # print(alleleData)
    if int(n) > len(organismData):
        n = len(organismData)
    random.shuffle(organismData)
    for org in organismData[:int(n)]:
        genome = ["|".join(org[1][i]) for i in range(len(org[1]))]
        genome = ";".join(genome)
        organism = "O>" + "|".join(org[0])
        organismData = ">".join([organism, genome])
        outputfile.write(organismData + "\n")
        # print(organismData)
    outputfile.close()
######################################################################
# Section 6: Command-line executor
######################################################################
if __name__ == '__main__':
    exposed_functions = {
        'combinepop': combine_populations,
        'gpop': generate_population,
        'random': randomly_select_population,
        'readpf': read_parameter_file,
        'readpop': read_population_file,
        'simulate': simulate_population,
        'tabulateCount': tabulate_allelic_counts
        }
    fire.Fire(exposed_functions)
