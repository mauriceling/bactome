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

    Usage:
        
        python island.py readpf --parameterfile=island_parameter.csv
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
                        ploidy=2, generation_count=0,
                        population_size=10):
    """!

    Usage:
    
        python island.py gpop --populationfile=test_pop.pop --ploidy=2 --generation_count=0 --population_size=10 --parameterfile=island_parameter.csv
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

def _simulate_one_generation(populationfile, gen_count, organisms, headerData):
    filename = '.'.join([populationfile, str(gen_count)]) + '.pop'
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
    _simulation_writeout(filename, new_organisms, headerData)
    return new_organisms

def simulate_population(populationfile, simulation_type='simple',
                        generations=10):
    """!

    Usage:
    
        python island.py simulate --populationfile=test_pop.pop --simulation_type=simple --generations=10
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
    for gen_count in range(int(generations)):
        if simulation_type.upper() == 'SIMPLE':
            organisms = _simulate_one_generation(populationfile, 
                                                 gen_count+1, 
                                                 organisms, 
                                                 headerData)


if __name__ == '__main__':
    exposed_functions = {
        'gpop': generate_population,
        'readpf': read_parameter_file,
        'simulate': simulate_population
        }
    fire.Fire(exposed_functions)
