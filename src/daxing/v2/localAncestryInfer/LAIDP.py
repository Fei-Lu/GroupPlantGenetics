import random
import collections
import math
import demes
import sys
import matplotlib.pyplot as plt
import msprime
import numpy as np
import os

graphFile = sys.argv[1]
outFile_simulateGenotype = sys.argv[2]
taxaPopInfo = sys.argv[3]
outFile_recombinationMap = sys.argv[4]
outDir_tracts = sys.argv[5]
sequence_length = sys.argv[6]
random_seed = sys.argv[7]
mutation_rate = sys.argv[8]
recombination_rate = sys.argv[9]


# graphFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/000_graphs/01_two_way_admixture_proportion0.2_genetration100_P1P2_divergence0.9Ne.yaml"
# outFile_simulateGenotype = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/temp_genotype/simulate.vcf"
# taxaPopInfo = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/temp/taxaInfo.txt"
# outFile_recombinationMap = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/temp_genotype/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_recombinationMap.txt"
# outDir_tracts = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/temp_tracts"
# sequence_length = 1_000_000
# random_seed = 1
# mutation_rate = 7.1e-9
# recombination_rate = 1e-8

graph = demes.load(graphFile, format="yaml")
demography = msprime.Demography.from_demes(graph)
print(demography.debug())

ts = msprime.sim_ancestry(recombination_rate=1e-8,
                          sequence_length=sequence_length,
                          samples=[
                              msprime.SampleSet(30, ploidy=1, population='C'),
                              msprime.SampleSet(30, ploidy=1, population='D'),
                              msprime.SampleSet(30, ploidy=1, population='E'),
                          ],
                          demography = demography,
                          record_migrations=True,  # Needed for tracking segments.
                          random_seed=random_seed)

def get_migrating_tracts(ts):
    neanderthal_id = [p.id for p in ts.populations() if p.metadata['name'] == 'Neanderthal'][0]
    migrating_tracts = []
    # Get all tracts that migrated into the neanderthal population
    for migration in ts.migrations():
        if migration.dest == neanderthal_id:
            migrating_tracts.append((migration.left, migration.right))
    return np.array(migrating_tracts)

def get_coalescing_tracts(ts, sourcePopID, targetIndividualID):
    coalescing_tracts = []
    tract_left = None
    for tree in ts.trees():
        sourceIndividuals = [indi for indi, pop in enumerate(ts.individuals_population) if pop == sourcePopID]
        mrca_pops = []
        for m in sourceIndividuals:
            mrca = tree.mrca(targetIndividualID, m)
            node = ts.node(mrca)
            mrca_pops.append(node.population)
        left = tree.interval[0]
        if sourcePopID in mrca_pops and tract_left is None:
            # Start a new tract
            tract_left = left
        elif sourcePopID not in mrca_pops and tract_left is not None:
            # End the last tract
            coalescing_tracts.append((tract_left, left))
            tract_left = None
    if tract_left is not None:
        coalescing_tracts.append((tract_left, ts.sequence_length))
    return np.array(coalescing_tracts)

print("start sim mutation")
mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=random_seed)
print(str(mts.num_mutations) + " mutations")
print(str(mts.sites_position.size)+ " sites")
print("sim mutation end, and start write recombination map")

rateMap = msprime.RateMap.uniform(sequence_length, recombination_rate)
sites_position = mts.sites_position
recombination_rate = rateMap.get_rate(sites_position)
geneticsDistance = rateMap.get_cumulative_mass(sites_position)

with open(outFile_recombinationMap, "wt") as recombinationMap:
    recombinationMap.write("Position\tRate(cM/Mb)\tMap(cM)\n")
    for pos,rate,distance in zip(sites_position, recombination_rate, geneticsDistance):
        recombinationMap.write(str(math.ceil(pos))+"\t"+str(format(rate*1e+8, '.5f'))+"\t"+str(format(distance*100, '.5f')+"\n"))

print("write recombination map completed, and start write simulated vcf")
with open(outFile_simulateGenotype, "wt") as vcf_file:
    mts.write_vcf(vcf_file)

individualPop = ts.individuals_population
with open(taxaPopInfo, "wt") as taxaPopInfoFile:
    taxaPopInfoFile.write("Taxon\tPopulation\tPopulationName\n")
    for indi, pop in enumerate(ts.individuals_population):
        popName = ts.tables.populations[pop].metadata['name']
        taxaPopInfoFile.write("tsk_"+str(indi)+"\t"+str(pop)+"\t"+popName+"\n")

print("simulated vcf and taxaInfo completed, and start write tracts")
sourceIndividuals = [indi for indi, pop in enumerate(ts.individuals_population) if pop == 4]
for i in sourceIndividuals:
    tracts = get_coalescing_tracts(ts, sourcePopID=2, targetIndividualID=i)
    with open(str(os.path.join(outDir_tracts, "individual"+str(i)+".tract.txt")), "w") as outFileTract:
        for row in tracts:
            outFileTract.write("\t".join([str(int(a)) for a in row]) + "\n")

print("tracts completed, Done")










