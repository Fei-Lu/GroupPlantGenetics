import random
import collections

import demes
import sys
import matplotlib.pyplot as plt
import msprime
import numpy as np
import os

# graph = demes.load(sys.argv[1], format="yaml")
# demography = msprime.Demography.from_demes(graph)
# print(demography.debug())
# sequence_length = sys.argv[2]
# random_seed = sys.argv[3]
# mutation_rate = sys.argv[4]
# outFile = sys.argv[5]

graphFile = sys.argv[1]
outFile_simulateGenotype = sys.argv[2]
outDir_tracts = sys.argv[3]
sequence_length = sys.argv[4]
random_seed = sys.argv[5]
mutation_rate = sys.argv[6]


# graphFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/001_deme/two_way_admixture_proportion0.1_genetration100.yaml"
# outFile_simulateGenotype = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/002_simulatedGenotype/simulate.vcf"
# outDir_tracts = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/003_simulatedTract"
# sequence_length = 100_000_000
# random_seed = 1
# mutation_rate = 7.1e-9

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
print(mts.num_mutations)
print("sim mutation end, and start write simulated vcf")
with open(outFile_simulateGenotype, "wt") as vcf_file:
    mts.write_vcf(vcf_file)

print("simulated vcf completed, and start write tracts")
sourceIndividuals = [indi for indi, pop in enumerate(ts.individuals_population) if pop == 4]
for i in sourceIndividuals:
    tracts = get_coalescing_tracts(ts, sourcePopID=2, targetIndividualID=i)
    with open(str(os.path.join(outDir_tracts, "individual"+str(i)+".tract.txt")), "w") as outFileTract:
        for row in tracts:
            outFileTract.write("\t".join([str(int(a)) for a in row]) + "\n")

print("tracts completed, Done")










