import random
import collections
import math
import demes
import sys
import demesdraw
import msprime
import numpy as np
import os
import matplotlib.pyplot as plt

graphFile = sys.argv[1]
outFile_model = sys.argv[2]
outFile_simulateGenotype = sys.argv[3]
taxaPopInfo = sys.argv[4]
outFile_recombinationMap = sys.argv[5]
outFile_tracts = sys.argv[6]
sequence_length = sys.argv[7]
random_seed = sys.argv[8]
mutation_rate = sys.argv[9]
recombination_rate = sys.argv[10]

## admixed population, native population, introgressed population
pops = sys.argv[11].split(",")
sampleSize = sys.argv[12].split(",")

# graphFile = "/Users/xudaxing/Desktop/LocalAncestry/demes/two_way_admixture_proportion0.1_genetration500.yaml"
# outFile_model = "/Users/xudaxing/Desktop/LocalAncestry/simulatedData/temp.pdf"
# outFile_simulateGenotype = "/Users/xudaxing/Desktop/LocalAncestry/simulatedData/simulate.vcf"
# taxaPopInfo = "/Users/xudaxing/Desktop/LocalAncestry/simulatedData/taxaInfo.txt"
# outFile_recombinationMap = "/Users/xudaxing/Desktop/LocalAncestry/simulatedData/recombinationMap.txt"
# outFile_tracts = "/Users/xudaxing/Desktop/LocalAncestry/simulatedData/localAncestryTract.txt"
# sequence_length = 10_000_000
# random_seed = 1
# mutation_rate = 7.1e-9
# recombination_rate = 1e-8

## admixed population, native population, introgressed population
# pops = "F,G,D,E".split(",")
# sampleSize = "10,10,10,10".split(",")
# pops = "E,D,C".split(",")
# sampleSize = "30,30,30".split(",")

graph = demes.load(graphFile, format="yaml")
demography = msprime.Demography.from_demes(graph)
# print(demography.debug())

demesdraw.tubes(graph)
# plt.show()
plt.savefig(fname=outFile_model,format="pdf")

sampleSetList = []
for pop,size in zip(pops, sampleSize):
    sampleSetList.append(msprime.SampleSet(int(size), ploidy=1, population=pop))

ts = msprime.sim_ancestry(recombination_rate=recombination_rate,
                          sequence_length=sequence_length,
                          # samples=[
                          #     msprime.SampleSet(3, ploidy=1, population='D'),
                          #     msprime.SampleSet(3, ploidy=1, population='E'),
                          #     msprime.SampleSet(3, ploidy=1, population='F'),
                          #     msprime.SampleSet(3, ploidy=1, population='G'),
                          # ],
                          samples=sampleSetList,
                          demography = demography,
                          record_migrations=True,  # Needed for tracking segments.
                          random_seed=random_seed)

def get_coalescing_tracts(ts, sourcePopID, admixedIndividualID):
    coalescing_tracts = []
    tract_left = None
    for tree in ts.trees():
        sourceIndividuals = [indi for indi, pop in enumerate(ts.individuals_population) if pop == sourcePopID]
        mrca_pops = []
        for m in sourceIndividuals:
            mrca = tree.mrca(admixedIndividualID, m)
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
        coalescing_tracts.append((tract_left, ts.sequence_length-1))
    return np.array(coalescing_tracts)

def getPopulationID(ts, populationName):
    populationID = None
    index = -1
    for row in ts.tables.populations:
        index += 1
        if (row.metadata["name"] == populationName):
            populationID = index
            break
    return populationID

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


admixedIndividuals = [indi for indi, pop in enumerate(ts.individuals_population) if pop == getPopulationID(ts, pops[0])]

introgressionDonorPopulation = pops[2:len(pops)]
with open(outFile_tracts, "w") as outFileTract:
    outFileTract.write("AdmixedIndividual\tIntrogressedPopulation\tstart\tEnd\n")
    for i in admixedIndividuals:
        tractsList = []
        for donor in introgressionDonorPopulation:
            tracts = get_coalescing_tracts(ts, sourcePopID=getPopulationID(ts, donor), admixedIndividualID=i)
            tractsList.append(tracts)
        for index,tract in enumerate(tractsList):
            for row in tract:
                outFileTract.write("tsk_"+str(i)+"\t"+introgressionDonorPopulation[index]+"\t"+"\t".join([str(int(a)+1) for a in row])+"\n")
print("tracts completed, Done")