library(tidyverse)
library(MOSAIC)
## model parameters, etc
# modelParameters <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/MOSAIC_RESULTS/E_2way_1-15_1-1_90_60_0.99_60.RData"
## local ancestry
# localAncestry <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/MOSAIC_RESULTS/localanc_E_2way_1-15_1-1_90_60_0.99_60.RData"
#
# mosaicInputData <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/mosaicData/"
#
# outDir <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/ancestry_0/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne_LAI.txt"

args <- commandArgs(trailingOnly = T)
modelParameters <- args[1]
localAncestry <- args[2]
mosaicInputData <- args[3]
outDir <- args[4]

load(modelParameters)
load(localAncestry)

## The above calculations and plots show local ancestry along evenly spaced gridpoints on recombination distances. You can get local ancestry estimates at the SNP positions using
local_pos=grid_to_pos(localanc,mosaicInputData,g.loc,chrnos)

## chr 1
chr1 = local_pos[[1]]

for(i in 1:dim(chr1)[1]){
    ancestry = t(chr1[i,,])
    write_tsv(as_tibble(ancestry), file.path(outDir, str_c("ancestry",i-1,".txt")))
}