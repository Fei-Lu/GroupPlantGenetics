library(tidyverse)
library(MOSAIC)
## model parameters, etc
# modelParameters <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/MOSAIC_RESULTS/E_2way_1-15_1-1_90_60_0.99_60.RData"
## local ancestry
# localAncestry <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/MOSAIC_RESULTS/localanc_E_2way_1-15_1-1_90_60_0.99_60.RData"
#
# mosaicInputData <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/mosaicData/"
#
# outFileAncestry0 <- "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne/ancestry_0/simulate_two_way_admixture_proportion0.4_genetration100_P1P2_divergence0.9Ne_LAI.txt"

args <- commandArgs(trailingOnly = T)
modelParameters <- args[1]
localAncestry <- args[2]
mosaicInputData <- args[3]
# outFileAncestry0 <- args[4]

load(modelParameters)
load(localAncestry)

## The above calculations and plots show local ancestry along evenly spaced gridpoints on recombination distances. You can get local ancestry estimates at the SNP positions using
local_pos=grid_to_pos(localanc,mosaicInputData,g.loc,chrnos)

## chr 1
chr1 = local_pos[[1]]

## matrix ancestry1 and ancestry2
chr1_ancestry0 = t(chr1[1,,])
chr1_ancestry1 = t(chr1[2,,])

## data_frame ancestry1 and ancestry2
# ancestry1 <- as_data_frame(chr1_ancestry1) %>%
#     rowwise() %>%
#     map_dfc(~ifelse(.x > 0.5, 0, 1))
#
# ancestry2 <- as_data_frame(chr1_ancestry2) %>%
#     rowwise() %>%
#     map_dfc(~ifelse(.x > 0.5, 1, 0))

# ancestry1 %>%
#     ggplot(aes(x=V1))+
#     geom_histogram()

# write_tsv(ancestry1, outFileAncestry0,col_names = F)
write_tsv(as_tibble(chr1_ancestry0),path=stdout())