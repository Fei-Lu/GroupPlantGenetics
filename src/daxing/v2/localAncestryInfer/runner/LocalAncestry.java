package daxing.v2.localAncestryInfer.runner;

public interface LocalAncestry {

    /**
     * @return local ancestry matrix,
     * dim1 is different run,
     * dim2 is different haplotype
     * dim3 is source population,  first all introgressed populations, then native population
     * dim4 is variants
     */
    double[][][][] extractLocalAncestry();
}
