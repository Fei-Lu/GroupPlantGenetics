package daxing.v2.localAncestryInfer.evaluation;

public interface LocalAncestry {

    /**
     * @return local ancestry matrix,
     * dim1 is different run, dim2 is different haplotype
     * dim3 is source population, dim4 is variants
     */
    double[][][][] extractLocalAncestry();
}
