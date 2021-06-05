package daxing;

import daxing.hybrid.dection.HybridDection;

public class Start {

    public static void main(String[] args) {
        String inputDir="/Users/xudaxing/Desktop/001_fromVMap2.0_singleChr0.001";
        String outDir_depth="/Users/xudaxing/Desktop/test/001_depth";
        String outDir_Heterozygosity="/Users/xudaxing/Desktop/test/002_heter";
//        String outFile="/Users/xudaxing/Desktop/depth.txt.gz";
        HybridDection.getDepthHeterozygosity(inputDir, outDir_depth, outDir_Heterozygosity);
    }

}