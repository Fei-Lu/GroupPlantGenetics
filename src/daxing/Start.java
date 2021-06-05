package daxing;

import daxing.hybrid.dection.HybridDection;

public class Start {

    public static void main(String[] args) {
        String inputDir="/Users/xudaxing/Desktop/test/001_vmap2";
        String outDir_depth="/Users/xudaxing/Desktop/test/002_outDir";
//        String outDir_Heterozygosity="/Users/xudaxing/Desktop/test/002_heter";
////        String outFile="/Users/xudaxing/Desktop/depth.txt.gz";
        HybridDection.getDepthHeterozygosity(inputDir, outDir_depth, 10_000_000, 5_000_000);
    }

}