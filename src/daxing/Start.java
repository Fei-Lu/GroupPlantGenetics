package daxing;

import daxing.common.VCF;

public class Start {

    public static void main(String[] args) {
//        String chrSNPNumFile=args[0];
//        String vcfInputDir=args[1];
//        int variantsNum=Integer.parseInt(args[2]);
//        String pgfFileGeneHC=args[3];
//        String ancestraDir=args[4];
//        String outDir=args[5];
//        VcfSubsetUtils.go(chrSNPNumFile, vcfInputDir, variantsNum, pgfFileGeneHC, ancestraDir, outDir);
        VCF.fastMergeVCFtoLineage(args[0],args[1]);
    }
}