package daxing;

import daxing.vmapII_1000.variantsAnnotation.VcfSubsetUtils;

public class Start {

    public static void main(String[] args) {
        String ancestralDir=args[0];
        String vcfDir=args[1];
        String outDir=args[2];
        String variantsNumSummaryFile=args[3];
        VcfSubsetUtils.extractVariantsWithAncestralState(ancestralDir, vcfDir, outDir, variantsNumSummaryFile);
    }
}