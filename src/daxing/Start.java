package daxing;

import daxing.applets.Ne;
import daxing.common.Group;
import daxing.common.SubgenomeCombination;

public class Start {

    public static void main(String[] args) {
//        String vcfDir=args[0];
//        String bedDir=args[1];
//        String outDir=args[2];
//        String groupFile=args[3];
//        Group group=Group.valueOf(args[4]);
//        SubgenomeCombination subgenomeCombination=SubgenomeCombination.valueOf(args[5]);
//        String outFileAB=args[6];
//        String outFileD=args[7];
//        String genomeFile=args[8];
//        String logDir=args[9];
//        ScriptMethods.bulidSMC(vcfDir, bedDir, outDir, groupFile, group, subgenomeCombination, outFileAB, outFileD,
//                genomeFile, logDir);

        String vmap2InputDir=args[0];
        String taxaInfo_depthFile=args[1];
        String pseudoDiploidDir=args[2];
        Group group=Group.valueOf(args[3]);
        String ancestralDir=args[4];
        String ancestralVCFDir=args[5];
        Ne.syntheticPseudoDiploid(vmap2InputDir, taxaInfo_depthFile, pseudoDiploidDir, group, SubgenomeCombination.AB,45);
        Ne.syntheticPseudoDiploid(vmap2InputDir, taxaInfo_depthFile, pseudoDiploidDir, group, SubgenomeCombination.D,45);
        Ne.transformRefToAncestral(pseudoDiploidDir, ancestralDir, ancestralVCFDir);
    }
}