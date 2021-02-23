package daxing;

import daxing.applets.Ne;
import daxing.common.vmap2Group.Group;
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

//        String vcfDir="/Users/xudaxing/Desktop/Ne/004_refToAncestral_vcf";
//        String bedDir="/Users/xudaxing/Desktop/Ne/005_complementChr";
//        String outDir="/Users/xudaxing/Desktop/Ne/006_smc";
//        String taxaInfo_depthFile="/Users/xudaxing/Desktop/Ne/taxaInfo_depth.txt";
//        Group group=Group.Subcontinent;
//        String outFileAB="/Users/xudaxing/Desktop/Ne/007_out/ab.sh";
//        String outFileD="/Users/xudaxing/Desktop/Ne/007_out/d.sh";
//        String genomeFile="/Users/xudaxing/Desktop/Ne/chrID.Size.txt";
//        String logDir="/Users/xudaxing/Desktop/Ne/log";
//        int sampleSize=45;
//        ScriptMethods.bulidSMC(vcfDir, bedDir, outDir, taxaInfo_depthFile, group, outFileAB, outFileD,
//                genomeFile, logDir, sampleSize);

        String vmap2InputDir=args[0];
        String taxaInfo_depthFile=args[1];
        String pseudoDiploidDir=args[2];
        Group group= Group.valueOf(args[3]);
        String ancestralDir=args[4];
        String ancestralVCFDir=args[5];
        Ne.syntheticPseudoDiploid(vmap2InputDir, taxaInfo_depthFile, pseudoDiploidDir, group, SubgenomeCombination.AB,45);
        Ne.syntheticPseudoDiploid(vmap2InputDir, taxaInfo_depthFile, pseudoDiploidDir, group, SubgenomeCombination.D,45);
        Ne.transformRefToAncestral(pseudoDiploidDir, ancestralDir, ancestralVCFDir);
    }
}