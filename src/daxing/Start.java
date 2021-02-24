package daxing;

import daxing.applets.Ne;

public class Start {

    public static void main(String[] args) {
//        String vcfDir=args[0];
//        String bedDir=args[1];
//        String outDir=args[2];
//        String taxaInfo_depthFile=args[3];
//        Group group=Group.valueOf(args[4]);
//        String outFileAB=args[5];
//        String outFileD=args[6];
//        String genomeFile=args[7];
//        String logDir=args[8];

//        String vcfDir="/Users/xudaxing/Desktop/Ne/004_refToAncestral_vcf";
//        String bedDir="/Users/xudaxing/Desktop/Ne/005_complementChr";
//        String outDir="/Users/xudaxing/Desktop/Ne/006_smc";
//        String taxaInfo_depthFile="/Users/xudaxing/Desktop/Ne/taxaInfo_depth.txt";
//        Group group= Group.Subcontinent;
//        String outFileAB="/Users/xudaxing/Desktop/Ne/007_out/ab.sh";
//        String outFileD="/Users/xudaxing/Desktop/Ne/007_out/d.sh";
//        String genomeFile="/Users/xudaxing/Desktop/Ne/chrID.Size.txt";
//        String logDir="/Users/xudaxing/Desktop/Ne/log";
//        int sampleSize=45;
//        ScriptMethods.bulidSMC(vcfDir, bedDir, outDir, taxaInfo_depthFile, group, outFileAB, outFileD,
//                genomeFile, logDir, 45);

        String smcDir=args[0];
        String refChrSMCDir=args[1];
        int popNumAB=Integer.parseInt(args[2]);
        int popNumD=Integer.parseInt(args[3]);
        Ne.mergeSMC(smcDir, refChrSMCDir, popNumAB, popNumD);
    }
}