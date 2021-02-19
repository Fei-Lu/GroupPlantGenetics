package daxing;

import daxing.applets.ScriptMethods;

public class Start {

    public static void main(String[] args) {
        String vcfDir="/Users/xudaxing/Desktop/Ne/004_refToAncestral_vcf";
        String bedDir="/Users/xudaxing/Desktop/Ne/005_complementChr";
        String outDir="/Users/xudaxing/Desktop/Ne/006_smc";
        String groupFile="/Users/xudaxing/Desktop/Ne/006_/groupSMC.txt";
        String distTaxon="/Users/xudaxing/Desktop/Ne/006_/popDistTaxon.txt";
        String outFileAB="/Users/xudaxing/Desktop/Ne/007_out/ab.sh";
        String outFileD="/Users/xudaxing/Desktop/Ne/007_out/d.sh";
        String genomeFile="/Users/xudaxing/Desktop/Ne/006_/genome.txt";
        String logDir="/Users/xudaxing/Desktop/Ne/log";
        ScriptMethods.bulidSMC(vcfDir,bedDir,outDir,groupFile,distTaxon,outFileAB,outFileD,genomeFile,logDir);
    }
}