package daxing;

import daxing.load.neutralSiteLoad.Go;

public class Start {

    public static void main(String[] args) {
//        String vcfDir="/Users/xudaxing/Desktop/test/001_vcfDir";
//        String ancestralDir="/Users/xudaxing/Desktop/test/002_ancestralDir";
//        String pgfFile="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/wheat_v1.1_Lulab.pgf";
//        String triadGeneNameFile="/Users/xudaxing/Desktop/test/BeiJing8.triad.txt.gz";
//        String outDir="/Users/xudaxing/Desktop/test/003_out";
//        String vmapIIGroupFile="/Users/xudaxing/Desktop/test/vmapGroup.txt";
//        Go.getDerivedCount(vcfDir, ancestralDir, pgfFile, triadGeneNameFile, outDir, vmapIIGroupFile);
        Go.getDerivedCount(args[0],args[1],args[2],args[3],args[4],args[5]);
    }
}