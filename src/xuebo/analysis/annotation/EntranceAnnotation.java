/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import format.dna.FastaByte;
import utils.IOFileFormat;

/**
 *
 * @author xuebozhao
 */
//public class EntranceAnnotation {    
//    public static void main (String[] args) { 
////        GeneFeature a = new GeneFeature (args[0]);
////        a.writeFile("lala.txt");
//////          RangeAttribute b = new RangeAttribute(arg[0]);
//////          b.writeTextFile("lalala");
//            String infileS = "D:\\FeiLub\\UniquenessScore\\maizeV3\\AP.bfthresh1.1.MNaseHS.Ranges.dat";
//            String outfileS = "D:\\FeiLub\\UniquenessScore\\maizeV3\\lalala.bed";
//            new DatToBed (infileS,outfileS);
////            c.readFile("lalalabed");
//    }
//} 


public class EntranceAnnotation {    
    public static void main (String[] args) { 
//        GeneFeature a = new GeneFeature (args[0]);
//        a.readFromMaizeGFF("/Users/xuebozhao/Documents/LuLab/cpScore/Zea_mays.AGPv4.36.gff3.gz");
//        a.writeFile("lalalaV4.txt");
//          RangeAttribute b = new RangeAttribute(arg[0]);
//          b.writeTextFile("lalala");

//            String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/temptest/GetCDS.txt";
//            
//            
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/temptest/GetCDSSingleGeneFeaturePos.txt";
//            
////            FastaByte genomef = new FastaByte(infileS);
//            new SingleGeneFeaturePos (infileS,outfileS);
//           String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/test";
//           String outfileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/test.CpG.bed.gz";
//           String outfileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/test.CHH.bed.gz";
//           String outfileS3 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/test.CHG.bed.gz";           
//           new MethylationAnnotation (infileS , outfileS1, outfileS2 , outfileS3 );
           
           

            String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/Zea_mays.AGPv4.38_longestTrans.txt";
                       
            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePosUpsream.txt";
            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePosDownsream.txt";
            String outfileS3 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePos5UTR.txt";
            String outfileS4 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePosCDS.txt";
            String outfileS5 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePosIntron.txt";
            String outfileS6 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/SingleGeneFeaturePos/biubiubiu/GeneFeaturePos3UTR.txt";
            
            new SingleGeneFeaturePos (infileS,outfileS1,outfileS2,outfileS3,outfileS4,outfileS5,outfileS6);

    }
} 