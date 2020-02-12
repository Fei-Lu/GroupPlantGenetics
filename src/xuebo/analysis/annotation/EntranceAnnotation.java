/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import pgl.infra.dna.FastaByte;
import pgl.infra.utils.IOFileFormat;

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

//           String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/temp_end6_no0.txt";
//           String outfileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/siteScoreCpG.txt";
//           String outfileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/siteScoreCHH.txt";
//           String outfileS3 = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/siteScoreCHG.txt";           
//           new MethylationAnnotation (infileS , outfileS1, outfileS2 , outfileS3 );
           
//           String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/siteScoreCHH.txt";
//           String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/methylation/V4/allsiteCHHscore.txt";
//           new MethylationAnnotation (infileS,outfileS);
           
           
           String infileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforest/opentraining/APopenchrV4sortedL.text2.txt";
           String outfileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforest/opentraining/APopenchrV4sortedadd01.text2.txt";
           new MethylationAnnotation (infileS,outfileS);


//            String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/Methy_8/GeneFeaturePostranscript.txt";
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/Methy_8/GeneFeaturePostranscriptRemove.txt";
//            new RemoveStreamOverlap (infileS,outfileS);
//           
           
           
//
//            String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/Zea_mays.AGPv4.38_longestTrans.txt";
//           
//            //String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/temptest.txt";
//            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePosUpsream.txt";
//            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePosDownsream.txt";
//            String outfileS3 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePos5UTR.txt";
//            String outfileS4 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePosCDS.txt";
//            String outfileS5 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePosIntron.txt";
//            String outfileS6 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePos3UTR.txt";
//            String outfileS7 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/tempS/GeneFeaturePostranscript.txt";
//            
//            new SingleSingleGeneFeaturePos (infileS,outfileS1,outfileS2,outfileS3,outfileS4,outfileS5,outfileS6,outfileS7);
             //new generateScripts();
    }
} 