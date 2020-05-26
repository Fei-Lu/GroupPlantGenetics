/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

/**
 *
 * @author xuebozhao
 */
public class TTest {
    public static void main (String[] args) {
////    public TTest () {
//        FastaByte fb = new FastaByte("/Users/xuebozhao/Documents/LuLab/cpScore/Zea_mays.AGPv4.dna.toplevel.fa.gz");
//        GeneFeature gf = new GeneFeature ("/Users/xuebozhao/Documents/LuLab/cpScore/Zea_mays.AGPv4.38.kgf");
//        gf.writeCDSSequence(fb, "/Users/xuebozhao/Documents/LuLab/cpScore/testtextlalala.txt");
////        String seq1 = fb.getSeq(1);
////        seq1.subSequence(0, 0);


//        try{             
//            BufferedReader br = IOUtils.getTextReader("/Users/xuebozhao/Documents/LuLab/cpScore/temptest/tempCpscore/chr5CpCore.txt");
//            List<Double> listab = new ArrayList<Double>();
//            String temp = null;
//            while((temp = br.readLine()) != null){
//                double tempaaa = Double.parseDouble(temp);
//                listab.add(tempaaa);
//            }
//            
//            
//            
//            
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//            
//        }


//        String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/CHH_10/SitemethyTT10CHH";
////        String inFilePos = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/GeneFeaturePos5UTR.txt";
//           String inFilePos = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/CHH_10/GeneFeaturePos5UTR.txt";
//        String chr = "10-10";
////        String outFileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/CHH_10/Pos5UTRCHH.txt";
//        String outFileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/CHH_10/Pos5UTRCHH.txt";
//        GenePattern test = new GenePattern(infileS,chr,inFilePos,100, outFileS);
        
        
        
        String infileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/Cpscore/Cpscorepatten/CpCore.txt";
//        String inFilePos = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/GeneFeaturePos5UTR.txt";
           String inFilePos = "/Users/xuebozhao/Documents/LuLab/cpScore/data/Cpscore/Cpscorepatten/GeneFeaturePos5UTR.txt";
        String chr = "10-10";
//        String outFileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/CHH_10/Pos5UTRCHH.txt";
        String outFileS = "/Users/xuebozhao/Documents/LuLab/cpScore/data/Cpscore/Cpscorepatten/Pos5UTRCpS.txt";
        GenePattern test = new GenePattern(infileS,chr,inFilePos,100, outFileS);
        
        
    }
}
