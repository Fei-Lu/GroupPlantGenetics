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
//public class Entrance4C {
//    public static void main (String[] args) { 
//            String infileS1 = "/Users/xuebozhao/Documents/4C_data/cleandata2/C4-Library-Index2_H3VFNCCXY_L3_1.clean.fq.gz";
//            String infileS2 = "/Users/xuebozhao/Documents/4C_data/cleandata2/C4-Library-Index2_H3VFNCCXY_L3_2.clean.fq.gz";
//            String outfileS1 = "/Users/xuebozhao/Documents/4C_data/test/test.out";
//            String outfileS2 = "/Users/xuebozhao/Documents/4C_data/test/test_2.out";
//           
//            new Containing(infileS1 ,infileS2 , outfileS1,outfileS2);
//    }
//}



public class Entrance4C{
    
    public static void main (String[] args) {
        
        String infileS = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/containning_results/C4_L3.bed";
        String outfileS = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/containning_results/C4_L3_Filtering.txt";
        new Filtering(infileS,outfileS);
        
    }
}
