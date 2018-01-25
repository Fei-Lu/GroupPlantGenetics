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
//            String infileS1 = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/test/1.txt";
//            String infileS2 = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/test/2.txt";
//            String outfileS1 = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/test/test1111.out.gz";
//            String outfileS2 = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/test/test2222.out.gz";
//           
//            new Containing(infileS1 ,infileS2 , outfileS1,outfileS2);
//    }
//}



//public class Entrance4C{
//    
//    public static void main (String[] args) {
        
//        String infileS = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/containning_results/C4_L3.bed";
//        String outfileS = "/Users/xuebozhao/Documents/4Cdata/4C_data_1/containning_results/C4_L3_Filtering.txt";
//        new Filtering(infileS,outfileS);


//        String a  = "agatct";
//        int b = a.indexOf("atct");
//        System.out.println( b );
        
//    }
//}

public class Entrance4C{
   
        public static void main (String[] args) {
        
            String infileS = "/Users/xuebozhao/Documents/4Cdata/Library/1.txt";
            String outfileS = "/Users/xuebozhao/Documents/4Cdata/Library/11111.txt";
            new ReducedLibrary(infileS,"TTAA","AGATCT",outfileS);
        }
        
}