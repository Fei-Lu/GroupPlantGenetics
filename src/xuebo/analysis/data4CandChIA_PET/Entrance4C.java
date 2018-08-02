/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4CandChIA_PET;

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



public class Entrance4C{
    
    public static void main (String[] args) {
//        
//        String infileS = "/Users/xuebozhao/Documents/4Cdata/Library/addL23567.sorted.bed";
//        String outfileS = "/Users/xuebozhao/Documents/4Cdata/Library/addL23567_Filtering.txt";
////        new FirstStatistics(infileS,outfileS);
//        new FilteringBed(infileS,outfileS);
////        new BedUnSorted(infileS,outfileS);
//      

//        String infileS = "/Users/xuebozhao/Documents/4Cdata/Library/C4_L7.bed";
//        String outfileS = "/Users/xuebozhao/Documents/4Cdata/Library/C4_L7_BedUnSorted.txt";
////        new FirstStatistics(infileS,outfileS);
////        new FilteringBed(infileS,outfileS);
//        new BedUnSorted(infileS,outfileS);

//        String infileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/TAIR10_chr_all.fas";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/SplitTAIR10chr4.fas";
////        new FirstStatistics(infileS,outfileS);
//        new Cheakdata (infileS,outfileS);
        
//        String infileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/TAIRGenome/Chr5_rm.long.fas";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/TAIRGenome/Chr5ReplacePericentromere_rm.long.fas";
////        new FirstStatistics(infileS,outfileS);
//        new GetArabMaskedGenome (infileS,outfileS);
//        
//        String infileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/perimask/perimask5/positiveK5_hits_allNoadjusted001.bed";
//        String outfileS = "/Users/xuebozhao/Documents/LuLab/4C/4Cdata/Library/realproduct/perimask/perimask5/K5Noadjusted001ForCircos.bed";       
//        new GetCircosBed (infileS,outfileS);
        
        
        String infileS = "/Users/xuebozhao/Documents/allSIFTpredictions.txt";
        String outfileS = "/Users/xuebozhao/Documents/123.txt";       
        new GetCircosBed (infileS,outfileS);
        
    }
}

//public class Entrance4C{
//   
//        public static void main (String[] args) {
//        
//            String infileS = "/Users/xuebozhao/Documents/4Cdata/Library/1.txt";
//            String outfileS = "/Users/xuebozhao/Documents/4Cdata/Library/11111.txt";
//            new ReducedLibrary(infileS,"TTAA","AGATCT",outfileS);
//        }
//        
//}