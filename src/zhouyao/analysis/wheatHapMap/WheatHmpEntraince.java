/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

/**
 *
 * @author feilu
 */
public class WheatHmpEntraince {
    /**
     * @param args the command line arguments
     * -m: model, LogS: statistic of log files;
     * -i: input file name
     * -o: output file name
     * -s: size, unit is K reads; 5 represent 5k reads;
     */
    public static void main(String[] args) {
        int len = args.length;
        String model = "LogS";
        String inFile = "";
        String outFile="";
        int ExSize = 5*1000;
        for (int i = 0; i < len; i++){
            if (null != args[i])switch (args[i]) {
                case "-m":
                    model = args[i+1];
                    i++;
                    break;
                case "-i":
                    inFile = args[i+1];
                    i++;
                    break;
                case "-o":
                    outFile = args[i+1];
                    i++;
                    break;
                case "-s":
                    double ExSize1 = Double.valueOf(args[i+1]);
                    ExSize = (int) (ExSize1 * 1000);
                    i++;
                    break;
                default:
                    break;
            }
        }
        if(outFile.equals("")){
            String[] temp = inFile.split("/");
            outFile = temp[temp.length -1].split("\\.")[0];
        }
        
        if(model.equals("LogS")){
            new LogReadStatistic(inFile);
        }
        if(model.equals("ExFastq")){
            new ExtractionOfFastq(inFile,outFile,ExSize);
        }
        if(model.equals("GenoStat")){
            new GenomeStatistic(inFile,outFile);
        }
        if(model.equals("GenoDist")){
            System.out.println("GenoDist is working...");
            new GenomeDistribute(inFile,outFile);
        }
        if(model.equals("testing")){
            String genome = "ATGGATGCATGCAGGGGGCATGCATGCATGGATGCATGGG";
            String cutter1 = "GGG";
            String cutter2 = "ATGG";
            new SubStringSearch(genome,cutter1,cutter2);
        }
        
    }
  
}
