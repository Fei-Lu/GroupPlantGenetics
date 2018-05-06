/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

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
     * -cutter: two REs, consist with -model ReducedLib
     */
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();
        int len = args.length;
        String model = "test";
        String inFile = "",keyFileName = "",inFile2 = null;
        String outFile="",names = "";
        int ExSize = 5*1000;
        String RE1="GGATCC",RE2="CCGG";
        String bed=null;
        boolean gzip = false,hmp = false,byChr=false;
        boolean plink = false, blink = false,depth = false;
        String size = "0";
        for (int i = 0; i < len; i++){
            if (null != args[i])switch (args[i]) {
                case "--model":
                    model = args[i+1];
                    i++;
                    break;
                case "--i":
                    inFile = args[i+1];
                    i++;
                    break;
                case "--o":
                    outFile = args[i+1];
                    i++;
                    break;
                case "--size":
                    double ExSize1 = Double.valueOf(args[i+1]);
                    ExSize = (int) (ExSize1 * 1000);
                    size = args[i+1];
                    i++;
                    break;
                case "--cutter":
                    RE1 = args[i+1];
                    RE2 = args[i+2];
                    i = i+2;
                    break;
                case "--make-plink":
                    plink = true;
                    break;
                case "--make-blink":
                    blink = true;
                    break;
                case "--file":
                    inFile = args[i+1];
                    i++;
                    break;
                case "--out":
                    outFile = args[i+1];
                    i++;
                    break;
                case "--inFile2":
                    inFile2 = args[i+1];
                    i++;
                    break;
                case "--gzip":
                    gzip = true;
                    break;
                case "--hmp":
                    //for model is GH
                    hmp = true;
                    break;
                case "--byChr":
                    //for model is GH
                    byChr = true;
                    break;
                case "--subname":
                    //for model is GH
                    names = args[i+1];
                    break;
                case "--key":
                    keyFileName = args[i+1];
                    i++;
                    break;  
                case "--bed":
                    bed = args[i+1];
                    i++;
                    break; 
                case "--depth":
                    depth = true;
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
            //
            new LogReadStatistic(inFile);
        }
        if(model.equals("ExFastq")){
            new ExtractionOfFastq(inFile,outFile,ExSize);
        }
        if(model.equals("GenoStat")){
            if(!size.equals("0")){
                new GenomeStatistic(inFile,outFile,size);
            }else{
                new GenomeStatistic(inFile,outFile);
            }
        }
        if(model.equals("GenoDist")){
            System.out.println("GenoDist is working...");
            if(byChr){
                new GenomeDistribute(inFile,outFile,byChr);
            }else{
                if(size.equals("0")){
                    new GenomeDistribute(inFile,outFile);
                }else{
                    new GenomeDistribute(inFile,outFile,size);
                }
            }
        }
        if(model.equals("test")){
           new SNPsStat(inFile,outFile); 
        }
        if(model.equals("ReducedLib")){
            //Generate a new fasta format file, only contain sequence between 2 REs;
            new ReducedLibrary(inFile,RE1,RE2,outFile);
        }
        
        if (model.equals("format")){
            new FormatTransform(inFile,plink,blink,outFile,gzip);
        }
        if (model.equals("RemoveDup")){
            new RemoveDuplication(inFile,outFile);
        }
        if (model.equals("GH")){
            new GenerateHybrid(inFile,outFile,hmp);
        }
        if (model.equals("gff3")){
            new gff3ToGeneList(inFile, outFile);
        }
        if (model.equals("vcf")){
            if(depth){
                new VcfTools(inFile,outFile,true);
            }
//            new VcfTools(inFile, outFile);
        }
        if (model.equals("V3ToGene")){
            new maizeV3toGenelist(inFile,outFile);
        }
        if (model.equals("rmpro")){
            new RemoveDupProtein(inFile,outFile);
        }
        if (model.equals("subFasta")){
            new SubstractFromFasta(inFile,names,outFile);
        }
        if (model.equals("orderTree")){
            new newickPattern(inFile);
        }
        if (model.equals("GBS")){
            new GBStassel(inFile,keyFileName);
        }
        if (model.equals("changeChr")){
            new changeChromosome(inFile,outFile,bed);
        }
        if (model.equals("merge")){
            new MergeGenome(inFile,outFile);
        }
        if (model.equals("ToLower")){
            new toLowerBase(inFile,inFile2,outFile);
        }
        if (model.equals("psl")){
            new modifyPsl(inFile);
        }
        if (model.equals("unmask")){
            new unmask(inFile);
        }
        long endTime = System.currentTimeMillis();
        int timeLast = (int) ((endTime-startTime)/1000);
        System.out.println("Process finished in: "+ timeLast + " seconds");
    }
  
}
