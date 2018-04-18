/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import com.koloboke.collect.map.hash.HashIntObjMap;
import com.koloboke.collect.map.hash.HashIntObjMaps;
import format.table.RowTable;
import htsjdk.samtools.BAMFileReader;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import utils.Benchmark;
import utils.IOFileFormat;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class ThreePrimeExpressionProfiler {
    //The directory of reference genome index
    String referenceGenomeDirS = null;
    //The SampleInformation file (with header), the format is Taxa\tBarcode\tPlateName\tFastqPath
    String sampleInformationFileS = null;
    //The gene annotation file (GTF format)
    String geneAnnotationFileS = null;
    //The path of STAR alignment software
    String starPath = null;
    //The directory of output
    String outputDirS = null;
    //The path of .pgf file
    String gtfPath=null;
    //The path of samtools
    String samtoolsPath=null;
    //The path of geneNewS
    String geneNewS=null;
    //The path of geneName file
    String geneName=null;
    
    int overhangLength = 150;
    
    int multiMapN = 10;
    
    float mismatchRate = (float)0.04;
    
    String[] subDirS = {"subFastqs", "sams", "geneCount","geneTranscript"};
    List<String> fqFileSList = null;
    
    int[] barcodeLengths = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeTaxaMaps = null;
    
    List<String>[] taxaLists = null;
    
    
    public ThreePrimeExpressionProfiler (String parameterFileS) throws IOException {
        this.parseParameters(parameterFileS);
        //mkGeneRange();
//        this.parseFq(); //Remove Ts?remove
        /////////this.mkIndexOfReference(); //one-time set, not requried for standard runs
 //       this.starAlignment(); 
//        this.filterReads();
//        this.cigarPerReads();       
        this.geneCount2(); 
//        this.geneCount();
        this.outputToTable();
//        this.getReadsTrancript();
//        this.transcriptOfGene();
    }
    private void getReadsTrancript(){//use the HashIntIntMap method
        long startTimePoint = System.nanoTime();
        GeneFeature xj = new GeneFeature(this.gtfPath);
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-CIGAR.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        HashIntObjMap[] intObj=ThreePrimeExpressionProfiler.hashIntObjOfExon();
        int numCores = Runtime.getRuntime().availableProcessors();
        nameSet.stream().forEach(p ->{
            String infile = new File (subBamDirS, p+"-CIGAR.txt").getAbsolutePath();
            String outfile=new File(new File(this.outputDirS,subDirS[3]).getAbsolutePath(), p+"-readsTranscript.txt").getAbsolutePath();
            String temp=null;
            String gene=null;
            int startsite=0;
            int endsite=0;
            int chro=0;  
            String index1=null;
            String index2=null;
            List<String> w = new ArrayList<>();
            int[] cont=new int[137117];
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);                             
                while((temp=br.readLine())!=null){
                    chro=Integer.parseInt(temp.split("\t")[0]);
                    if(chro==0){continue;}
                    else{    
                        startsite=Integer.parseInt(temp.split("\t")[1]);
                        endsite=Integer.parseInt(temp.split("\t")[2]);
                        if(intObj[chro-1].get(startsite)==null|intObj[chro-1].get(endsite)==null){}
                        else{
                            index1=intObj[chro-1].get(startsite).toString();                           
                            index2=intObj[chro-1].get(endsite).toString();                        
                            if(index1==index2){ 
                                if(!(w.contains(index1))){
                                    w.add(index1);
                                }
                                int index=w.indexOf(index1);
                                cont[index]++;
                            }
                        }                                                
                    }
                }
                for(int i=0;i<w.size();i++){
                    bw.write(w.get(i)+"\t"+cont[i]);
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) { 
               System.out.println(index1+"\t"+index2);
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Get transcript of each read using HashIntObjMap.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }  
    private void transcriptOfGene() throws IOException{
        long startTimePoint = System.nanoTime();
        String subTranscriptDirS = new File (this.outputDirS,subDirS[3]).getAbsolutePath();
        File[] fs = new File(subTranscriptDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-readsTranscript.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }        
        nameSet.stream().forEach(p ->{
            String infile = new File (subTranscriptDirS, p+"-readsTranscript.txt").getAbsolutePath(); 
            String outfile=new File(new File(this.outputDirS,subDirS[3]).getAbsolutePath(), p+"-transcriptOfGene.txt").getAbsolutePath();                       
            List<String> w = new ArrayList<>();
            List<String>[] transcriptS=new ArrayList[137117];
            for(int i=0;i<transcriptS.length;i++){
                transcriptS[i]=new ArrayList<>();
            }
            String temp=null;
            String transcript=null;
            String gene=null;
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while((temp= br.readLine()) != null){
                    transcript=temp.split("\t")[0];
                    if(transcript.contains("_T")){
                        gene=transcript.split("_T")[0];
                    }else{
                        gene=transcript.split("-T")[0];                        
                    } 
                    if(!(w.contains(gene))){
                            w.add(gene);                           
                    }
                    int index=w.indexOf(gene);
//                    List<String> trans=transcriptS[index];
//                    trans.add(temp);
                    transcriptS[index].add(temp);
                }
                for(int i=0;i<w.size();i++){
                    bw.write(w.get(i)+"\t"+"\t");
                    for(int j=0;j<transcriptS[i].size();j++){
                        bw.write(transcriptS[i].get(j)+"\t");
                    }
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Output the transcript of each gene.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private void geneCount1(){//use the HashIntIntMap method       
        GeneFeature xj=new GeneFeature(this.gtfPath);
 //       GeneFeature xj=new GeneFeature("/Users/xujun/Desktop/Zea_mays.AGPv4.38.pgf");
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-CIGAR.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        HashIntIntMap[] intint=ThreePrimeExpressionProfiler.hashIntIntOfGeneNew();
        xj.sortGeneByStartPosition();
        long startTimePoint = System.nanoTime();
        nameSet.stream().forEach(p ->{
            String infile = new File (subBamDirS, p+"-CIGAR.txt").getAbsolutePath();
            String outfile=new File(new File(this.outputDirS,subDirS[2]).getAbsolutePath(), p+"_1-geneCount.txt").getAbsolutePath();
            String temp=null;
            String name=null;
            int startsite=0;
            int endsite=0;
            int index1=0;
            int index2=0;
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                List w=new ArrayList();
                int count []=new int[20000];
                int chro=0;                
                while((temp=br.readLine())!=null){
                    chro=Integer.parseInt(temp.split("\t")[0]);
                    if(chro==0){continue;}
                    else{    
                        startsite=Integer.parseInt(temp.split("\t")[1]);                       
                        index1=xj.getGeneIndex(chro, startsite);
                        if (index1 < 0) continue;
                        endsite=Integer.parseInt(temp.split("\t")[2]);
                        index2=intint[index1].get(endsite);
                        if(index2<0) continue;
                        if(index1==index2){
                            name=xj.getGeneName(index1);
                            if(!w.contains(name)){
                               w.add(name);                    
                            }
                            int index=w.indexOf(name);
                            count[index]+=1;                                                                               
                        }
                                               
                    }
                }
                bw.write(xj.getGeneName(0)+"\t"+xj.getGeneStart(0)+"\t"+xj.getGeneEnd(0));
                bw.newLine();
                bw.write(w.size());
                bw.newLine();
                for(int i=0;i<w.size();i++){
                   bw.write((String)w.get(i)+"\t"+count[i]);
                   bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) { 
               System.out.println(temp);
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("GeneCount of each sample using HashIntIntap.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }     
     private void geneCount2(){//use the HashIntIntMap method       
        GeneFeature xj=new GeneFeature(this.gtfPath);
 //       GeneFeature xj=new GeneFeature("/Users/xujun/Desktop/Zea_mays.AGPv4.38.pgf");
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-CIGAR.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
//        HashIntIntMap[] intIntOfTrans=ThreePrimeExpressionProfiler.hashIntIntOfTrans();
        HashIntIntMap[] intIntOfGene=ThreePrimeExpressionProfiler.hashIntIntOfGene();
        
        List<String>[] transcriptS=new ArrayList[xj.getTotalTranscriptNumber()];
        for(int i=0;i<transcriptS.length;i++){
            transcriptS[i]=new ArrayList<>();
        }
        long startTimePoint = System.nanoTime();
        nameSet.stream().forEach(p ->{
            String infile = new File (subBamDirS, p+"-CIGAR.txt").getAbsolutePath();
            String outfileG=new File(new File(this.outputDirS,subDirS[2]).getAbsolutePath(), p+"-geneCount.txt").getAbsolutePath();
//            String outfileT=new File(new File(this.outputDirS,subDirS[3]).getAbsolutePath(), p+"-geneTrans.txt").getAbsolutePath();
            String temp=null;
            String name=null;
            String transName=null;
            int startsite=0;
            int endsite=0;
            int index1=0;int index2=0;
//            int index3=0;int index4=0;
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bwG = utils.IOUtils.getTextWriter(outfileG);
//                BufferedWriter bwT = utils.IOUtils.getTextWriter(outfileT);
                List gene=new ArrayList();
                List trans=new ArrayList();
                int countG []=new int[xj.getGeneNumber()];
//                int countT []=new int[transcriptS.length];
                int chro=0;                
                while((temp=br.readLine())!=null){
                    chro=Integer.parseInt(temp.split("\t")[0]);
                    if(chro==0){continue;}
                    else{    
                        startsite=Integer.parseInt(temp.split("\t")[1]);                        
                        index1=intIntOfGene[chro-1].get(startsite);
                        endsite=Integer.parseInt(temp.split("\t")[2]);
                        index2=intIntOfGene[chro-1].get(endsite);                        
                        if(index1==index2){
                            if(index1==0){
                                if(chro==1&startsite>=44289&endsite<=49837){
                                     name=xj.getGeneName(index1);                                    
                                     if(!gene.contains(name)){
                                        gene.add(name);                    
                                     }
                                     int index=gene.indexOf(name);
                                     countG[index]+=1;   
/*                                     index3=intIntOfTrans[index1].get(startsite);
                                     index4=intIntOfTrans[index1].get(endsite);                                    
                                     if(index3==index4){
                                         transName=xj.genes[index1].ts.get(index3).getTranscriptName();
                                         if(!trans.contains(transName)){
                                             trans.add(transName);
                                         }
                                         int indexT=trans.indexOf(transName);
                                         countT[indexT]+=1;
                                     }*/
                                }     
                            }else{
                                name=xj.getGeneName(index1);
                                if(!gene.contains(name)){
                                   gene.add(name);                    
                                }
                                int indexG=gene.indexOf(name);
                                countG[indexG]+=1; 
/*                                index3=intIntOfTrans[index1].get(startsite);
                                index4=intIntOfTrans[index1].get(endsite);                                    
                                if(index3==index4){
                                    transName=xj.genes[index1].ts.get(index3).getTranscriptName();
                                    if(!trans.contains(transName)){
                                        trans.add(transName);
                                    }
                                    int indexT=trans.indexOf(transName);
                                    countT[indexT]+=1;
                                }*/
                            }                                                                                                          
                        }
                                               
                    }
                }
                for(int i=0;i<gene.size();i++){
                   bwG.write((String)gene.get(i)+"\t"+countG[i]);
                   bwG.newLine();
                }
/*                for(int i=0;i<trans.size();i++){
                   bwT.write((String)trans.get(i)+"\t"+countT[i]);
                   bwT.newLine();
                }*/
                br.close();
                bwG.flush();bwG.close();
//                bwT.flush();bwT.close();
            }
            catch (Exception ex) { 
               System.out.println(temp);
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("GeneCount of each sample using HashIntIntap.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private static HashIntIntMap[] hashIntIntOfGene(){//key site value index.establish the HashIntIntMap of the whole chro.
//       String inputFile="/Users/xujun/Desktop/TEP/GeneName.txt";
//        String lengthOfGene="/home/aoyue/xujun/TEP/GeneLengthOfChro.txt";
        String inputFile="/home/aoyue/xujun/TEP/GeneName.txt";   
//        String lengthOfGene="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneLengthOfChro.txt";
        String gtfFile="/home/aoyue/xujun/TEP/Zea_mays.AGPv4.38.modified.gtf";
//        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        HashIntIntMap[] posGeneMaps = new HashIntIntMap[12]; 
        int endindex=0;int startindex=0;
        int length=0;int geneLengthOfChro=0;
        int[][] key = new int[12][];
        int[][] value = new int[12][];
        for(int i=0;i<12;i++){
            startindex=xj.getStartIndexOfChromosome(i+1);
            endindex=xj.getEndIndexOfChromosome(i+1); 
            for(int j=startindex;j<endindex;j++){
                length=xj.getGeneLength(j);
                geneLengthOfChro+=length;
            }
            key[i]=new int[geneLengthOfChro];
            value[i]=new int[geneLengthOfChro];
            geneLengthOfChro=0;
        }
        String map=null;
        String temp=null;        
        int index=0;
        int startsite=0;
        int endsite=0;    
        int chro=0;
        int[] cont=new int[12];
        xj.sortGeneByName();
        try{
            BufferedReader br = utils.IOUtils.getTextReader(inputFile);
            br.readLine();
            while((temp=br.readLine())!=null){               
                index=xj.getGeneIndex(temp);
                chro=xj.getGeneChromosome(index);
                startsite=xj.getGeneStart(index);
                endsite=xj.getGeneEnd(index);
                for(int i=startsite;i<=endsite;i++){
                    key[chro-1][cont[chro-1]]=i;
                    value[chro-1][cont[chro-1]]=index;
                    cont[chro-1]++;
                }
            }            
            br.close();
        }
        catch (Exception ex) { 
            System.out.println(temp+"\t"+chro+"\t"+startsite+"\t"+endsite);
            ex.printStackTrace();              
        }
        for(int i=0;i< 12;i++){
            posGeneMaps[i] = HashIntIntMaps.newImmutableMap(key[i],value[i]);
        }
        return posGeneMaps;        
    }
      private static HashIntIntMap[] hashIntIntOfTrans(){//key site value index.establish the HashIntIntMap of the whole chro.
//       String inputFile="/Users/xujun/Desktop/TEP/GeneName.txt";
//        String lengthOfGene="/home/aoyue/xujun/TEP/GeneLengthOfChro.txt";
        String inputFile="/home/aoyue/xujun/TEP/GeneName.txt";   
//        String lengthOfGene="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneLengthOfChro.txt";
        String gtfFile="/home/aoyue/xujun/TEP/Zea_mays.AGPv4.38.modified.gtf";
//        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        int geneNumber=xj.getGeneNumber();
        HashIntIntMap[] posGeneMaps = new HashIntIntMap[12]; 
        int endindex=0;int startindex=0;
        int length=0;int transLength=0;
        int[][] key = new int[geneNumber][];
        int[][] value = new int[geneNumber][];
        for(int i=0;i<geneNumber;i++){
            int transcriptNumber=xj.genes[i].getTranscriptNumber();
            for(int j=0;j<transcriptNumber;j++){
                length=xj.genes[i].ts.get(j).getTranscriptlength();
                transLength+=length;
            }
            key[i]=new int[transLength];
            value[i]=new int[transLength];
            transLength=0;
        }
        String map=null;
        String temp=null;        
        int index=0;
        int startsite=0;
        int endsite=0;    
        int chro=0;
        int[] cont=new int[geneNumber];
        try{
            BufferedReader br = utils.IOUtils.getTextReader(inputFile);
            br.readLine();
            while((temp=br.readLine())!=null){               
                index=xj.getGeneIndex(temp);
                chro=xj.getGeneChromosome(index);
                int transNumber=xj.genes[index].getTranscriptNumber();
                for(int i=0;i<transNumber;i++){
                    startsite=xj.genes[index].ts.get(i).getTranscriptStart();
                    endsite=xj.genes[index].ts.get(i).getTranscriptEnd();
                    for(int j=startsite;j<endsite;j++ ){
                        key[index][cont[index]]=j;
                        value[index][cont[index]]=i;
                        cont[index]++;
                    }
                }
            }            
            br.close();
        }
        catch (Exception ex) { 
            System.out.println(temp+"\t"+chro+"\t"+startsite+"\t"+endsite);
            ex.printStackTrace();              
        }
        for(int i=0;i< geneNumber;i++){
            posGeneMaps[i] = HashIntIntMaps.newImmutableMap(key[i],value[i]);
        }
        return posGeneMaps;        
    }
      private static HashIntIntMap[] hashIntIntOfGeneNew(){//key site value index.establish the HashIntIntMap of the whole chro.
        String inputFile="/Users/xujun/Desktop/TEP/GeneName.txt";
//        String lengthOfGene="/home/aoyue/xujun/TEP/GeneLengthOfChro.txt";
//        String inputFile="/home/aoyue/xujun/TEP/GeneName.txt";   
//        String lengthOfGene="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneLengthOfChro.txt";
//        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
//          String gtfFile="/home/aoyue/xujun/TEP/Zea_mays.AGPv4.38.modified.gtf";
        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        
 //       xj.sortGeneByName();
        int mapNumber=xj.getGeneNumber();
        HashIntIntMap[] posGeneMaps = new HashIntIntMap[mapNumber]; 
        int size=0;
        int[][] key = new int[mapNumber][size];
        int[][] value = new int[mapNumber][size];
        String temp=null;
        int index=0;
        int startsite=0;
        int endsite=0;
        int[] cont=new int[mapNumber];
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            br.readLine();            
            while((temp=br.readLine())!=null){                             
                index=xj.getGeneIndex(temp);
                startsite=xj.getGeneStart(index);
                endsite=xj.getGeneEnd(index); 
                size=xj.getGeneLength(index);
                key[index]=new int[size];
                value[index]=new int[size];
                for(int i=startsite;i<endsite+1;i++){
                    key[index][cont[index]]=i;
                    value[index][cont[index]]=index;
                    cont[index]++;
                }
            }
        }
        catch (Exception ex) { 
            System.out.println(temp+"\t"+index+"\t"+size+"\t"+cont[index]);
            ex.printStackTrace();              
        }
        for(int i=0;i< mapNumber;i++){
            posGeneMaps[i] = HashIntIntMaps.newMutableMap(key[i],value[i]);
        }
        return posGeneMaps;        
    }
    private HashIntObjMap hashIndexGeneMap(){
        String inputFile="/Users/xujun/Desktop/GeneNewS.txt";
        String temp=null;
        int size =45805;
        int[] key=new int[size];
        String [] value=new String[size];
        int cont=0;
        try{
            BufferedReader br = utils.IOUtils.getTextReader(inputFile);
            while((temp=br.readLine())!=null){
                key[cont]=Integer.valueOf(temp.split("\t")[1]);
                value[cont]=temp.split("\t")[2];
            }
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }
        HashIntObjMap indexGeneMap = HashIntObjMaps.newImmutableMap(key, value);
        return indexGeneMap;
    }
    private static HashIntObjMap[] hashIntObjOfExon(){//key site value index.establish the HashIntIntMap of the whole chro.
//        String inputFile="/home/aoyue/xujun/TEP/ExonInfor.txt";
//        String lengthOfGene="/home/aoyue/xujun/TEP/ExonLengthOfChro.txt";
        String inputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/ExonInfor.txt";
        String lengthOfGene="/Users/xujun/Desktop/RNA_seq/twice/STAR/ExonLengthOfChro.txt";
        HashIntObjMap[] posExonMaps = new HashIntObjMap[12];   
        int[][] key = new int[12][];
        String[][] value = new String[12][];
        String map=null;
        String temp=null;
        int chro=0;
        String transcript=null;
        int index=0;
        String positionNewS=null;
        int startsite=0;
        int endsite=0;
        int[] cont=new int[12];
//        List<String> w = new ArrayList<>();
        try{
            BufferedReader br = utils.IOUtils.getTextReader(inputFile);
            BufferedReader br1 = utils.IOUtils.getTextReader(lengthOfGene);
            HashMap<Integer,Integer> lengthChroMap = new HashMap<>();
            while ((map=br1.readLine())!=null) {
                lengthChroMap.put( Integer.parseInt(map.split("\t")[0]),Integer.parseInt(map.split("\t")[1]));
            }
            for(int i=0;i<12;i++){
                key[i]=new int[lengthChroMap.get(i+1)];
                value[i]=new String[lengthChroMap.get(i+1)];
            }
            while((temp=br.readLine())!=null){                     
                    chro=Integer.parseInt(temp.split("\t")[1]);
                    positionNewS=temp.split("\t")[2];
                    transcript=temp.split("\t")[0];
//                    w.add(transcript);
//                    index=w.indexOf(transcript);
                    for(int j=0;j<positionNewS.split(";").length;j++){
                       startsite=Integer.parseInt(positionNewS.split(";")[j].split(":")[0]);
                       endsite=Integer.parseInt(positionNewS.split(";")[j].split(":")[1]);    
                       for (int i = startsite; i < endsite; i++) {
                         key[chro-1][cont[chro-1]] = i;
                         value[chro-1][cont[chro-1]] = transcript;
                         cont[chro-1]++;
                       } 
                    }                                                                  
            }            
            br.close();
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }
        for(int i=0;i< 12;i++){
            posExonMaps[i] = HashIntObjMaps.newMutableMap(key[i],value[i]);
        }
        return posExonMaps;        
    }
    private void outputToTable() throws IOException{
        long startTimePoint = System.nanoTime();
        RowTable<String> t = new RowTable<>(this.geneName);
        HashMap<String, Integer> geneIndexMap = new HashMap<>();
        List<String> geneList = new ArrayList<>();
        List<String> w = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            geneIndexMap.put(t.getCell(i, 0), i);
            geneList.add(t.getCell(i, 0));
            w.add("0");
        }
        String subBamDirS = new File (this.outputDirS,subDirS[2]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-geneCount.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        int numCores = Runtime.getRuntime().availableProcessors();
        nameSet.stream().forEach(p ->{
            String infile = new File (subBamDirS, p+"-geneCount.txt").getAbsolutePath(); 
            String outfile=new File(new File(this.outputDirS,subDirS[2]).getAbsolutePath(), "CountTable.txt").getAbsolutePath();
            
            String temp=null;
            String name=null;
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while((temp= br.readLine()) != null){
                    name=temp.split("\t")[0];
                    int index=geneIndexMap.get(temp.split("\t")[0]);
                    w.set(index, temp.split("\t")[1]);                                                   
                }
                t.addColumn(p, w);
                t.writeTextTable(outfile, IOFileFormat.Text);
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Output as a table.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private void geneCount(){
        long startTimePoint = System.nanoTime();
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-CIGAR.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        GeneFeature pgf=new GeneFeature(this.gtfPath);
        nameSet.parallelStream().forEach(p -> {
            String infile = new File (subBamDirS, p+"-CIGAR.txt").getAbsolutePath();
            String outfile=new File(new File(this.outputDirS,subDirS[2]).getAbsolutePath(), p+"-geneCount.txt").getAbsolutePath();
            BufferedReader br = utils.IOUtils.getTextReader(infile);
            BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
            String a=null;
            int index=0;
            String name=null;
            String yon=null;
            List w=new ArrayList();
            int count []=new int[20000];
            try{
                while((a=br.readLine())!=null){
                    index=pgf.getGeneIndex(Integer.parseInt(a.split("\t")[0]),Integer.parseInt(a.split("\t")[1]));
                    if (index < 0) continue;
                    name=pgf.getGeneName(index);                                                          
                    if(!w.contains(name)){
                        w.add(name);                    
                    }               
                    int position=w.indexOf(name);
                    count[position]+=1;                   
                }
                for(int i=0;i<w.size();i++){
                    bw.write((String)w.get(i)+"\t"+count[i]);
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            catch(Exception ex){
                System.out.println(a);
                ex.printStackTrace();
            } 
        });
        StringBuilder time = new StringBuilder();
        time.append("Gene count of each sample.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
     }
    private void cigarPerReads(){
        long startTimePoint = System.nanoTime();
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-news.txt");
        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        int numCores = Runtime.getRuntime().availableProcessors();
        nameSet.parallelStream().forEach(p -> {
            String infile = new File (subBamDirS, p+"-news.txt").getAbsolutePath();
            String outfile=new File(subBamDirS, p+"-CIGAR.txt").getAbsolutePath();
            String a=null;
            String seq=null;
            String q=null;
            String seq1=null;                   
            String q1=null;
            String startpos=null;
            String startpos1=null;
            String chro=null;
            String cigar=null;
            try{
                BufferedReader br = utils.IOUtils.getTextReader(infile);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while ((a= br.readLine()) != null){
                    chro=a.split("\t")[0];
                    startpos=a.split("\t")[1];                   
                    cigar=a.split("\t")[2];
                    seq=a.split("\t")[3];
                    q=a.split("\t")[4];
                    int cont=0;
                    int cont2=0;
                    int cont3=0;
                    int cont4=0;
                    for(int i=0;i<cigar.length();i++){                           
                            if (!Character.isDigit(cigar.charAt(i))){ 
                                if(cigar.charAt(i)=='M'){
                                    startpos1=String.valueOf(Integer.parseInt(startpos)+cont3);
                                    seq1=seq.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));  
                                    q1=q.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));
                                    bw.write(chro+"\t");bw.write(startpos1+"\t");bw.write(Integer.parseInt(startpos1)+seq1.length()+"\t");bw.write(seq1+"\t");bw.write(q1+"\t");
                                    bw.newLine();
                                    cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                    cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                    cont=i+1;
                                }else{
                                    if(cigar.charAt(i)=='D'|cigar.charAt(i)=='N'){
                                        cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                    }else{ 
                                        if(cigar.charAt(i)=='I'){}
                                        else{
                                            cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        }
                                        cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                    }
                                }                                
                            }else{
                                continue;
                            }   
                }                              
            }
            br.close();
            bw.flush();bw.close();
                
            }
            catch (Exception ex) {
               System.out.println(p+"\t"+a+"\t1234");  
               ex.printStackTrace();
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Parse each read with CIGAR.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private void filterReads(){
        long startTimePoint = System.nanoTime();
        String subBamDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
        File[] fs = new File(subBamDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.sortedByCoord.out.bam");
        List<File> fList = Arrays.asList(fs);
        int numCores = Runtime.getRuntime().availableProcessors();
        fList.stream().forEach(f -> {            
            StringBuilder sb = new StringBuilder();           
            sb.append(this.samtoolsPath).append(" view -F 0x100 ").append(f).append(" >> ").append(f.getName().replaceFirst("Aligned.sortedByCoord.out.bam", ".sam"));
            String command = sb.toString();
            System.out.println(command);
            try {
                BufferedWriter bw1 = utils.IOUtils.getTextWriter(new File(subBamDirS,f.getName().replaceFirst("Aligned.sortedByCoord.out.bam", ".sam")).getAbsolutePath());
                File dir = new File(subBamDirS);
                String []cmdarry ={"/bin/bash","-c", command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                BufferedReader br1 = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                while ((temp = br1.readLine()) != null) {
                    System.out.println(temp);
                }                
                p.waitFor();
                String command2="cut -f 3,4,6,10,11 "+new File(subBamDirS,f.getName().replaceFirst("Aligned.sortedByCoord.out.bam", ".sam")).getAbsolutePath()+" >> "+f.getName().replaceFirst("Aligned.sortedByCoord.out.bam", "-news.txt");//+"/"+f.getName().replaceFirst("Aligned.sortedByCoord.out.bam", ".sam"
                System.out.println(command2);
                String[] cmdarry2={"/bin/bash","-c", command2};
                Process pb=Runtime.getRuntime().exec(cmdarry2,null,dir);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
        StringBuilder time = new StringBuilder();
        time.append("Get the uniquely mapped and multimapped number < 10 reads and get the usefull column.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    private void starAlignment () throws IOException { //private int[][] getIndividualGeneCount ()
        long startTimePoint = System.nanoTime();
        String subFqDirS = new File (this.outputDirS, subDirS[0]).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
            BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
            br.readLine();br.readLine();br.readLine();br.readLine();
            if(br.readLine()!=null){
                File file=fs[i];
                fList.add(file);
            }
        }
        int numCores = Runtime.getRuntime().availableProcessors();
        //int[][] geneCount = new int[taxaNumber][];
        fList.stream().forEach(f -> {
            //int[] indiGeneCount = new int[geneNumber];
            StringBuilder sb = new StringBuilder();
            sb.append(this.starPath).append(" --runThreadN ").append(numCores).append(" --genomeDir ").append(this.referenceGenomeDirS);
            sb.append(" --sjdbGTFfile ").append(this.geneAnnotationFileS).append(" --sjdbOverhang ").append(this.overhangLength);
            sb.append(" --readFilesIn ").append(f);
            sb.append(" --outFileNamePrefix ").append(new File(new File(this.outputDirS, subDirS[1]).getAbsolutePath(), f.getName().replaceFirst(".fq", ""))
                .getAbsolutePath()).append(" --outFilterMultimapNmax ").append(this.multiMapN);
            sb.append(" --outFilterMismatchNoverLmax ").append(this.mismatchRate)
                .append(" --outFilterIntronMotifs RemoveNoncanonicalUnannotated ");
            //sb.append(" --outSAMtype SAM");
            sb.append(" --outSAMtype BAM SortedByCoordinate");
            String command = sb.toString();
            System.out.println(command);
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                p.waitFor();
 //               command = ">";// pipe to sam samtools view > sam
 //               p = rt.exec(command);
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            //int index = Collections.binarySearch(taxaList, taxonName);
            //geneCount[index] = indiGeneCount;
            System.out.println("Finished"+f);
        });
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
        
    }
    
    private void mkIndexOfReference () {
        int numCores = Runtime.getRuntime().availableProcessors();
        String referenceGenomeFileS = "/Users/feilu/Documents/database/maize/reference/AGPv4/maizeAGPv4.fa";
        String outputDirS = "/Users/feilu/Documents/database/maize/reference/starLib";
        try {
            StringBuilder sb = new StringBuilder("/Users/feilu/Software/STAR-2.5.4b/bin/MacOSX_x86_64/STAR");
            sb.append(" --runThreadN ").append(numCores).append(" --runMode genomeGenerate --genomeDir ").append(outputDirS);
            sb.append(" -- genomeFastaFiles ").append(referenceGenomeFileS);
            String command = sb.toString();
            System.out.println(command);
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            
        }
        catch (Exception ee) {
            ee.printStackTrace();
        }
    }
    
    private void parseFq () {
        long startTimePoint = System.nanoTime();
        fqFileSList.parallelStream().forEach(f -> {
            int fqIndex = Collections.binarySearch(this.fqFileSList, f);
            String subFqDirS = new File (this.outputDirS, subDirS[0]).getAbsolutePath();
            List<String> barcodeList = barcodeLists[fqIndex];
            String[] subFqFileS = new String[barcodeList.size()];
            HashMap<String, String> btMap = barcodeTaxaMaps[fqIndex];
            Set<String> barcodeSet = btMap.keySet();
            BufferedWriter[] bws = new BufferedWriter[subFqFileS.length];
            HashMap<String, BufferedWriter> barcodeWriterMap = new HashMap<>();
            for (int i = 0; i < subFqFileS.length; i++) {
                String taxon = btMap.get(barcodeList.get(i));
                subFqFileS[i] = new File(subFqDirS, taxon+".fq").getAbsolutePath();
                bws[i] = IOUtils.getTextWriter(subFqFileS[i]);
                barcodeWriterMap.put(barcodeList.get(i), bws[i]);
            }
            int barcodeLength = this.barcodeLengths[fqIndex];
            try {
                BufferedReader br = null;
                if (f.endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f);
                }
                else {
                    br = IOUtils.getTextReader(f);
                }
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                BufferedWriter tw = null;
                int cnt = 0;
                int cnt2 = 0;
                while((temp = br.readLine())!=null){
                    cnt2++;
                    seq = br.readLine();
                    currentBarcode = seq.substring(0, barcodeLength);
                    int cutIndex = 0;
                    if (barcodeSet.contains(currentBarcode)) {
                        tw = barcodeWriterMap.get(currentBarcode);
                        tw.write(temp);
                        tw.newLine();
                        tw.write(seq.substring(8));
                        tw.newLine();
                        tw.write(br.readLine());
                        tw.newLine();
                        tw.write(br.readLine().substring(8));
                        tw.newLine();
/*                        int cont=0;
                        for(int i= barcodeLength;i<seq.length();i++){
                            if(i==seq.length()-1){
                                br.readLine();br.readLine();
                                continue;
                            }
                            if(seq.charAt(i)=='T'){}
                            else{
                                cont++;
                                if(cont>3){
                                    cutIndex = i;
                                    tw.write(temp);
                                    tw.newLine();
                                    tw.write(seq.substring(cutIndex-4));
                                    tw.newLine();
                                    tw.write(br.readLine());
                                    tw.newLine();
                                    tw.write(br.readLine().substring(cutIndex-4));
                                    tw.newLine();
                                    cnt++;
                                    if (cnt%1000000 == 0) System.out.println("Wrote " + String.valueOf(cnt) + " sequences. " + f);
                                    break;
                                }
                            }
                            
                        }*/
//                        byte[] seqB = seq.getBytes();
//                        for(int i=barcodeLength; i<seqB.length; i+=4){
//                            if (i+4 > seq.length()) {
//                                br.readLine();br.readLine();
//                                continue;
//                            }
//                            int value = this.getsum(seqB, i, i+4);
//                            //317 = TTAA
//                            if (value < 317) {
//                                cutIndex = i;
//                                tw.write(temp);
//                                tw.newLine();
//                                tw.write(seq.substring(cutIndex));
//                                tw.newLine();
//                                tw.write(br.readLine());
//                                tw.newLine();
//                                tw.write(br.readLine().substring(cutIndex));
//                                tw.newLine();
//                                cnt++;
//                                if (cnt%1000000 == 0) System.out.println("Wrote " + String.valueOf(cnt) + " sequences. " + f);
//                                break;
//                            }              
//                        }
                    }
                    else {
                        br.readLine();br.readLine();
                        continue;
                    }
                }
                StringBuilder sb = new StringBuilder();
                sb.append(cnt).append(" out of ").append(cnt2).append(", ").append(((float)(double)cnt/cnt2)).append(" of total reads were parsed from " + f);
                System.out.println(sb.toString());
                for (int i = 0; i < subFqFileS.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    
    private int getsum(byte[] seqB, int startIndex, int endIndex){
        int sum=0;
        for(int i = startIndex; i < endIndex; i++){           
            sum+=seqB[i];
        }
        return sum;
    }
    
    private void parseParameters(String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("Three' Expression Profiler (TEP)")) ifOut = true;
            if (!(temp = br.readLine()).equals("Author: Jun Xu, Fei Lu")) ifOut = true;
            if (!(temp = br.readLine()).equals("Email: liuzhongxujun@163.com; flu@genetics.ac.cn")) ifOut = true;
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) ifOut = true;
            if (ifOut) {
                System.out.println("Thanks for using TEP.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.referenceGenomeDirS = pLineList.get(0);
        this.sampleInformationFileS = pLineList.get(1);
        this.geneAnnotationFileS = pLineList.get(2);
        this.starPath = pLineList.get(3);
        this.outputDirS = pLineList.get(4);
        this.gtfPath=pLineList.get(5);
        this.samtoolsPath=pLineList.get(6);
        this.geneNewS=pLineList.get(7);
        this.geneName=pLineList.get(8);
        this.processTaxaAndBarcode();
        for (int i = 0; i < this.subDirS.length; i++) {
            new File(this.outputDirS, subDirS[i]).mkdir();
        }
    }
    
    private void processTaxaAndBarcode () {
        RowTable<String> t = new RowTable<>(this.sampleInformationFileS);
        Set<String> fqSet = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            fqSet.add(t.getCell(i, 3));
        }
        fqFileSList = new ArrayList<>(fqSet);
        Collections.sort(fqFileSList);
        barcodeLengths = new int[fqFileSList.size()];
        barcodeLists = new ArrayList[fqFileSList.size()];
        taxaLists = new ArrayList[fqFileSList.size()];
        barcodeTaxaMaps = new HashMap[fqFileSList.size()];
        int[] cnts = new int[fqFileSList.size()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            cnts[index]++;
        }
        for (int i = 0; i < cnts.length; i++) {
            barcodeLists[i] = new ArrayList<>();
            taxaLists[i] = new ArrayList<>();
            barcodeTaxaMaps[i] = new HashMap<>();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
            cnts[index]++;
            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);
            taxaLists[index].add(taxon);
            barcodeLists[index].add(t.getCell(i, 1));
            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
            barcodeLengths[index] = t.getCell(i, 1).length();
        }
    }
    
    public static void main(String args[]) throws IOException {
        new ThreePrimeExpressionProfiler(args[0]);
    }
    
}
