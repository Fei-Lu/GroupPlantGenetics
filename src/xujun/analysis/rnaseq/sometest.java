/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

import pgl.infra.utils.IOUtils;

/**
 *
 * @author xujun
 */
public class sometest {
    public sometest(){
//        this.test1();
//        this.test2();
//            this.getExonInfor();
//        this.getExonLengthOfChro();
//        this.getReadsTransscript();
//            this.getGeneName();
//        this.transfer();
//        this.sortByName();
//        this.mappedNumber();
        this.find();
    }
    public void find(){
        String inputfile="";
        String inputDirS="";
        RowTable<String> t = new RowTable<>(inputfile);
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            nameList.add(t.getCell(i, 0));
        }        
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> sampleSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            sampleSet.add(fs[i].getName().split(".html")[0]);
        }
        for(int i=0;i<nameList.size();i++){
            if(!sampleSet.contains(nameList)){
                System.out.println(nameList.get(i));
            }
        }
        
    }
    public void mappedNumber(){
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/sams/TEPWithoutB.sam";
        String outputFile="/Users/xujun/Desktop/TEP/TEPOut/sams/morethan10.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            String temp=null;
            int cont =0;
            while((temp= br.readLine()) != null){
                if(Integer.parseInt(temp.split("\t")[11].split(":")[2])>10){
                    bw.write(temp);bw.newLine();
                    cont++;
                }                
            }
            System.out.print(cont);
            bw.write(cont);
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
        }
    }
    public void sortByName(){
        ArrayList<String> al = new ArrayList<String>();
        al.add("Zm00001d027231_T001");
        al.add("Zm00001d027231_T002");
        al.add("Zm00001d027231_T031");
        al.add("Zm00001d027231_T004");
        al.add("Zm00001d027231_T005");
        Collections.sort(al);
        System.out.println("List after the use of" +
                           " Collection.sort() :\n" + al);
    }
    public void transfer(){
        String inputFile="/Users/xujun/Desktop/TEP/GeneName.txt";
        String outputFile="/Users/xujun/Desktop/TEP/GeneName1.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            String temp=null;
            String geneName=null;
            int chro=0;
            while((temp= br.readLine()) != null){
                if(temp.contains("zma-mir")){
                    temp=temp.replaceAll("zma-mir", "zma-MIR");
                }
                bw.write(temp);bw.newLine();
            }
            
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
        } 
    }
    public void getGeneName(){
        String inputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        String outputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneName.txt";
        HashSet<String> geneNameSet = new HashSet();
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            String temp=null;
            String geneName=null;
            int chro=0;
            while((temp= br.readLine()) != null){
                geneName=temp.split("\t")[8].split(";")[1].split(":")[1].substring(0,temp.split("\t")[8].split(";")[1].split(":")[1].length()-1);
                chro=Integer.valueOf(temp.split("\t")[0]);
                if(!(geneNameSet.contains(geneName))){
                    geneNameSet.add(geneName);
                    bw.write(geneName+"\t"+chro);
                    bw.newLine();
                }                
            }
            String[] GeneName = geneNameSet.toArray(new String[geneNameSet.size()]);
            
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
        } 
    }
    public void getReadsTransscript(){
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/sams/TEP_CIGAR.txt";
        String outputFile="/Users/xujun/Desktop/TEP/TEPOut/sams/TEP_Exon.txt";
        String pgfFile="/Users/xujun/Desktop/Zea_mays.AGPv4.38.pgf";
        GeneFeature xj = new GeneFeature(pgfFile);
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            String temp=null;
            String name=null;
            int chro=0;
            int startsite=0;
            int endsite=0;
            while((temp= br.readLine()) != null){
                chro=Integer.parseInt(temp.split("\t")[0]);
                startsite=Integer.parseInt(temp.split("\t")[1]);
                endsite=Integer.parseInt(temp.split("\t")[1]);
                xj.getTranscriptName(startsite, endsite);
                bw.write(name+"\t"+chro);
                bw.newLine();
            }
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
        } 
    }
    public void getExonInfor(){//输出的文件是 transcript chro position
        String gtfFile="/home/aoyue/xujun/Exon/Zea_mays.AGPv4.38.modified.gtf";
        String outputFile="/home/aoyue/xujun/Exon/ExonInfor.txt";
//        String output="/Users/xujun/Desktop/RNA_seq/twice/STAR/ExonLengthInfor.txt";
        List<String> w=new ArrayList();
        int index=0;
        String[] position=new String[400000];
        for(int i=0;i<position.length;i++){
            position[i]="";
        }
        try{
            BufferedReader br = IOUtils.getTextReader(gtfFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
 //           BufferedWriter bw1 = utils.IOUtils.getTextWriter(output);
            String temp=null;
            String genotype=null;
            String transcriptnews=null;
            String transcript=null;
            int chro=0;
 //           int length=0;
 //           int [] chroLength=new int[12];
            while((temp= br.readLine()) != null){        
                genotype=temp.split("\t")[2];
                if(genotype.equals("exon")){
                    chro=Integer.parseInt(temp.split("\t")[0]);
                    transcriptnews=temp.split("\t")[8];
                    transcript=transcriptnews.split(" ")[1].substring(12,transcriptnews.split(" ")[1].length()-2);
                    if(w.contains(transcript)){
                        index=w.indexOf(transcript);
                        position[index]=position[index].concat(temp.split("\t")[3]+":"+temp.split("\t")[4]+";");
                    }
                    else{
                        w.add(transcript);
                        index=w.indexOf(transcript);
                        position[index]=position[index].concat(chro+"\t"+temp.split("\t")[3]+":"+temp.split("\t")[4]+";");
                    }
                    
                }
            }
            
            for(int i=0;i<w.size();i++){
                bw.write(w.get(i)+"\t"+position[i]);
                bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
         }  
    }
    public void getExonLengthOfChro(){//输出的文件是 transcript chro position
        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        String outputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/ExonLengthOfChro.txt";
//        String output="/Users/xujun/Desktop/RNA_seq/twice/STAR/ExonLengthInfor.txt";
        int index=0;
        try{
            BufferedReader br = IOUtils.getTextReader(gtfFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
 //           BufferedWriter bw1 = utils.IOUtils.getTextWriter(output);
            String temp=null;
            String genotype=null;
            int chro=0;
           int length=0;
           int [] chroLength=new int[12];
            while((temp= br.readLine()) != null){        
                genotype=temp.split("\t")[2];
                if(genotype.equals("exon")){
                    chro=Integer.parseInt(temp.split("\t")[0]);
                  length=Integer.parseInt(temp.split("\t")[4])-Integer.parseInt(temp.split("\t")[3])+1;
                  chroLength[chro-1]+=length;
                }
            }
            for(int i=0;i<12;i++){
                bw.write(i+1+"\t"+chroLength[i]);
                bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception ex){
             ex.printStackTrace();
         }  
    }
    public void test1(){//Zm00001d031168 对这个超高表达的基因进行统计 第一列是read的长度 第二列是质量值 第三列是这条read的质量均值
        String inputfile="/Users/xujun/Desktop/notSIDname.text";
        String outputfile="/Users/xujun/Desktop/Zm00001d031168.text"; 
        String news=null;
        String phred=null;
        
        try{
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                phred=news.split("\t")[2];
                int phred1=0;
                if(news.split("\t")[0].equals("Zm00001d031168")){
                    for(int i=0;i<phred.length();i++){
                        phred1=phred1+(int)phred.charAt(i)-33;
                    }
                    bw.write(news.split("\t")[1]/*+"\t"+phred*/+"\t"+phred1/phred.length()); 
                    bw.newLine();
                }
                
             }
             br1.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }           
        
    }
    public void test2(){//找出大于10000的gene
        String inputfile="/Users/xujun/Desktop/notSDIgenecount.text";
        String outputfile="/Users/xujun/Desktop/high.text";
        String news=null;
        String number=null;
        String length=null;
        String name=null;
        String phred=null;
        try{
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                number=news.split("\t")[1];               
                if(Integer.parseInt(number)==1){
                    name=news.split("\t")[0];
                    length=news.split("\t")[2];
                    phred=news.split("\t")[3];
                    if(length.length()>=3&Integer.parseInt(phred)>25){
                        bw.write(name+"\t"+length+"\t"+phred);             
                        bw.newLine();
                    }
                    
                } 
                    
                }
             br1.close();
             bw.flush();bw.close();                     
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }           
        
    }
    public void test3(){//
        String inputfile="/Users/xujun/Desktop/notSDIgenecount.text";
        String outputfile="/Users/xujun/Desktop/high.text";
        String news=null;
        String number=null;
        try{
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                number=news.split("\t")[1]; 
                for(int i=1;i<=10;i++){
                    if(1<=Integer.parseInt(number)&Integer.parseInt(number)<100){
                        
                    }
                }
                if(Integer.parseInt(number)==1){
                                 
                    bw.newLine();
                } 
                    
                }
             br1.close();
             bw.flush();bw.close();                     
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         } 
    }

    public static void main(String[] args){
        new sometest();
    }
}
