/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.Range;
import pgl.infra.range.RangeValStr;
import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import rcaller.RCaller;
import rcaller.RCode;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import xuebo.analysis.annotation.FStringUtils;
import static xujun.analysis.rnaseq.DistinguishSample.getKeyList;

/**
 *
 * @author xujun
 */
public class WheatRNASeq20181107 {
    int[] barcodeLengths = null;
    List<String> indexList = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeMethodMaps = null;
//    HashMap<Integer,Integer> 
    HashMap barcodeStrain = new HashMap();
    HashMap indexStrain = new HashMap();
    List<String>[] methodLists = null;
    
    List<String> allMethodList = new ArrayList<String>();
    String inputDirS=null;
    String outputDirS=null;
    int geneNumber=92;
    
     public WheatRNASeq20181107(){//String parameterFileS
//         this.parseParameters(parameterFileS);
//        this.processMethodAndBarcode();
//        this.parseFq();
//        this.parseFqByBarcodePlateseq();//两端文件都读入 但是只分出了不带ployT的那一端！
//            this.checkFq();
//        this.ployTLengthAndQ();
//        this.averageQ();
//        this.HTSeqCountPair();
//        this.parseFqByIndexAndBarcode();
//        this.dataValumn();
//        this.HTSeqCountSingle();
//        this.HTSeqCountSingle();
//        this.HTSeqCountMerge();
//        this.HTSeqCountMergeRPKM();
//        this.parseFqByIndex();
//        this.findIndex();
//        this.removeUnmappedStartEnd();
//        this.avergeTranscriptLengthInWheat();
//        this.splitSample();//双端文件读入，R1和R2都输出，用于双端比对
//        this.findTriad();
//        this.geneSequence();
//        this.shScoreParse();
//        this.getList();
//        this.highDiverganceGene();
//        this.isSame();
//        this.geneNumber();
//        this.overlapSNP();
//        this.NPos();
//        this.changPos();
        
    }
     public void NPos(){
         BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/Axiom_wheat660.na34.annot.35N35.fa").getAbsolutePath());
         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/pos.txt").getAbsolutePath());
         try{
             String tem =null;String temp=null; 
             while((temp=br.readLine())!=null){
                 temp=br.readLine();
                 for(int i=0;i<temp.length();i++){
                     if(temp.charAt(i)=='N'){
                         System.out.println(i);
                         bw.write(i);
                         bw.newLine();
                     }
                 }

             }
             br.close();bw.flush();bw.close();
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
     }
     public void changPos(){
         BufferedReader br =IOUtils.getTextReader(new File("/data1/home/junxu/tree/chr1D.filt.vcf").getAbsolutePath());
         BufferedWriter bw =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/changechr1D.vcf").getAbsolutePath());
//         BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/AA358.unique.mapped.ok.hmp.txt").getAbsolutePath());
//         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/changePosAA358.unique.mapped.ok.hmp.txt").getAbsolutePath());
         try{
             String [] tem =null;String temp=null; 
             while((temp=br.readLine())!=null){
                if(!temp.startsWith("#")){
                    tem=temp.split("\t");
                    for(int i=0;i<tem.length;i++){
                        if(i!=0){
                            bw.write(tem[i]+"\t");
                        }else{
                            bw.write("15"+"\t");
                        }
                    }
                    bw.newLine();
                }else{
                    bw.write(temp);bw.newLine();
                }  
             }
             br.close();bw.flush();bw.close();
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
         
     }
     public void overlapSNP(){
         HashIntIntMap posline = null;
         posline = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
         ArrayList infor = new ArrayList();List<Integer> h1= new ArrayList();
         BufferedReader br1 =IOUtils.getTextReader(new File("/data1/home/junxu/tree/changechr1B.vcf").getAbsolutePath());
         BufferedReader br2 =IOUtils.getTextReader(new File("/data1/home/junxu/tree/zhengchr1B.vcf").getAbsolutePath());
         BufferedWriter bw1 =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/overjiaochr1B.vcf").getAbsolutePath());
         BufferedWriter bw2 =IOUtils.getTextWriter(new File("/data1/home/junxu/tree/overzhengchr1B.vcf").getAbsolutePath());
//         BufferedReader br1 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/zhengchr1A.vcf").getAbsolutePath());
//         BufferedReader br2 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/tree/overzhengchr1A.vcf").getAbsolutePath());
//         BufferedWriter bw1 =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/overjiaochr1A.vcf").getAbsolutePath());
//         BufferedWriter bw2 =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/tree/overzhengchr1A.vcf").getAbsolutePath());
         
         try{
//            for (int i=0;i<rt.getRowNumber();i++){
//                if(rt.getCell(i,0).equals("1B")){
//                    pos.add(rt.getCell(i, 1));
//                    posline.put(Integer.valueOf(rt.getCell(i, 1).toString()).intValue(),i);
//                }
//            }
            String [] temt =null;String tempt=null;int i=0;
            while((tempt=br1.readLine())!=null){
                temt=tempt.split("\t");
                if(temt[0].equals("8")){
                    i++;
                    infor.add(tempt);
                    posline.put(i,Integer.valueOf(temt[1].toString()).intValue());
                }
            }
            String [] tem =null;String temp=null; 
            while((temp=br2.readLine())!=null){
                tem=temp.split("\t");
                if(!tem[0].startsWith("#")){                    
                    if(posline.containsValue(Integer.valueOf(tem[1]).intValue())){
                        h1=getKeyList(posline,Integer.valueOf(tem[1]).intValue());
                        for(int j=0;j<h1.size();j++){
                            bw1.write(infor.get(h1.get(j).intValue()-1).toString());
                            bw1.newLine();
                        }
                        bw2.write(temp);bw2.newLine();
                    }
                    
                }
                    
            }
            bw1.flush();bw1.close();br2.close();bw2.flush();bw2.close();
//            for(int i=1;i<geneNumber.length;i=i+2){
//                for(int j=0;j<geneNumber[i].length;j++){
//                    System.out.println(i+"\t"+geneNumber[i][j]);
//                }
//            }
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
         
     }
     public static List<Integer> getKeyList(HashIntIntMap map,int value){//the method to get key from value
         List<Integer> keyList = new ArrayList();
         for(int getKey: map.keySet()){
             if(map.get(getKey)==value){
                 keyList.add(getKey);
             }
         }
         return keyList;
     }
     public void geneNumber(){
         RowTable rt = new RowTable("/Users/xujun/Desktop/tempory/readme.txt");
         HashMap subChrLength=new HashMap();
         for( int i=0;i< rt.getRowNumber();i++){
             subChrLength.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 2));
         }
         BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3").getAbsolutePath());
         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneDis.txt").getAbsolutePath());
         int [][] geneNumber = new int [45][10];
//         GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
         int erro=0;
         try{
            String [] tem =null;String temp=null; 
            while((temp=br.readLine())!=null){
                tem=temp.split("\t");
                if(tem[2].startsWith("gene")){
                    if(Integer.valueOf(tem[0])%2!=0 && Integer.valueOf(tem[0])!=0){
                        int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])+1).toString());
                        double pos = (double)Integer.valueOf(tem[3])/chrLength;
                        bw.write(pos+"");bw.newLine();
                    }else{
                        if(Integer.valueOf(tem[0])!=0){
                            int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString());
                            double pos =(double) (Integer.valueOf(tem[3])+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString()))/chrLength;
//                            geneNumber[Integer.valueOf(tem[0])-1][pos]++;
                            bw.write(pos+"");bw.newLine();
                        }
                    }
                }
                    
            }
            bw.flush();bw.close();br.close();
//            for(int i=1;i<geneNumber.length;i=i+2){
//                for(int j=0;j<geneNumber[i].length;j++){
//                    System.out.println(i+"\t"+geneNumber[i][j]);
//                }
//            }
         }
         catch (Exception e) {
             System.out.println(erro);
            e.printStackTrace();
         } 
     }
     
     public void isSame(){
         BufferedReader br1 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAB.txt").getAbsolutePath());
         BufferedReader br2 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAD.txt").getAbsolutePath());
         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/samA.txt").getAbsolutePath());
         List<String> nameList1 = new ArrayList();
         List<String> nameList2 = new ArrayList();
         try{
            String [] tem =null;String temp=null;
            while((temp=br1.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList1.add(tem[1].split("\t")[1]);
                }
            }
            while((temp=br2.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList2.add(tem[1].split("\t")[1]);
                }
            }
            for(int i=0;i<nameList2.size();i++){
                if(nameList1.contains(nameList2.get(i))){
                    bw.write(nameList2.get(i));bw.newLine();
                }
            }
            br1.close();br2.close();
            bw.flush();bw.close();
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
     }
     public void highDiverganceGene(){
         RowTable rt = new RowTable("/Users/xujun/Desktop/tempory/readme.txt");
         HashMap chrStart=new HashMap();
         HashMap chrEnd = new HashMap();
         GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
         int geneNumber =gf.getGeneNumber();
         for( int i=0;i< rt.getRowNumber();i++){
             chrStart.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 4));
             chrEnd.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 5));
         }
         String temp=null;String position []= null;
         BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAB.txt").getAbsolutePath());
         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/AB.txt").getAbsolutePath());
         try{
            String [] tem =null;
            while((temp=br.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    if(Integer.valueOf(tem[2])%2==0){
                        int pos = Integer.valueOf(tem[3])+Integer.valueOf(chrStart.get(Integer.valueOf(tem[2])).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[2])).toString());
//                        System.out.println(pro);
                        bw.write(pro+"\t");
                    }else{
                        int pos=Integer.valueOf(tem[3]);
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[2])+1).toString());
//                        System.out.println(pro);
                        bw.write(pro+"\t");
                    }
                    if(Integer.valueOf(tem[6])%2==0){
                        int pos = Integer.valueOf(tem[7])+Integer.valueOf(chrStart.get(Integer.valueOf(tem[6])).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])).toString());
                        bw.write(pro+"");bw.newLine();
                    }else{
                        int pos=Integer.valueOf(tem[7]);
//                        int end = Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])+1).toString());
                        double pro=(double)pos/Integer.valueOf(chrEnd.get(Integer.valueOf(tem[6])+1).toString());
                        bw.write(pro+" ");bw.newLine();
                    }
//                    bw.write(tem[2]+"\t"+tem[6]);bw.newLine();
                }
            }
            bw.flush();bw.close();
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
     }
     public void getList(){//to get the data from the distance matrix
         ColumnTable ct = new ColumnTable("/Users/xujun/Desktop/Tempory/g.txt");
         BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneticList.txt").getAbsolutePath());
         try{
            for(int i =0 ; i < ct.getColumnNumber() ; i++){
                for(int j =i ; j < ct.getRowNumber();j++){
                   bw.write(ct.getCellAsString(i, j));bw.newLine();     
                }
            }
            bw.flush();bw.close();
         }
         catch (Exception e) {
            e.printStackTrace();
         } 
     }
     public void shScoreParse(){
        String inputFileS = "/Users/xujun/summarys.txt";
        String outputFileS = "/Users/xujun/shScore.txt";
        String temp=null;String shtemp=null;
        double sh=0;double shM=0;int k=0;int tempK=0;String num=null;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            while((temp=br.readLine())!=null){
                if(temp.startsWith("logical")){
                    bw.write(String.valueOf(k));bw.newLine();
                    sh=0;shM=0;k=0;tempK=0;
                }
                if(temp.startsWith("Silhouette")){
                    tempK=Integer.valueOf(temp.split(" ")[5]);
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();br.readLine();
                    shtemp=br.readLine();
                    if(!shtemp.contains("Individual")){
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }else{
                        br.readLine();shtemp=br.readLine();
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }
                     
                    if(sh>shM){
                        shM=sh;k=tempK;
                        if(k==6){
                            
                        }
                    }
                }
            }
            br.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            System.out.println(sh);
            e.printStackTrace();
        }
     }
     public void shScore (){
        String inputfsK="/Users/xujun/Desktop/eQTL/total/Homology/shScore";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".txt");
        for(int i=0;i<fsK.length;i++){
            if(!fsK[i].getName().contains("temp")){
                fListK.add(fsK[i]);
            }
        }
        RCaller callerL = new RCaller();
        callerL.setRscriptExecutable("/usr/local/bin/R");
        RCode rCodeL = new RCode();
        rCodeL.addRCode("library(cluster);library(factoextra);library(vegan);install.packages(\"Runiversal\");");
        callerL.setRCode(rCodeL);
        callerL.runOnly();
        List KPool =new ArrayList();
        fListK.stream().forEach(f -> {
            double sh=0;double shM=0;int k=0;
            for(int i=1;i<11;i++){
                RCaller caller = new RCaller();
                caller.setRscriptExecutable("//usr//local//bin//R");
                RCode rCode = new RCode();
                rCode.addRCode("a<-read.csv(\""+f+"\",sep=\"\\t\",header = F)");
                rCode.addRCode("results=kmeans(a,");
                rCode.addRCode(i+");dis = dist(a^2);");
                rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
                caller.setRCode(rCode);
                caller.runAndReturnResult("sh");
                sh=caller.getParser().getAsDoubleArray("sh")[0];
                if(sh>shM){
                    shM=sh;
                }else{
                    k=i;
                    break;
                }
            }
            KPool.add(k);
        });
        for(int i=0;i<KPool.size();i++){
            System.out.println(KPool.get(i));
        }
     }


     public void countTPM(){
        RCaller caller = new RCaller();
        caller.setRscriptExecutable("/usr/local/bin/R");
        RCode rCode = new RCode();
        rCode.addRCode("library(cluster);library(factoextra);library(vegan);results=kmeans(");
        rCode.addRCode(""+",");
        for(int i=1;i<11;i++){
            rCode.addRCode(i+");dis = dist(a");
            rCode.addRCode(""+")^2;");
        }
        rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
        caller.setRCode(rCode);
        caller.runAndReturnResult("sh");
    }
     public void geneSequence(){
//         String inputFile="/Users/xujun/Desktop/IGVmaterial/wheat/abd_iwgscV1.fa";
         SiPASHomo ss = new SiPASHomo();
     }
     public void findTriad(){
         String inputFileS="/Users/xujun/Desktop/eQTL/total/DESeq/TaS-D/750MSiPAS-TPM.txt";
         RowTable rt = new RowTable(inputFileS);
         HashMap AExpre = new HashMap();
         HashMap BExpre = new HashMap();
         HashMap DExpre = new HashMap();
         for (int i=0;i<rt.getRowNumber();i++){
             if(rt.getCellAsDouble(i, 1)>0.3 ){
                 if(rt.getCellAsString(i, 0).charAt(8)=='A'){
                     AExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                 }else{
                     if(rt.getCellAsString(i, 0).charAt(8)=='B'){
                         BExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                     }else{
                         DExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                     }
                 }               
             }
         }
        Iterator iterA = AExpre.entrySet().iterator(); 
        while (iterA.hasNext()) { 
            HashMap.Entry entryA = (HashMap.Entry)iterA.next(); 
            Object keyA = entryA.getKey();
            Object valA = entryA.getValue(); 
            Iterator iterB = BExpre.entrySet().iterator(); 
            while(iterB.hasNext()){
                HashMap.Entry entryB = (HashMap.Entry)iterB.next(); 
                Object keyB = entryB.getKey();
                if(keyB.toString().charAt(7)==keyA.toString().charAt(7)){
                    Object valB = entryB.getValue(); 
                    if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valB.toString()))){
                        Iterator iterD = DExpre.entrySet().iterator();
                        HashMap.Entry entryD = (HashMap.Entry)iterD.next();    
                        Object keyD = entryD.getKey();
                        if(keyA.toString().charAt(7)==keyD.toString().charAt(7)){
                            Object valD = entryD.getValue();
                            while(iterD.hasNext()){
                                if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valD.toString()))){
                                    System.out.println(keyA.toString()+"\n"+keyB.toString()+"\n"+keyD.toString());
                                }
                            }
                        }              
                    }
                } 
            }
        } 
     }
     public boolean isEqual(double a,double b){
        double eps = 1E-1;
        if (Math.abs(a - b) / a <= eps) {
            return true;
        }else{
            return false;
        }
     }
     public void splitSample(){
         String inputDirS="/Users/xujun/Desktop/Undetermined_K19Pool1030_L008";
         String outputDirS ="/Users/xujun/Desktop/Undetermined_K19Pool1030_L008/split";
         RowTable rt = new RowTable("/Users/xujun/Desktop/RNA-seq.txt");
         HashMap barcodeName = new HashMap();
         List<String> nameList = new ArrayList();
         for (int i=0;i<rt.getRowNumber();i++){
             barcodeName.put(rt.getCellAsString(i, 0).substring(1), rt.getCellAsString(i, 1));
             nameList.add(rt.getCellAsString(i, 1));
         }
         File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("L008_")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"L008_R1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"L008_R2.fq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[nameList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
                for (int i = 0; i < nameList.size(); i++) {
                   bw[i]=IOUtils.getTextWriter(new File(outputDirS, nameList.get(i)+"_R1.fq").getAbsolutePath());
                   bw1[i]=IOUtils.getTextWriter(new File(outputDirS, nameList.get(i)+"_R2.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile2);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                String temp = null;String seq=null;
                String currentBarcode=null;
                while ((temp = br.readLine()) != null) {
                    seq=br.readLine();
                    if(barcodeName.get(seq.substring(1, 8))!=null){
                        int index = nameList.indexOf(barcodeName.get(seq.substring(1, 8)));
                        bw[index].write(temp);bw[index].newLine();
                        bw[index].write(seq);bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                    }else{
                        br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                }
                br.close();br1.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                    bw1[i].flush();bw1[i].close();
                }
                
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
     }
     public void avergeTranscriptLengthInWheat(){
         String inputFile="/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3";
         GeneFeature gf = new GeneFeature(inputFile);
         long totalTranscriptLength = 0;
         int transcriptNumber = 0;
         for(int i=0;i<gf.genes.length;i++){
             for(int j=0;j<gf.genes[i].ts.size();j++){
                 totalTranscriptLength+=gf.genes[i].ts.get(j).getTranscriptEnd()-gf.genes[i].ts.get(j).getTranscriptStart()+1;
                 transcriptNumber++;
             }
         }
         System.out.println(totalTranscriptLength/transcriptNumber);
     }
     private void removeUnmappedStartEnd () {
        List<String> diffValumnMethodList = new ArrayList<String>();
        String inputDirS = new File("/data1/home/junxu/analysis0215/test-result/test-withoutERCC").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace(".sam", "")); 
        }
        Collections.sort(diffValumnMethodList);
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(inputDirS+"/change"+f.getName());
                String temp = null;
                String [] tem = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("[")) continue;
                    if(temp.startsWith("@")){
                        bw.write(temp);bw.newLine();
                    }else{
                         List<String> tList= FStringUtils.fastSplit(temp);
                         tem = tList.toArray(new String[tList.size()]);
                         if(!tem[2].equals("*")){
                            bw.write(temp);bw.newLine();
                         }
                    }   
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void findIndex(){
        String inputFile="/data1/home/junxu/wheat/splitByIndex/indexSiPAS-all.fq";
        List <String> indexList=new ArrayList<>();       
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextWriter(inputFile.replace(".fq", ".txt"));
            String temp=null;
            while((temp = br.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                if(!indexList.contains(index)){
                    indexList.add(index);
                }
            }
            BufferedReader br1 = IOUtils.getTextReader(inputFile);
            int [] indexNumber = new int[indexList.size()];
            while((temp = br1.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                indexNumber[indexList.indexOf(index)]++; 
            }
            for(int i=0;i<indexList.size();i++){
                bw.write(indexList.get(i)+"\t"+indexNumber[i]);
                bw.newLine();
            }
            br.close();br1.close();
            bw.flush();bw.close();
        }
        catch (Exception ex) { 
            ex.printStackTrace();
        }
    }
    public void HTSeqCountMergeRPKM(){
        List <String> nameList=new ArrayList<>();
        List<String> fileList=new ArrayList<>();
//        String subFqDirS = new File (this.outputDirS).getAbsolutePath();
        String subFqDirS = new File ("/Users/xujun/Desktop/TEP/TEPOut/HTSeqLibrary").getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();   
        fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
        List<File> fList = Arrays.asList(fs);
        for(int i=0;i < fList.size();i++){
            fileList.add(fList.get(i).getName().replace("Count.txt", ""));
        }
        Collections.sort(fileList);
        int[][] geneCount = new int[this.geneNumber][]; // chr, geneByChr, taxon
        double [][] RPKM=new double[this.geneNumber][];
        int [] allReadsNumber = new int [fList.size()];
        for (int i = 0; i < geneNumber; i++) {
            geneCount[i] = new int[fList.size()];   
            RPKM[i] = new double[fList.size()]; 
        }
        fList.stream().forEach(f -> {
            String temp=null;String[] tem = null;
            try{           
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                while((temp = br.readLine()) != null){
                    List<String> tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    if(tem[0].startsWith("ERCC")){
//                    if(tem[0].startsWith("TraesCS")){
                        if(!nameList.contains(tem[0])){
                            nameList.add(tem[0]);
                        }
                        int index=nameList.indexOf(tem[0]);
                        geneCount[index][fileList.indexOf(f.getName().replace("Count.txt", ""))]=Integer.parseInt(tem[1]);
                        allReadsNumber[fileList.indexOf(f.getName().replace("Count.txt", ""))]+=Integer.parseInt(tem[1]);
                    }      
                }

            }
            catch (Exception ex) {
                System.out.println(tem[0]+"\t1234");  
                ex.printStackTrace();

            }
        });
        String row = null;String [] tem =null;
        HashMap<String,String> geneNameAndLength = new HashMap();
        try{           
//                BufferedReader br = IOUtils.getTextReader("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/ERCC92.gtf");
                BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/wheat/ERCC92/ERCC92.gtf");
                while((row = br.readLine()) != null){
                    List<String> tList= PStringUtils.fastSplit(row);
                    tem = tList.toArray(new String[tList.size()]);
                    geneNameAndLength.put(tem[0], tem[4]);
                }

        }
        catch (Exception ex) { 
            ex.printStackTrace();
        }
        for(int i=0;i<fList.size();i++){
            for(int j=0;j<this.geneNumber;j++){
                RPKM[j][i]=geneCount[j][i]*1000000000/allReadsNumber[i]/Integer.parseInt(geneNameAndLength.get(nameList.get(j)));
            }
        }
        String outputFileS = new File (this.outputDirS,"RowAndRPKM-countResult.txt").getAbsolutePath();
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            sb.append("Gene"+"\t");
            for(int i=0;i<fileList.size();i++){            
                sb.append(fileList.get(i).replace("Count.txt", "")+"\t"+"RPKM"+"\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for(int i=0;i<geneCount.length;i++){
                sb = new StringBuilder();  
                for(int j=0;j<fileList.size();j++){
                    if(j==0){
                        sb.append(nameList.get(i)+"\t");
                    }
                    sb.append(geneCount[i][j]+"\t"+RPKM[i][j]+"\t");           
                }
                bw.write(sb.toString());
                bw.newLine();
            }

            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
 public void HTSeqCountMerge(){//顺序是HTSeq里面的顺序

    List <String> nameList=new ArrayList<>();
    List<String> fileList=new ArrayList<>();
    String subFqDirS = new File (this.outputDirS).getAbsolutePath();
    File[] fs = new File(subFqDirS).listFiles();   
    fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
    List<File> fList = Arrays.asList(fs);
    for(int i=0;i < fList.size();i++){
        fileList.add(fList.get(i).getName().replace("Count.txt", ""));
    }
    Collections.sort(fileList);
//    int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
//    double [][][] RPKM=new double[chrs.length][][];
//    int [][] count=new int[92][fList.size()];
    int [][] count=new int[110790][fList.size()];
    fList.stream().forEach(f -> {
        String temp=null;String[] tem = null;
        try{           
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            while((temp = br.readLine()) != null){
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
//                if(tem[0].startsWith("ERCC")){
                  if(tem[0].startsWith("TraesCS")){
                    if(!nameList.contains(tem[0])){
                        nameList.add(tem[0]);
                    }
                    int index=nameList.indexOf(tem[0]);
                    count[index][fileList.indexOf(f.getName().replace("Count.txt", ""))]=Integer.parseInt(tem[1]);
                }      
            }
            
        }
        catch (Exception ex) {
            System.out.println(tem[0]+"\t1234");  
            ex.printStackTrace();

        }
    });
    
    String outputFileS = new File (this.outputDirS,"countResult.txt").getAbsolutePath();
    try{
        StringBuilder sb = new StringBuilder();
        BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
        sb.append("Gene"+"\t");
        for(int i=0;i<fileList.size();i++){            
            sb.append(fileList.get(i).replace("Count.txt", "")+"\t");
        }
        bw.write(sb.toString());
        bw.newLine();
        for(int i=0;i<count.length;i++){
            sb = new StringBuilder();  
            for(int j=0;j<fileList.size();j++){
                if(j==0){
                    sb.append(nameList.get(i)+"\t");
                }
                sb.append(count[i][j]+"\t");           
            }
            bw.write(sb.toString());
            bw.newLine();
        }
        
        bw.flush();
        bw.close();
    }
    catch (Exception e) {
        e.printStackTrace();
    }
    
}
     public void HTSeqCountSingle(){
        String inputDirS = this.inputDirS;
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {  
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s no ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat/rightchangewheat.gtf").append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File (this.outputDirS).getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
     public void HTSeqCountPair(){
        String inputDirS = "/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/Out/diffValumnSams-D";
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {  
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s reverse ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/ERCC92.gtf").append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File ("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/Out/HTSeqCount-D").getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
    public void dataValumn() {
        String subFqDirS = new File ("/data1/home/junxu/wheat/R2").getAbsolutePath();
//        String subFqDirS = "xujun/TEP/TEPOut/subFastqs";
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, "R2.fq");
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fList.add(fs[i]);
        }
//        String command =null;
        int numCores = Runtime.getRuntime().availableProcessors();
        fList.stream().forEach(f -> {
                try { 
                    File dir = new File(new File ("/data1/home/junxu/wheat/singleR2All/subfq").getAbsolutePath());
                    for( int i = 1 ;i<=12;i++){
                        StringBuilder sb = new StringBuilder();
                        sb.append("/data1/home/junxu/wheat/seqtk/seqtk sample -s100 ").append(f+" ").append(1000000*i+" ").append(" >> ").append(" "+1000000*i+f.getName().replace(".fq", "_R2.fq"));
                        String command = sb.toString();
                        System.out.println(command);  
                        String []cmdarry ={"/bin/bash","-c",command};
                        Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                        p.waitFor();
                    } 
                }
                catch (Exception e) {
                    e.printStackTrace();
                }            
        });
    }
    public void ployTLengthAndQ () {
//        String barcodeFileS = "/data1/home/junxu/wheat/RNA-seq20181030.txt";
        String inputDirS ="/data1/home/junxu/wheat/doubleAll/subFastqs";     
//        String barcodeFileS = "/Users/xujun/Desktop/RNA-seq20181030.txt";
        String outputDirS ="/data1/home/junxu/wheat/doubleAll/lengthPloyTAndQ";      
//        RowTable<String> rt = new RowTable<>(barcodeFileS);
//        int rowNumber = rt.getRowNumber();
//        int columnNumber = rt.getColumnNumber();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R2.fq");
        fs = IOUtils.listFilesStartsWith(fs, "5000000");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(p -> {
            String seq = null;
            String temp=null;
            String phred=null;
            try {
                BufferedReader br = IOUtils.getTextReader(p.getAbsolutePath());                              
                BufferedWriter bwL = IOUtils.getTextWriter(new File(outputDirS,p.getName().replace(".fq", "L.txt")).getAbsolutePath());
                BufferedWriter bwQ = IOUtils.getTextWriter(new File(outputDirS,p.getName().replace(".fq", "Q.txt")).getAbsolutePath());
                List L=new ArrayList();
                List Q=new ArrayList();
                while ((temp = br.readLine()) != null){
                    int a=0;int allPhred = 0;                 
                    seq=br.readLine(); br.readLine();
                    phred=br.readLine();
                    for(int i=8;i<seq.length();i++){
                        if(seq.charAt(i)!='T'){
                            a++; 
                        }
                        if(a>2){
                            L.add(i-8);
                            for(int j=i-2;j<phred.length();j++){
                                allPhred+=(int)phred.charAt(j)-33;
                            }
                            Q.add(allPhred/(150-i));
                            break;
                        }   
                    }
                    if(a<=2){
                        L.add(0);
                        for(int j=8;j<phred.length();j++){
                            allPhred+=(int)phred.charAt(j)-33;
                        }
                        Q.add(allPhred/142);
                    } 
                }
                Collections.shuffle(L); 
                Collections.shuffle(Q);
                int randomSeriesLength = 10000;
                List<Integer> randomSeriesL = L.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriesQ = Q.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      bwL.write(randomSeriesL.get(k)+ " ");bwL.newLine();
                      bwQ.write(randomSeriesQ.get(k)+ " ");bwQ.newLine();
                   } else {
                       bwL.write(randomSeriesL.get(k)+ " ");     
                       bwQ.write(randomSeriesQ.get(k)+ " ");
                   }
                }
                br.close();
                bwL.flush();bwL.close();
                bwQ.flush();bwQ.close();
            }               
            catch (Exception ex) {
                   System.out.println(seq+"\t1234");  
                   ex.printStackTrace();

            }
        });      
     }
    public void averageQ () {
        String outputDirS ="/data1/home/junxu/wheat/doubleAll/lengthPloyTAndQ";
        String inputDirS ="/data1/home/junxu/wheat/doubleAll/subFastqs";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R2.fq");
        fs = IOUtils.listFilesStartsWith(fs, "6000000");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(p -> {
            String seq = null;
            String temp=null;
            String phred=null;
            try {
                BufferedReader br = IOUtils.getTextReader(p.getAbsolutePath());                              
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS,p.getName().replace(".fq", "Q.txt")).getAbsolutePath());
                List Q=new ArrayList();
                while ((temp = br.readLine()) != null){
                    int allPhred=0;br.readLine();br.readLine();phred=br.readLine();
                    for(int j=0;j<phred.length();j++){
                        allPhred+=(int)phred.charAt(j)-33;
                    }
                    Q.add(allPhred/150);
                }
                Collections.shuffle(Q); 
                int randomSeriesLength = 10000;
                List<Integer> randomSeries = Q.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      bw.write(randomSeries.get(k)+ " ");bw.newLine();
                   } else {
                      bw.write(randomSeries.get(k)+ " ");                     
                   }
                }
                br.close();
                bw.flush();bw.close();
            }               
            catch (Exception ex) {
                   System.out.println(seq+"\t1234");  
                   ex.printStackTrace();

            }
        });      
     }
    public void checkFq(){
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".fastq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {  
            try{
                BufferedReader br = null;
                if (f.getAbsolutePath().endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                System.out.println(f);
                String temp = null;
                int count =0;
                while((temp = br.readLine())!=null){
                    count++; 
                    if(temp.startsWith("@")){
                        br.readLine();br.readLine();br.readLine();
                    }
                    else{
                        System.out.println(temp+"\t"+count);  
                        break;
                    }    
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    private void parseFq () {
        long startTimePoint = System.nanoTime();
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/splitByIndex";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {            
            try {
                BufferedWriter[][] bws = null;
                bws = new BufferedWriter[indexList.size()][];
                BufferedWriter [] bw = new BufferedWriter[indexList.size()];
                for (int i = 0; i < indexList.size(); i++) {
                    bws[i]=new BufferedWriter[barcodeLists[i].size()];
                    for(int j=0;j<barcodeLists[i].size();j++){
                       String fileName=barcodeMethodMaps[i].get(barcodeLists[i].get(j)); 
                       bws[i][j]=IOUtils.getTextWriter(new File(outputDirS, fileName+".fq").getAbsolutePath());
                    }
                    if(barcodeLists[i].size()==0){
                        String fileName=barcodeMethodMaps[i].get(indexList.get(i));
                        bw[i]=IOUtils.getTextWriter(new File(outputDirS, fileName+".fq").getAbsolutePath());
                    }
                }
                BufferedReader br = null;
                if (f.getAbsolutePath().endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                String index1 = null ;String index2 = null ;
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    index1=temp.split(" ")[1].split(":")[3].substring(1, 6);
                    if(indexList.contains(index1)){
                        int pos1 = indexList.indexOf(index1);
                        seq = br.readLine();
                        if(barcodeLists[pos1].size()!=0){
                            currentBarcode = seq.substring(1,8);
                            if(barcodeLists[pos1].contains(currentBarcode)){
                                int pos2 = barcodeLists[pos1].indexOf(currentBarcode);
                                bws[pos1][pos2].write(temp);bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(seq);bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(br.readLine());bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(br.readLine());bws[pos1][pos2].newLine();
                            }else{
                                br.readLine();br.readLine();
                            }
                            
                        }else{
                            bw[pos1].write(temp);bw[pos1].newLine();
                            bw[pos1].write(seq);bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                        }   
                    }else{
                        index2=temp.split(" ")[1].split(":")[3].substring(1);                       
                        if(indexList.contains(index2)){
                            seq=br.readLine();
                            int pos1 = indexList.indexOf(index2);
                            bw[pos1].write(temp);bw[pos1].newLine();
                            bw[pos1].write(seq);bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                        }else{
                            br.readLine();br.readLine();br.readLine();
                        }
                    }     
                }
                br.close();
                for(int i=0;i<indexList.size();i++){
                    for(int j=0;j<barcodeLists[i].size();j++){
                        bws[i][j].flush();bws[i][j].close();
                    }
                    if(barcodeLists[i].size()==0){
                        bw[i].flush();bw[i].close();
                    }    
                }
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
    private void parseFqByIndex () {
        long startTimePoint = System.nanoTime();
        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq.txt");
//       RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030.txt");
        List<String> indexList = new ArrayList() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
           barcodeStrain.put(t.getCell(i, 2), t.getCell(i, 0));
           indexList.add(t.getCell(i, 2));
        }
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/Novoseq";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, "R2_001.fq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {      
            try {
                BufferedWriter [] bw = new BufferedWriter[indexList.size()];
                for (int i = 0; i < indexList.size(); i++) {
//                   bw[i]=IOUtils.getTextWriter(new File(outputDirS, f.getName()+"_"+barcodeMethod.get(barcodeList.get(i))+".fq").getAbsolutePath());
                    bw[i]=IOUtils.getTextWriter(new File(outputDirS,barcodeStrain.get(indexList.get(i))+".fq").getAbsolutePath());
                }
                BufferedReader br = null;
//                BufferedReader br1 = null;
                if (f.getAbsolutePath().endsWith(".fq.gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                int index1 = -1 ;String index2 = null ;
                String temp = null;String temp1 = null;
                String seq = null;String seq1 = null;
                String currentIndex = null;
                while((temp = br.readLine())!=null){
                    seq=br.readLine();//seq1=br1.readLine();temp1=br1.readLine();
                    currentIndex=temp.split(":")[9];
                    if(barcodeStrain.get(currentIndex)!=null){
                        index1=indexList.indexOf(currentIndex);
                        bw[index1].write(temp);bw[index1].newLine();
                        bw[index1].write(seq);bw[index1].newLine();
                        bw[index1].write(br.readLine());bw[index1].newLine();
                        bw[index1].write(br.readLine());bw[index1].newLine();
                    }else{
                        br.readLine();br.readLine();
                    }
                }          
                for(int i=0;i<indexList.size();i++){
                      bw[i].flush();bw[i].close();
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
    private void parseFqByIndexAndBarcode () {
        long startTimePoint = System.nanoTime();
        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq20181030.txt");
//       RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030.txt");
        List<String> indexList = new ArrayList() ;
        List<String> barcodeList = new ArrayList() ;
//        for (int i = 0; i < t.getRowNumber(); i++) {
//           barcodeMethod.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
//           barcodeList.add(t.getCell(i, 1).substring(1));
//        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0).split("-")[0]);
            barcodeStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            barcodeList.add(t.getCell(i, 1).substring(1));
        }
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/doubleAll";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, "_001.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {      
            try {
                String infile1 = new File (inputDirS, p+"-R1_001.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"-R2_001.fq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[barcodeList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[barcodeList.size()];
                for (int i = 0; i < barcodeList.size(); i++) {
                   bw[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R1.fq").getAbsolutePath());
                   bw1[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R2.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                int pos = -1 ;
                String temp = null;
                String seq = null;String index = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    index=temp.split(":")[9].substring(1,6);
                    if(indexStrain.get(index)!=null){
                        seq=br.readLine();
                        currentBarcode=seq.substring(1,8);
                        if(barcodeStrain.get(currentBarcode)!=null){
                            pos=barcodeList.indexOf(currentBarcode);
                            bw[pos].write(temp);bw[pos].newLine();
                            bw[pos].write(seq);bw[pos].newLine();
                            bw[pos].write(br.readLine());bw[pos].newLine();
                            bw[pos].write(br.readLine());bw[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                        }else{
                            br.readLine();br.readLine();
                            br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                        }
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                    
                }          
                for(int i=0;i<barcodeList.size();i++){
                      bw[i].flush();bw[i].close();
                      bw1[i].flush();bw1[i].close();
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
    private void parseFqByBarcodePlateseq () {
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq20181030.txt");
        RowTable<String> t = new RowTable<>("/data1/home/junxu/analysis0319/3/3RNA-seq.txt");
//        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC.txt");
//        RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030.txt");
        List<String> barcodeList = new ArrayList() ;
        List<String> strainList = new ArrayList() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
           barcodeStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 1));
           barcodeList.add(t.getCell(i, 2).substring(1));
           strainList.add(t.getCell(i, 1));
        }
        int [] count =new int [strainList.size()];
        String inputDirS = "/data1/home/junxu/analysis0319/total/3";
        String outputDirS = "/data1/home/junxu/analysis0319/total/3-double";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {      
            try {
                String infile1 = new File (inputDirS, p+"_R1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_R2.fq.gz").getAbsolutePath();
                BufferedWriter [] bwR1 = new BufferedWriter[barcodeList.size()];
                BufferedWriter [] bwR2 = new BufferedWriter[barcodeList.size()];
                BufferedWriter bw1 = IOUtils.getTextWriter(new File (outputDirS,"ratio.txt").getAbsolutePath());
                BufferedWriter bw2 = IOUtils.getTextWriter(new File (outputDirS,"error.txt").getAbsolutePath());
                for (int i = 0; i < barcodeList.size(); i++) {
                   bwR1[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R1.fq").getAbsolutePath());
                   bwR2[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R2.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                
                int index1 = -1 ;int index2 = -1 ;
                String temp = null;String temp1 = null;
                String seq = null;String seq1 = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    seq=br.readLine();temp1=br1.readLine();seq1=br1.readLine();
                    currentBarcode=seq.substring(1,8);
                    if(barcodeStrain.get(currentBarcode)!=null){
                        index1=barcodeList.indexOf(currentBarcode);
                        index2=strainList.indexOf(barcodeStrain.get(currentBarcode));
                        count[index2]++;
//                        if(String.valueOf(barcodeMethod.get(currentBarcode)).contains("SiPAS")){
                            bwR1[index1].write(temp);bwR1[index1].newLine();
                            bwR1[index1].write(seq);bwR1[index1].newLine();
                            bwR1[index1].write(br.readLine());bwR1[index1].newLine();
                            bwR1[index1].write(br.readLine());bwR1[index1].newLine();
                            bwR2[index1].write(temp1);bwR2[index1].newLine();
                            bwR2[index1].write(seq1);bwR2[index1].newLine();
                            bwR2[index1].write(br1.readLine());bwR2[index1].newLine();
                            bwR2[index1].write(br1.readLine());bwR2[index1].newLine();
                            
//                        }else{
//                            br1.readLine();br1.readLine();
//                        }
//                        br.readLine();br.readLine();
//                        br1.readLine();br1.readLine();
//                        bw[index1].write(temp);bw[index1].newLine();
//                        bw[index1].write(seq);bw[index1].newLine();
//                        bw[index1].write(br.readLine());bw[index1].newLine();
//                        bw[index1].write(br.readLine());bw[index1].newLine();
                    }else{
                        bw2.write(temp);bw2.newLine();
                        bw2.write(seq);bw2.newLine();
                        bw2.write(br.readLine());bw2.newLine();
                        bw2.write(br.readLine());bw2.newLine();
                        br1.readLine();br1.readLine();
                    }
                }  
                for(int i=0;i<count.length;i++){
                    bw1.write(strainList.get(i)+"\t"+count[i]);
                    bw1.newLine();
                }
                for(int i=0;i<barcodeList.size();i++){
                    bwR1[i].flush();bwR1[i].close(); 
                    bwR2[i].flush();bwR2[i].close();                     
                }
                bw1.flush();bw1.close();
                bw2.flush();bw2.close();
                br.close();br1.close();
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
    private void processMethodAndBarcode () {
        RowTable<String> t = new RowTable<>("/data1/home/junxu/rnaseq20181107/RNA-seq20181030.txt");
//        RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030-2.txt");
        Set<String> indexSet = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
           indexSet.add(t.getCell(i, 2).substring(1));
        }
        indexList = new ArrayList<>(indexSet);
        Collections.sort(indexList);
        barcodeLengths = new int[indexList.size()];//不同的样本可能带有不同长度的barcode
        barcodeLists = new ArrayList[indexList.size()];
        barcodeMethodMaps = new HashMap[indexList.size()];
        int[] cnts = new int[indexList.size()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = Collections.binarySearch(indexList, t.getCell(i, 2).substring(1));
            cnts[index]++;
            
        }
        for (int i = 0; i < cnts.length; i++) {
            barcodeLists[i] = new ArrayList<>();
            barcodeMethodMaps[i] = new HashMap<>();
        }
        int index = -1;
        for (int i = 0; i < t.getRowNumber(); i++) {
            index = Collections.binarySearch(indexList, t.getCell(i, 2).substring(1));
            allMethodList.add(t.getCell(i, 0));
            if(!t.getCell(i, 1).equals("")){
                barcodeLists[index].add(t.getCell(i, 1).substring(1));
                barcodeMethodMaps[index].put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
                barcodeLengths[index] = t.getCell(i, 1).length();
            }else{
                barcodeMethodMaps[index].put(t.getCell(i, 2).substring(1), t.getCell(i, 0));  
            }
            
        }
        Collections.sort(allMethodList);
    }
    public void parseParameters(String infileS) {
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
        this.inputDirS = pLineList.get(0);
        this.outputDirS = pLineList.get(1);
//        this.processTaxaAndBarcode();
    }
}
