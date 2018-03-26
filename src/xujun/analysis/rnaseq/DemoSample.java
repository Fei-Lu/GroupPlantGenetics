/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

/**
 *
 * @author Jun Xu
 */
import com.koloboke.collect.map.hash.HashByteByteMap;
import format.genomeAnnotation.GeneFeature;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import jxl.read.biff.BiffException;
import utils.IOUtils;
public class DemoSample{
    public DemoSample() throws IOException, FileNotFoundException, BiffException {

      //    this.testL33();
//           this.mappedread();
//            this.phred();
            this.getgenename();
//            this.geneposition();
//            this.sort();
    }      
public void test5 () throws IOException {
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/twice/clean_data";
        String outputDirS = "/Users/xujun/Desktop/RNA_seq/twice/test-twice";        
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
            nameList.add(rt.getCell(i, 0));
        }
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        nameSet.parallelStream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_1.clean.fq").getAbsolutePath();
            String infile2 = new File (inputDirS, p+"_2.clean.fq").getAbsolutePath();
 //           String outfile=new File(outputDirS,p+".fq").getAbsolutePath();
            String seq = null;
            String seq1=null; 
            String temp=null;
            String seq2=null;
            String temp2=null;
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);
                BufferedReader br2 = utils.IOUtils.getTextReader(infile2);
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
  //              BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)+".fq").getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                
                 
                
                //String seq=null;
                while ((temp = br1.readLine()) != null){
                    int i=0;                 
                    seq=br1.readLine(); 
                    Integer index = barcodeIndexMap.get(seq.substring(0, 8));
                    if (index == null) {
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                        br1.readLine();br1.readLine();
                        continue;
                    }else{
                        br1.readLine();br1.readLine();
                        temp2=br2.readLine();
                        seq2=br2.readLine();
                        for(int a=0;a<seq2.length();a++){
                            if(seq2.charAt(a)=='A'){
                                i++;
                                if(i>10){
                                    seq1=seq2.substring(0, a+1-i);
                                    break;
                                }
                            }else{
                                i=0;
                            }

                        }
                        if(i<=0){
                            seq1=seq2;
                        }
                        if(seq1.length()==0){
                            br2.readLine();br2.readLine();
                        }else{
                            bws[index].write(temp2 + "\n");
                            bws[index].write(seq1 + "\n");
                            bws[index].write(br2.readLine() + "\n");
                            bws[index].write(br2.readLine().substring(0,seq1.length())+"\n");
                        }
                        
                    }
                    
                }
                br1.close();
                br2.close();
//                bw.flush();bw.close();
                for (int i = 0; i < bws.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                
                    
                       
        }
                
        catch (Exception ex) {
               System.out.println(seq+"\t1234");  
               System.out.println(seq1+"\t1234"); 
               ex.printStackTrace();
               
        }
        });      
     }
     public void testL35 () throws IOException {
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/L3-1/clean_data/L3-1";
        String outputDirS = "/Users/xujun/Desktop/RNA_seq/L3-1/clean_data/test";        
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
            nameList.add(rt.getCell(i, 0));
        }
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        nameSet.parallelStream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_R1.fq").getAbsolutePath();
            String outfile2 = new File (inputDirS, "555.fq").getAbsolutePath();
            String outfile=new File(outputDirS,p+".fq").getAbsolutePath();
            String seq = null;
            String seq1=null; 
            String temp=null;
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);                               
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)+".fq").getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                
                 
                List w=new ArrayList();
                //String seq=null;
                while ((temp = br1.readLine()) != null){
                    int i=0;                 
                    seq=br1.readLine(); 
                    //用charat的方法 运行31s
                    for(int a=0;a<seq.length();a++){
                        if(seq.charAt(a)=='A'){
                            i++;
                            if(i>10){
                                seq1=seq.substring(0, a+1-i);
                                w.add(a+1-i);
                                break;
                            }
                        }else{
                            i=0;
                        }
                        
                    }
                    if(i<=0){
                        seq1=seq;
                        w.add(150);
                    }
                    if(seq1.length()==0){
                        br1.readLine();br1.readLine();
                    }else{
                        bw.write(temp+ "\n");bw.write(seq1+ "\n");bw.write(br1.readLine() + "\n");
                        bw.write(br1.readLine().substring(0, seq1.length()) + "\n");
                    }
                    
                }
                br1.close();
                bw.flush();bw.close();
/*                for (int i = 0; i < bws.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }*/
                Collections.shuffle(w); 
                int randomSeriesLength = 1000;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      System.out.println(randomSeries.get(k)+ " "); 
                      if(k==randomSeriesLength){
                           System.out.println("******************************************************** ");
                       } 
                   } else {
                       System.out.print(randomSeries.get(k)+ " ");
                           
                   }
                }
                
                    
                       
        }
                
        catch (Exception ex) {
               System.out.println(seq+"\t1234");  
               System.out.println(seq1+"\t1234"); 
               ex.printStackTrace();
               
        }
        });      
     }
    public void testL33 () throws IOException {
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/L3-1/clean_data/L3-1";
        String outputDirS = "/Users/xujun/Desktop/RNA_seq/L3-1/clean_data/test";        
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
            nameList.add(rt.getCell(i, 0));
        }
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        nameSet.parallelStream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_R2.fq").getAbsolutePath();
            String outfile2 = new File (outputDirS, "333.fq").getAbsolutePath();
            String seq = null;
            String seq1=null; 
            String temp=null;
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);                              
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile2);
                List w=new ArrayList();
                //String seq=null;
                while ((temp = br1.readLine()) != null){
                    int i=0;                 
                    seq=br1.readLine(); 
                    for(int a=0;a<=seq.length();a++){
                        if(seq.charAt(a)=='T'){
                            i++;
                        }else{
                            if(i<8){
                                i=0;
                            }else{
                                seq1=seq.substring(a);
                                w.add(a);
                                break;
                            }
                        }
                    
                        
                    }
                    if(i==0){
                        seq1=seq;
                    }
                        bw.write(temp+ "\n");bw.write(seq1+ "\n");bw.write(br1.readLine() + "\n");
                        bw.write(br1.readLine().substring(seq1.length()) + "\n");
                    
                    
                    
                }
                br1.close();
                bw.flush();bw.close();
                Collections.shuffle(w); 
                int randomSeriesLength = 1000;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      System.out.println(randomSeries.get(k)+ " "); 
                      if(k==randomSeriesLength){
                           System.out.println("******************************************************** ");
                       } 
                   } else {
                       System.out.print(randomSeries.get(k)+ " ");
                           
                   }
                }
                
                    
                       
        }
                
        catch (Exception ex) {
               System.out.println(seq+"\t1234");  
               System.out.println(seq1+"\t1234"); 
               ex.printStackTrace();
               
        }
        });      
     }
    public void mappedread() {
        String inputFile1 ="/Users/xujun/Desktop/aligned.text";
        String inputFile2 ="/Users/xujun/Desktop/normal.fq";
        String outputFile = "/Users/xujun/Desktop/mappedreads.fq"; 
        String start=null;
        String news=null;
        String seq=null;
        String q=null;
        String a=null;        
        try {
                BufferedReader br1 = utils.IOUtils.getTextReader(inputFile1);   
                BufferedReader br2 = utils.IOUtils.getTextReader(inputFile2);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outputFile);
                List w=new ArrayList();
                while ((a= br1.readLine()) != null){                     
                        String seq1="";                   
                        String q1="";
                        start=br2.readLine();
                        seq=br2.readLine();
                        news=br2.readLine();
                        q=br2.readLine();
                        int cont=0;
                        int cont2=0;
                        for(int i=0;i<a.length();i++){                           
                            if (!Character.isDigit(a.charAt(i))){ 
                                if(a.charAt(i)=='M'){
                                        seq1=seq1.concat(seq.substring(cont2,cont2+Integer.parseInt(a.substring(cont,i))));  
                                        q1=q1.concat(q.substring(cont2,cont2+Integer.parseInt(a.substring(cont,i)))); 
                                        cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                        cont=i+1;
                                        w.add(seq1.length());
                                }else{
                                    if(a.charAt(i)=='D'|a.charAt(i)=='N'){
                                        cont=i+1;
                                    }else{
                                        cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                        cont=i+1;
                                    }
                                }                                
                            }else{
                                continue;
                            } 
                             
                        }
                        bw.write(start+"\n");
                        bw.write(seq1+"\n");
                        bw.write(news+"\n");    
                        bw.write(q1+ "\n");                              
                    }
                                                                               
                
                br1.close();
                br2.close();
                bw.flush();bw.close();
                Collections.shuffle(w); 
                int randomSeriesLength = 1000;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      System.out.println(randomSeries.get(k)+ " "); 
                      if(k==randomSeriesLength){
                           System.out.println("******************************************************** ");
                       } 
                   } else {
                       System.out.print(randomSeries.get(k)+ " ");
                           
                   }
                }
                System.out.println(w.size());
        }
                
        catch (Exception ex) { 
               System.out.println(seq+"长度"+seq.length());
               System.out.println(a+"12345");
               ex.printStackTrace();
               
        }      
        
        
    
            
    }
     public void phred() {
        String inputFile1 ="/Users/xujun/Desktop/mappedreads.fq";
        String outputFile = "/Users/xujun/Desktop/phredscore.text"; 
        String news=null; 
        String phred=null;
        int phred1=0;
        try {
                BufferedReader br1 = utils.IOUtils.getTextReader(inputFile1);   
                BufferedWriter bw = utils.IOUtils.getTextWriter(outputFile);
                
                while ((news= br1.readLine()) != null){  
                    List w=new ArrayList();
                    br1.readLine();br1.readLine();
                    phred=br1.readLine();
                    for(int i=0;i<phred.length();i++){
                        phred1=(int)phred.charAt(i)-33;
                        bw.write(phred1+"\t");
                    }                      
                    bw.newLine();
                         
                }
                                                                               
                
                br1.close();
                bw.flush();bw.close();
        }
                
        catch (Exception ex) { 
               ex.printStackTrace();
               
        }
     }
     public void getgenename(){
         String kfffile="/Users/xujun/Desktop/Zea_mays.AGPv4.38.kgf";
         GeneFeature fs=new GeneFeature(kfffile);
         String inputfile="/Users/xujun/Desktop/geneposition.text"; 
         String outputfile="/Users/xujun/Desktop/index.text";
         BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
         BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
         String a=null;
         int index=0;
         String name=null;
         String yon=null;
         try{
            while((a=br1.readLine())!=null){
                index=fs.getGeneIndex(Integer.parseInt(a.split("\t")[0]),Integer.parseInt(a.split("\t")[1]) );
                name=fs.getGeneName(index);
                bw.write(name);               
                bw.newLine();                                              
            }
            br1.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             ex.printStackTrace();
         }
         
         
         
     }
     public void geneposition() {
        String inputFile1 ="/Users/xujun/Desktop/aligned.text";
        String inputFile2 ="/Users/xujun/Desktop/index.text";
        String inputFile3 ="/Users/xujun/Desktop/chromosome.text";
        String outputFile = "/Users/xujun/Desktop/geneposition.text"; 
        String start=null;
        String position=null;
        int end=0;
        int start1=0;
        int count=0;
        int count2=0;
        String a=null;
        try {
                BufferedReader br1 = utils.IOUtils.getTextReader(inputFile1);   
                BufferedReader br2 = utils.IOUtils.getTextReader(inputFile2);
                BufferedReader br3 = utils.IOUtils.getTextReader(inputFile3);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outputFile);
                List w=new ArrayList();
                while ((a= br1.readLine()) != null){ 
                        count++;
                        int cont=0;
                        int cont2=0;
                        start =br2.readLine();
                        position=br3.readLine();
                        for(int i=0;i<a.length();i++){                         
                            if (!Character.isDigit(a.charAt(i))){ 
                                if(a.charAt(i)=='N'){                                   
                                   bw.write(position+"\t"+start);
                                   bw.newLine();                                                                     
                                   cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                   cont=i+1;
                                   start=String.valueOf(cont2);
                                   count2++;
                                }else{                                    
                                   cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                   cont=i+1;
                                   if(i==a.length()-1){
                                       bw.write(position+"\t"+start);
                                       bw.newLine();
                                       count2++;
                                   }
                                }                                
                            }else{
                                continue;
                            } 
                             
                        }                             
                    }
                                                                               
                
                br1.close();
                br2.close();
                br3.close();
                bw.flush();bw.close();
                System.out.println(count);
                System.out.println(count2);
        }
                
        catch (Exception ex) { 
               System.out.println();
               System.out.println(a+"12345");
               ex.printStackTrace();
               
        }      
        
        
    
            
    }
     public void sort(){
         String inputfile="/Users/xujun/Desktop/geneposition.text";
         String outputfile="/Users/xujun/Desktop/sortedgeneposition.text";
         try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
             ArrayList<String> list = new ArrayList<String>();
             String news=null;
             while((news= br1.readLine()) != null){                
                 list.add(news);                
             }
             Collections.sort(list.subList(0, 0));
             bw.write(Arrays.toString(list.toArray())+"\n");
         }
         catch(Exception ex){
             ex.printStackTrace();
         }
     }
     public void geneindex(){
         String kfffile="/Users/xujun/Desktop/Zea_mays.AGPv4.38.kgf";
         GeneFeature fs=new GeneFeature(kfffile);
     }
    
    public static void main(String[] args) throws IOException, FileNotFoundException, BiffException{
        new DemoSample();
    }
}




