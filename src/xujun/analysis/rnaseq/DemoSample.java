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
            this.phred();
//            this.getgenename();
//            this.geneposition();
//            this.sort();
//           this.genenumber();
//          this.merge();
//          this.notSDI();
            this.test5();
    }      
public void test5 () throws IOException {
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/TEP/fastq";
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
                            bws[index].write(br2.readLine() + "\n");
                            bws[index].write(br2.readLine() + "\n");
                            bws[index].write(br2.readLine() + "\n");
                            bws[index].write(br2.readLine()+"\n");
                        
                        
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
         String inputfile="/Users/xujun/Desktop/news.fq"; 
         String outputfile="/Users/xujun/Desktop/notSIDname.text";
         BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
         BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
         String a=null;
         int index=0;
         String name=null;
         String yon=null;
         try{
            while((a=br1.readLine())!=null){
                index=fs.getGeneIndex(Integer.parseInt(a.split("\t")[0]),Integer.parseInt(a.split("\t")[1]) );
                if (index < 0) continue;
                name=fs.getGeneName(index);
                bw.write(name+"\t"+a.split("\t")[2]+"\t"+a.split("\t")[3]);               
                bw.newLine();
                }
            br1.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             System.out.println(a);
             ex.printStackTrace();
         }
         
         
         
     }
     public void geneposition() { //get the gene chr and start site  and get the mapped length
        String inputFile1 ="/Users/xujun/Desktop/aligned.text";
        String inputFile2 ="/Users/xujun/Desktop/pos.text";
        String inputFile3 ="/Users/xujun/Desktop/chromosome.text";
        String inputFile4="/Users/xujun/Desktop/normal.fq";
        String outputFile = "/Users/xujun/Desktop/containSDI/geneposition.text"; 
        String start=null;
        String position=null;
        int end=0;
        int start1=0;
        int count=0;
        int count2=0;
        String a=null;
        String seq=null;
        String seq1=null;
        String phred=null;
        String phred1=null;
        try {
                BufferedReader br1 = utils.IOUtils.getTextReader(inputFile1);   
                BufferedReader br2 = utils.IOUtils.getTextReader(inputFile2);
                BufferedReader br3 = utils.IOUtils.getTextReader(inputFile3);
                BufferedReader br4 = utils.IOUtils.getTextReader(inputFile4);
                BufferedWriter bw = utils.IOUtils.getTextWriter(outputFile);
                List w=new ArrayList();
                while ((a= br1.readLine()) != null){ 
                     count++;
                     int cont=0;
                     int cont2=0;
                     int cont3=0;
                     start =br2.readLine();
                     position=br3.readLine();
                     for(int i=0;i<a.length();i++){                         
                        if (!Character.isDigit(a.charAt(i))){ 
                            if(a.charAt(i)=='N'){ 
                                bw.write(position+"\t"+start+"\t"+cont2);
                                bw.newLine();                                                                     
                                cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                cont3=cont2;
                                cont=i+1;
                                start=String.valueOf(Integer.parseInt(start)+cont2);
                                count2++;
                             }else{                                    
                                   cont2=cont2+Integer.parseInt(a.substring(cont,i));
                                   cont=i+1;
                                   if(i==a.length()-1){
                                       bw.write(position+"\t"+start+"\t"+(cont2-cont3));
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
     public void genenumber(){//对map到的每个基因上的reads进行计数 第一列基因名 第二列计数 第三列平均长度 再加一个第四列测序质量值的均值
         String inputfile="/Users/xujun/Desktop/notSIDname.text";
         String outputfile="/Users/xujun/Desktop/notSDIgenecount.text";
         String news=null;
         List w=new ArrayList();
         int count []=new int[40000];
         int averagelength[]=new int[40000];
         int basenumber[]=new int[40000];
         int averagephred[]=new int[40000];
         String phred=null;         
         int location=0;
         int sum =0;
         try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);            
             while((news= br1.readLine()) != null){
                int phred1=0;
                if(!w.contains(news.split("\t")[0])){
                    w.add(news.split("\t")[0]);                    
                }
                location=w.indexOf(news.split("\t")[0]);
                count[location]=count[location]+1;
                averagelength[location]=averagelength[location]+Integer.parseInt(news.split("\t")[1]);
                phred=news.split("\t")[2];
                for(int i=0;i<phred.length();i++){
                    phred1=phred1+(int)phred.charAt(i)-33;
                }
                basenumber[location]=basenumber[location]+phred.length();
                averagephred[location]=averagephred[location]+phred1;
             }
             for(int i=0;i<w.size();i++){
                bw.write((String)w.get(i)+"\t"+count[i]+"\t"+averagelength[i]/count[i]+"\t"+averagephred[i]/basenumber[i]);
                bw.newLine();
             }
             for(int i=1;i<count.length;i++){
                 sum=sum+count[i];
             }
             System.out.println(sum);
             br1.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }
     }
     public void merge(){
        String inputFile1 ="/Users/xujun/Desktop/chro(MD).text";
        String inputFile2 ="/Users/xujun/Desktop/start(MD).text";
        String inputFile3 ="/Users/xujun/Desktop/CIGAR(MD).text";
        String inputFile4="/Users/xujun/Desktop/MD.text";
        String outputfile="/Users/xujun/Desktop/news(MD).text";
        String news=null;
        try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputFile1);
             BufferedReader br2 = utils.IOUtils.getTextReader(inputFile2);
             BufferedReader br3 = utils.IOUtils.getTextReader(inputFile3);
             BufferedReader br4 = utils.IOUtils.getTextReader(inputFile4);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);            
             while((news= br1.readLine()) != null){
                bw.write(news+"\t"+br2.readLine()+"\t"+br3.readLine()+"\t"+br4.readLine().substring(5)+"\n");
             }
             
             br1.close();br2.close();br3.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }
     }
     public void notSDI(){
         String inputfile1="/Users/xujun/Desktop/news(MD).text";
         String inputfile2="/Users/xujun/Desktop/normal.fq";
         String outputfile="/Users/xujun/Desktop/news.fq";
         String a=null;
         String start=null;
         String seq=null;
         String news=null;
         String MD=null;
         String q=null;
         String seq1=null;                   
         String q1=null;
         String startpos=null;
         String chro=null;
         String cigar=null;
         String [][] base=new String[13][];
         String SNPinsertion [][]=new String[13][];
         try{
            BufferedReader br1 = utils.IOUtils.getTextReader(inputfile1);
            BufferedReader br2 = utils.IOUtils.getTextReader(inputfile2);
            BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile); 
            while ((a= br1.readLine()) != null){ 
                startpos=a.split("\t")[1];
                chro=a.split("\t")[0];
                cigar=a.split("\t")[2];
                MD=a.split("\t")[3];
                start=br2.readLine();
                seq=br2.readLine();
                news=br2.readLine();
                q=br2.readLine();
                int cont=0;
                int cont2=0;
                int cont3=0;
                int cont4=0;
                for(int i=0;i<cigar.length();i++){                           
                            if (!Character.isDigit(cigar.charAt(i))){ 
                                if(cigar.charAt(i)=='M'){
                                    startpos=String.valueOf(Integer.parseInt(startpos)+cont3);
                                        seq1=seq.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));  
                                        q1=q.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));
                                        bw.write(chro+"\t");bw.write(startpos+"\t");bw.write(seq1.length()+"\t");bw.write(q1);
                                        bw.newLine();
                                        cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                        cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                }else{
                                    if(cigar.charAt(i)=='D'|cigar.charAt(i)=='N'){
                                        cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                        cont=i+1;
                                    }else{   
                                        if(cigar.charAt(i)=='I'){
//                                            base[Integer.parseInt(chro)][cont4]=seq.substring(cont2,cont2+Integer.parseInt(cigar.substring(cont,i)));
//                                            SNPinsertion[Integer.parseInt(chro)][cont4]=String.valueOf(Integer.parseInt(startpos)+cont3);
//                                            cont4++;
                                            cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                            cont=i+1;
                                        }else{
                                            cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
                                            cont3=cont3+Integer.parseInt(cigar.substring(cont,i));
                                            cont=i+1;
                                        }                                                                                 
                                        
                                    }
                                }                                
                            }else{
                                continue;
                            } 
                             
                        }                              
            }
            br1.close();br2.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             System.out.println(cigar);
             System.out.println(MD);
             ex.printStackTrace();
         }
     }
    public void genenumber2(){//对map到的每个基因上的reads进行计数 第一列基因名 第二列计数 第三列平均长度 再加一个第四列测序质量值的均值
         String inputfile="/home/aoyue/xujun/notSIDname.text";
         String outputfile="/home/aoyue/xujun/notSDIgenecount2.text";
         String news=null;
         List w=new ArrayList();
         int count []=new int[40000];
         int length[][]=new int[40000][55000];
         int basenumber[][]=new int[40000][55000];
         int allphred[][]=new int[40000][4000000];
         String phred=null;         
         int location=0;
         int sum =0;
         try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);            
             while((news= br1.readLine()) != null){
                int phred1=0;
                if(!w.contains(news.split("\t")[0])){
                    w.add(news.split("\t")[0]);                    
                }
                location=w.indexOf(news.split("\t")[0]);                
                length[location][count[location]]=Integer.parseInt(news.split("\t")[1]);
                count[location]=count[location]+1;
                phred=news.split("\t")[2];
                for(int i=0;i<phred.length();i++){
                    phred1=phred1+(int)phred.charAt(i)-33;
                }
                basenumber[location][count[location]]=phred.length();
                allphred[location][count[location]]=phred1;
             }
             for(int i=0;i<w.size();i++){
                bw.write((String)w.get(i)+"\t"+count[i]+"\t"+length[i]+"\t"+allphred[i]);
                bw.newLine();
             }
             br1.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }
     }
    public static void main(String[] args) throws IOException, FileNotFoundException, BiffException{
        new DemoSample();
    }
}




