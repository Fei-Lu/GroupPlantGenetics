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
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.lang.ProcessBuilder.Redirect;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xujun
 */
public class RNAseq {
    public RNAseq() throws IOException, FileNotFoundException {
//            this.getgenename();
//          this.genenumber();
//         this.notSDI();
//           this.test5();
//            this.parseFile()           
//            this.parseReads();
//            this.getGeneInfor();
 //           this.getGTFInfor();
//            this.getEachGeneInfor();
//            this.contGeneLength();
//            this.compare();
//                this.compareBarcode();
//            this.RNAGnene();
//            this.getGeneNameFromGTF();
//            this.caculateRPKM();
            this.diffEx();
//            this.caculateTPM();
    }
    public void diffEx(){
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/CountDifferent.txt";
        String outputfile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/expressionDiff.txt";
        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        xj.sortGeneByName();
        String temp=null;int index=0;int geneLength=0;
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputfile);
            bw.write(br.readLine());
            bw.newLine();
            while((temp=br.readLine())!=null){
                if(!temp.split("\t")[1].equals("0.0000")|!temp.split("\t")[2].equals("0.0000")){
                    index=xj.getGeneIndex(temp.split("\t")[0]);
                    geneLength=xj.getGeneLength(index);
                    bw.write(temp+"\t"+geneLength);bw.newLine();
                }
            }
            bw.flush();bw.close();
//            bw1.flush();bw1.close();  
            
            
        }
        catch (Exception ex) { 
            System.out.println(temp);
            System.out.print(index);
           ex.printStackTrace();
        }
    }
    public void caculateRPKM(){
        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        xj.sortGeneByName();
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/CountTable.txt";
        String outputFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/CountDifferent.txt";
        String output="/Users/xujun/Desktop/TEP/TEPOut/geneCount/lengthOfEachGene.txt";
        RowTable<String> rt = new RowTable<>(inputFile);
        int exonLength=0;int index=0;double TPM=0.00;
        double TPMwB=0.000;int cont=0;int wBcont=0;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            BufferedWriter bw1 = IOUtils.getTextWriter(output);
            bw.write("\t"+"WithoutB"+"\t"+"barcode");
            bw.newLine();
            for(int i=0;i<rt.getRowNumber();i++){
                index=xj.getGeneIndex(rt.getCell(i, 0));
                exonLength=xj.getGeneLength(index);
                bw1.write(rt.getCell(i, 0)+"\t"+exonLength);
                bw1.newLine();
                wBcont=rt.getCellAsInteger(i, 1);cont=rt.getCellAsInteger(i, 2);
                TPMwB=wBcont/exonLength/5.181;TPM=cont/exonLength/2.223;
                bw.write(rt.getCell(i, 0)+"\t"+TPMwB+"\t"+TPM);
                bw.newLine();
            }
            bw.flush();bw.close();
            bw1.flush();bw1.close();
        }
        catch (Exception ex) { 
           ex.printStackTrace();
        }
        
    }
    public void caculateTPM(){
        String gtfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        GeneFeature xj=new GeneFeature(gtfFile);
        xj.sortGeneByName();
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/CountTable.txt";
        String outputFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/TPM.txt";
        RowTable<String> rt = new RowTable<>(inputFile);
        int exonLength=0;int index=0;double TPM=0.00;
        double TPMwB=0.000;int cont=0;int wBcont=0;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            bw.write("\t"+"WithoutB"+"\t"+"barcode");
            bw.newLine();
            for(int i=0;i<rt.getRowNumber();i++){
                index=xj.getGeneIndex(rt.getCell(i, 0));
                exonLength=xj.getGeneLength(index);
                wBcont=rt.getCellAsInteger(i, 1);cont=rt.getCellAsInteger(i, 2);
                TPMwB=wBcont/exonLength;TPM=cont/exonLength;
                bw.write(rt.getCell(i, 0)+"\t"+TPMwB+"\t"+TPM);
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch (Exception ex) { 
           ex.printStackTrace();
        }
        
    }
    public void getGeneNameFromGTF(){
        String inputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        String outputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneName.txt";
        RowTable<String> rt = new RowTable<>(inputFile);
        HashSet<String> geneNameSet = new HashSet();
        String geneName=null;
        int startsite=0;
        int endsite=0;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            for(int i=0;i<rt.getRowNumber();i++){
                if(rt.getCell(i, 2).equals("exon")){
                    geneName=rt.getCell(i, 8).split(" ")[3].substring(6);
                    if(!(geneNameSet.contains(geneName))){
                        geneNameSet.add(geneName);
                        
                        if(endsite!=0){bw.write(endsite);}
                        startsite=Integer.valueOf(rt.getCell(i, 3));
                        bw.write(rt.getCell(i, 0)+"\t"+geneName+"\t"+startsite+"\t");
                        if(i!=0){
                            endsite=Integer.valueOf(rt.getCell(i-1, 4));
                            bw.write(endsite);
                        }
                    }
                }            
            }
        }
        catch (Exception ex) { 
           ex.printStackTrace();
        }
        
    }
    public void RNAGnene(){
        String inputfile="/Users/xujun/Desktop/TEP/TEPOut/geneTranscript/TEP-transcriptOfGene.txt";
        String input="/Users/xujun/Desktop/TEP/TEPOut/geneTranscript/TEPWithoutB-transcriptOfGene.txt";
        String output="/Users/xujun/Desktop/TEP/TEPOut/geneTranscript/TEPWithoutB-RNAGene.txt";
        String outputfile="/Users/xujun/Desktop/TEP/TEPOut/geneTranscript/TEP-RNAGene.txt";
        String temp=null;
        String name=null;
        int cont=0;
        try{
                BufferedReader br = IOUtils.getTextReader(inputfile);
                BufferedWriter bw = IOUtils.getTextWriter(outputfile);
                while((temp= br.readLine()) != null){
                    name=temp.split("\t")[0];
                    if(name.contains("ENSRNA")){
                        bw.write(temp);
                        bw.newLine();
                        cont++;
                    }                                                   
                }
                bw.write(cont);
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
    }
    public void getGeneName(){
        String gff3="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
        String outputfile="/Users/xujun/Desktop/TEP/GeneNewS.txt";
        String output="/Users/xujun/Desktop/TEP/GeneName.txt";
        GeneFeature xj=new GeneFeature(gff3);
        String temp=null;
        int cont=0;
        String news=null;
        String name=null;
        try{
                BufferedReader br = IOUtils.getTextReader(gff3);
                BufferedWriter bw = IOUtils.getTextWriter(outputfile);
                BufferedWriter bw1 = IOUtils.getTextWriter(output);
                bw1.write("SampleName");bw1.newLine();
                int chro=0;int startsite=0;
                while((temp= br.readLine()) != null){
                    if(temp.equals("###")){
                        if((news=br.readLine())!=null){
                            chro=Integer.valueOf(news.split("\t")[0]);startsite=Integer.valueOf(news.split("\t")[3]);
                            name=news.split("\t")[8].split(";")[0].substring(8);
                            bw.write(news.split("\t")[0]+"\t"+cont+"\t"+name+"\t"+news.split("\t")[3]+"\t"+news.split("\t")[4]);
                            bw.newLine();
                            cont++;                            
                            bw1.write(news.split("\t")[8].split(";")[0].substring(8));
                            bw1.newLine();
                        }
                       
                    }                                                  
                }
                br.close();
                bw.flush();bw1.flush();bw.close();bw1.close();
            }
            catch (Exception ex) { 
               System.out.println(news);
               ex.printStackTrace();
            }
    }
    public void compareGene(){
        String inputfile="/Users/xujun/Desktop/TEP/GeneName.txt";
        String input="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneName.txt";
        String outputfile="/Users/xujun/Desktop/TEP/compare.txt";
        RowTable<String> t = new RowTable<>(inputfile);
        List<String> nameList=new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            nameList.add(t.getCell(i, 0));
        }
        String temp=null;
        String name=null;
        try{
                BufferedReader br = IOUtils.getTextReader(inputfile);
                BufferedReader br1 = IOUtils.getTextReader(input);
                BufferedWriter bw = IOUtils.getTextWriter(outputfile);
                int cont=0;
                while((temp= br1.readLine()) != null){
                    name=temp.split("\t")[0];
                    if(nameList.contains(name)){
                        int index=nameList.indexOf(name);
                        bw.write(index+"\t"+name+"\t"+"both"+"\t"+cont);
                        bw.newLine();
                    }else{
                        bw.write("null");
                        bw.newLine();
                    }  
                    cont++;
                }
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
    }
    public void compare(){
        String inputfile="/Users/xujun/Desktop/TEP/GeneName.txt";
        String input="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneName.txt";
        String outputfile="/Users/xujun/Desktop/TEP/compare.txt";
        String output="/Users/xujun/Desktop/TEP/TEPOut/geneTranscript/Compare.txt";
        RowTable<String> t = new RowTable<>(inputfile);
        List<String> nameList=new ArrayList<>();
        List<String> w = new ArrayList<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            nameList.add(t.getCell(i, 0));
            w.add("without");
        }
        String temp=null;
        String name=null;
        int cont=0;
            try{
                BufferedReader br = IOUtils.getTextReader(input);
                BufferedWriter bw = IOUtils.getTextWriter(outputfile);
                BufferedWriter bw1 = IOUtils.getTextWriter(output);
                while((temp= br.readLine()) != null){
                    name=temp.split("\t")[0];
                    if(nameList.contains(name)){
                        int index=nameList.indexOf(name);
                        w.set(index, "both");
                    }else{
                        int index=nameList.indexOf(name);
                        w.set(index, "null");
                        cont++;
                        System.out.println(name);
                    }                                                   
                }
                t.addColumn("GTF", w);
                t.writeTextTable(outputfile, IOFileFormat.Text);
                bw1.write(cont);
                br.close();
                bw.flush();bw1.flush();bw.close();bw1.close();
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
    }
     public void compareBarcode(){
        String inputfile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/CountTable.txt";
        String output="/Users/xujun/Desktop/TEP/over.txt";
        RowTable<String> t = new RowTable<>(inputfile);
        int cont=0;
        int cont1=0;
            try{
                BufferedWriter bw = IOUtils.getTextWriter(output);
                for(int i=0;i<t.getRowNumber();i++){
                    if(Integer.parseInt(t.getCell(i, 1))>Integer.parseInt(t.getCell(i, 2))){
                       if(Integer.parseInt(t.getCell(i, 2))==0){
                           bw.write(t.getCell(i, 0)+"\t"+"Without");
                           cont1++;
                       }else{
                           bw.write(t.getCell(i, 0));
                            bw.newLine();
                            cont++;
                       }
                       
                    }
                }
                bw.flush();bw.close();
                System.out.println(cont);
                System.out.println(cont1);
            }
            catch (Exception ex) { 
               ex.printStackTrace();
            }
    }
    public void contGeneLength(){
        String inputFile="/Users/xujun/Desktop/TEP/GeneNewS.txt";
        String outputfile="/Users/xujun/Desktop/TEP/lengthOfEachChro.txt";
        String temp=null;
        int[] length=new int[12];
        int chro=0;
        int startsite=0;
        int endsite=0;
        int onelength=0;
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputfile);
            while((temp=br.readLine())!=null){  
                chro=Integer.parseInt(temp.split("\t")[1]);
                startsite=Integer.parseInt(temp.split("\t")[3]);
                endsite=Integer.parseInt(temp.split("\t")[4]);
                onelength=endsite-startsite+1;
                length[chro-1]+=onelength;
            }
            
            for(int i=0;i<length.length;i++){
                bw.write((i+1)+"\t"+length[i]);
                bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }
    }
    private HashIntIntMap hashIntIntOfGene1(){//key site value index .to establish the HashIntIntMap of one chrom(1)   
        String geneNewS="";
        int size = 306921361;
        int[] key = new int[size];
        int[] value = new int[size];
        String temp=null;
        int chro=0;
        int index=0;
        int startsite=0;
        int endsite=0;
        try{
            BufferedReader br = IOUtils.getTextReader(geneNewS);
            int cont=0;
            while((temp=br.readLine())!=null){ 
                    if(Integer.parseInt(temp.split("\t")[1])==1){
                        index=Integer.parseInt(temp.split("\t")[2]);
                        startsite=Integer.parseInt(temp.split("\t")[3]);
                        endsite=Integer.parseInt(temp.split("\t")[4]);
                        for (int i = startsite; i < endsite; i++) {
                            key[cont] = i;
                            value[cont] = index;
                            cont++;
                        }
                    }
                                                                                  
            }            
            br.close();
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }
        HashIntIntMap posGeneMaps = HashIntIntMaps.newImmutableMap(key,value);
        return posGeneMaps;        
    }
    public void getGTFInfor(){//read frome pgf and get information of each gene(name chro index startsite endsite)
        String pgfFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/Zea_mays.AGPv4.38.modified.gtf";
        String outputFile="/Users/xujun/Desktop/TEP/GeneNewS.txt";
        GeneFeature fs=new GeneFeature(pgfFile);
        String temp=null;
        int GeneStart=0;
        int GeneEnd=0;
        int index=0;
        int chro=0;
        try{
            BufferedReader br = IOUtils.getTextReader(pgfFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            while((temp=br.readLine())!=null){  
                if(temp.split("\t")[0].equals("Gene")){
 //                   index=fs.getGeneIndex(Integer.parseInt(temp.split("\t")[2]),Integer.parseInt(temp.split("\t")[3]));
                    index=fs.getGeneIndex(temp.split("\t")[1]);
                    bw.write(temp.split("\t")[1]+"\t"+temp.split("\t")[2]+"\t"+index+"\t"+temp.split("\t")[3]+"\t"+temp.split("\t")[4]);
                    bw.newLine();
                }                
            }
            br.close();
            bw.flush();bw.close();
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }      
    }
     public void getEachGeneInfor(){//read frome pgf and get information of each gene(name chro index startsite endsite)
        String inputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/GeneNewS.txt";
        RowTable<String> rt = new RowTable<>(inputFile);
        String outputFile="/Users/xujun/Desktop/RNA_seq/twice/STAR/EachGeneNewS.txt";
        String temp=null;
        String next=null;
        int GeneStart=0;
        int GeneEnd=0;
        int index=0;
        String chro="1";
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            for(int i=0;i<rt.getRowNumber();i++){
                if(rt.getCell(i, 1).equals(chro)){}
                else{
                    chro=rt.getCell(i, 1);
                    for(int j=0;j<rt.getRow(i).size();j++){
                        bw.write(rt.getCell(i-1, j)+"\t");
                    }
                    bw.newLine();
                }
                if(i==rt.getRowNumber()-1){
                    for(int j=0;j<rt.getRow(i).size();j++){
                        bw.write(rt.getCell(i, j)+"\t");
                    }
                }
                    
                
            }

            bw.flush();bw.close();
        }
        catch (Exception ex) { 
               ex.printStackTrace();              
        }      
    }
    public void test5 () throws IOException {
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/twice/clean1";
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
 //           String outfile=new File(outputDirS,p+".fq").getAbsolutePath();
            String seq = null;
            String seq1=null; 
            String temp=null;
            String seq2=null;
            String temp2=null;
            int cutIndex=0;
            try {
                BufferedReader br1 = IOUtils.getTextReader(infile1);
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
  //              BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)+".fq").getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                                                 
                //String seq=null;
                int cnt = 0;
                while ((temp = br1.readLine()) != null){
                    seq=br1.readLine(); 
                    Integer index = barcodeIndexMap.get(seq.substring(0, 8));
                    if (index == null) {
                        br1.readLine();br1.readLine();
                        continue;
                    }
                    int i=0;//用charat的方法 运行31s 试了两次都是30s左右
                    for(int a=9;a<seq.length();a++){
                        if (a+4 > seq.length()) {
                                br1.readLine();br1.readLine();
                                continue;
                            }
                            byte[] seqB = seq.substring(a, a+4).getBytes();
                            int value = this.getsum(seqB);
                            //317 = TTAA
                            if (value < 317) {
                                cutIndex = a;
                                bws[index].write(temp + "\n");
                                bws[index].write(seq + "\n");
                                bws[index].write(br1.readLine() + "\n");
                                bws[index].write(br1.readLine()+"\n");
                                cnt++;
                                if (cnt%100000 == 0) System.out.println("Wrote " + String.valueOf(cnt) + " sequences. " );
                                break;
                            }      
                    }
                    
                }
                br1.close();
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
    private int getsum(byte[] seqB){
        int sum=0;
        for(int i=0;i<seqB.length;i++){           
            sum+=seqB[i];
        }
        return sum;
    }
    public void getgenename(){
         String pfffile="/Users/xujun/Desktop/Zea_mays.AGPv4.38.pgf";
         GeneFeature fs=new GeneFeature(pfffile);         
         String inputfile="/Users/xujun/Desktop/TEP.text"; 
         String outputfile="/Users/xujun/Desktop/TEPname.text";
         BufferedReader br1 = IOUtils.getTextReader(inputfile);
         BufferedWriter bw = IOUtils.getTextWriter(outputfile);
         String a=null;
         int index=0;
         String name=null;
         String yon=null;
         try{
            while((a=br1.readLine())!=null){
                index=fs.getGeneIndex(Integer.parseInt(a.split("\t")[0]),Integer.parseInt(a.split("\t")[1]) );
                if (index < 0) continue;
                name=fs.getGeneName(index);
                bw.write(name+"\t"+a.split("\t")[0]+"\t"+a.split("\t")[1]+"\t"+(Integer.parseInt(a.split("\t")[1])+a.split("\t")[2].length())+"\t"+a.split("\t")[2]+"\t"+a.split("\t")[3]);               
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
    public void genenumber(){//对map到的每个基因上的reads进行计数 第一列基因名 第二列计数 第三列平均长度 再加一个第四列测序质量值的均值
         String inputfile="/Users/xujun/Desktop/TEPname.text";
         String outputfile="/Users/xujun/Desktop/TEPgenecount.text";
         String news=null;
         List w=new ArrayList();
         int count []=new int[20000];
         List<String>[] baseseqs=null;
         int averagelength[]=new int[20000];
         int basenumber[]=new int[20000];
         int averagephred[]=new int[20000];
         int chro[]=new int[20000];
         String phred=null;         
         int location=0;
         int sum =0;
         baseseqs = new List[20000];
         for (int i = 0; i < baseseqs.length; i++) {
             baseseqs[i] = new ArrayList<>();
         }
         try{
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);              
             while((news= br1.readLine()) != null){
                int phred1=0;
                if(!w.contains(news.split("\t")[0])){
                    w.add(news.split("\t")[0]);                    
                }               
                location=w.indexOf(news.split("\t")[0]);
                chro[location]=Integer.parseInt(news.split("\t")[1]);
                count[location]+=1;
                averagelength[location]+=news.split("\t")[4].length();
                List<String> baseseq=baseseqs[location];
                baseseq.add(news.split("\t")[2]+"\t"+news.split("\t")[3]+"\t"+news.split("\t")[4]);
                phred=news.split("\t")[5];
                for(int i=0;i<phred.length();i++){
                    phred1=phred1+(int)phred.charAt(i)-33;
                }
                basenumber[location]=basenumber[location]+phred.length();
                averagephred[location]+=phred1;
             }
             
             for(int i=0;i<w.size();i++){
                bw.write((String)w.get(i)+"\t"+chro[i]+"\t"+count[i]+"\t"+averagelength[i]/count[i]+"\t"+averagephred[i]/basenumber[i]);
                bw.newLine();
                for(int j=0;j<baseseqs[i].size();j++){
                    bw.write(baseseqs[i].get(j)+"\n");                              
                }
                
             }
             for(int i=0;i<count.length;i++){
                 sum=sum+count[i];
             }
             System.out.println(sum);
             br1.close();
             bw.flush();bw.close();           
          }
         catch(Exception ex){
             System.out.println();
             ex.printStackTrace();
         }
     }
    public void notSDI(){
         String inputfile1="/Users/xujun/Desktop/TEPnew.text";
         String outputfile="/Users/xujun/Desktop/TPE.text";
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
            BufferedReader br1 = IOUtils.getTextReader(inputfile1);
            BufferedWriter bw = IOUtils.getTextWriter(outputfile); 
            while ((a= br1.readLine()) != null){ 
                startpos=a.split("\t")[1];
                chro=a.split("\t")[0];
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
                                    bw.write(chro+"\t");bw.write(startpos1+"\t");bw.write(seq1+"\t");bw.write(q1+"\t");
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
//                                            base[Integer.parseInt(chro)][cont4]=seq.substring(cont2+cont,cont2+i);
//                                            SNPinsertion[Integer.parseInt(chro)][cont4]=String.valueOf(Integer.parseInt(startpos)+cont3);
//                                            cont4++;
//                                              cont2=cont2+Integer.parseInt(cigar.substring(cont,i));
//                                             cont=i+1;
                                        }else{
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
            br1.close();
            bw.flush();bw.close();
         }
         catch(Exception ex){
             System.out.println(cigar);
             ex.printStackTrace();
         }
     }
     private void getGeneInfor(){
        String outputDirS="/Users/xujun/Desktop";
        String subFqDirS = new File (outputDirS).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        int numCores = Runtime.getRuntime().availableProcessors();
        fList.stream().forEach(f -> { 
//            BufferedWriter bw = utils.IOUtils.getTextWriter(new File(subFqDirS,f.getName().replaceFirst(".sam", ".txt")).getAbsolutePath()); 
            StringBuilder sb = new StringBuilder();           
            sb.append("cut -f 3,4,6,10,11 ").append(f).append(" > ").append(new File(new File(subFqDirS,f.getName().replaceFirst(".sam", ".txt")).getAbsolutePath()));
            String command = sb.toString();
            System.out.println(command);
/*            try {
                Runtime rt = Runtime.getRuntime();                
                Process p = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    bw.write(temp);
                }
                p.waitFor();
            }*/
            try{
                Runtime rt = Runtime.getRuntime();
//                Process p = rt.exec(command);
                
                ProcessBuilder builder = new ProcessBuilder("cut","-f 3,4,6,10,11");
//                builder.redirectOutput();
                builder.redirectOutput(new File(new File(subFqDirS,f.getName().replaceFirst(".sam", ".txt")).getAbsolutePath()));
                builder.redirectError(new File(new File(subFqDirS,f.getName().replaceFirst(".sam", ".txt")).getAbsolutePath()));
                builder.start();
                builder.wait();
               
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
    private void parseReads(){//get the uniquely mapped and the multimapped rate <10 reads-get f.sam file
            String subFqDirS="/Users/xujun/Desktop/TEP";
            String f="/Users/xujun/Desktop/TEP.bam";
            String samtoolsPath="/usr/local/bin/samtools";
            StringBuilder sb = new StringBuilder();           
            sb.append(samtoolsPath).append(" view -F 0x100 ").append(f);//.append(" >> ").append("TEP.sam");
//            sb.append("ls ");
            String command = sb.toString();
            System.out.println(command);
            try {
//                File dir = new File(subFqDirS).getAbsolutePath();
//                File dir = new File("/Users/xujun/Desktop/TEP");
//                String []cmdarry ={"/bin/bash","-c", command+" >> out.txt"};               
//                String []cmdarry ={"cut", " -c "," -f 3 ",f, " >> ","output.txt"};                //String []cmdarry ={samtoolsPath, " view -F 0x100", f};

/*                for(int i=0;i<cmdarry.length;i++){
                    System.out.println(cmdarry[i]);
                }*/
                
                Runtime rt = Runtime.getRuntime();
 //               Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
//                p.waitFor();
 //               Process p=rt.exec(cmdarry,null.dir);
 //               Runtime rt = Runtime.getRuntime();              
                Process p = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    System.out.println(temp);
                }
                p.waitFor();
/*                ProcessBuilder builder = new ProcessBuilder(p.toString());
//                builder.redirectOutput();
                builder.redirectOutput(new File(new File(subFqDirS,f.replaceFirst(".bam", ".sam")).getAbsolutePath()));
                builder.redirectError(new File(new File(subFqDirS,f.replaceFirst(".bam", ".sam")).getAbsolutePath()));
                builder.start();
                builder.wait();
                builder.wait();*/
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);

        
    }
    public void parseFile(){
        String samFile="/Users/xujun/Desktop/TEP.sam";
        String output="/Users/xujun/Desktop/TEPnew.text";
        RowTable<String> t = new RowTable<>(samFile);
        try{            
            BufferedWriter bw = IOUtils.getTextWriter(output);
            for(int i=17;i<t.getRowNumber();i++){
                bw.write(t.getCell(i, 2)+"\t"+t.getCell(i, 3)+"\t"+t.getCell(i, 5)+"\t"+t.getCell(i, 9)+"\t"+t.getCell(i, 10));
                bw.newLine();
            }
            bw.flush();bw.close();
            
        }
        catch(Exception ex){
//             System.out.println(cigar);
             ex.printStackTrace();
        }
        
        
    }
    public static void main(String[] args) throws IOException, FileNotFoundException{
        new RNAseq();
    }
    
}
