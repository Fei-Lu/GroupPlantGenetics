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
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.*;
import java.util.*;

import static java.lang.Math.abs;
public class DemoSample{
    public DemoSample()  {

//          this.testL33();
//           this.mappedread();
//            this.phred();
            this.getgenename();
            this.geneposition();
//            this.sort();="
//           this.genenumber();
//          this.merge();
//          this.notSDI();
//            this.test5();
//            this.CompareTwoCountMethod();
//            this.LibraryCompare();
//            this.HTSeqCount();
            this.countCompare();
//            this.TPMHTSeq();
//            this.CoDectededGene();
//            this.testLexo();
//            this.pattern();
//            this.HTSeqCountDouble();
            this.splitBwaScript();
    }  
    public void splitBwaScript(){
        String infileS = "/Users/xujun/sh_htseq.sh";
        String outfileDirS = "/Users/xujun/Documents/splitScript/";
        String shfileS = "/Users/xujun/Documents/sh_htseq.sh";
        try{
        String[] outfileS= new String[32];
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter[] bw = new BufferedWriter[32];
        for(int i=0;i<outfileS.length;i++){
            String num = PStringUtils.getNDigitNumber(3, i+1);
            outfileS[i]=new File(outfileDirS,"htseq_20190709_" + num + ".sh").getAbsolutePath();
            bw[i] = IOUtils.getTextWriter(outfileS[i]);
            String temp;
            for(int j=0; j<6;j++){
                if((temp = br.readLine()) != null){
                    bw[i].write(temp);
                    bw[i].newLine();
                }
            }
            bw[i].flush();bw[i].close();
        }
        br.close();
    }
    catch(Exception e){
        e.printStackTrace();
        System.exit(1);
    }
     try{
        File[] fs = new File(outfileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sh");
        Arrays.sort(fs);
        BufferedWriter bw = IOUtils.getTextWriter(shfileS);
        for(int i=0; i<fs.length;i++){
            bw.write("sh " + fs[i].getName() + " &");
            bw.newLine();
        }
        bw.flush();bw.close();
    }
    catch(Exception e){
        e.printStackTrace();
        System.exit(1);
    }
   }
    public void HTSeqCountDouble(){
//        String inputDirS = new File ("/data1/home/junxu/wheat/doubleAll/sams").getAbsolutePath();
////        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
//        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
//        List<File> fList = Arrays.asList(fs);
//        fList.stream().forEach(f -> {  
//            StringBuilder ob = new StringBuilder();
//            ob.append("samtools view -h -q 255 -f 0x2 "+f+" > ");
//            ob.append(f.getName().replace("Aligned.out.sam", "UniqPE.sam"));
//            String commandob = ob.toString();
//            System.out.println(commandob);
//            try {
//                File dir = new File(new File ("/data1/home/junxu/wheat/doubleAll/UniqPE").getAbsolutePath());
//                String []cmdarryob ={"/bin/bash","-c",commandob};
//                Process pob=Runtime.getRuntime().exec(cmdarryob,null,dir);
//                pob.waitFor();
//                pob.destroy();
//            }
//            catch (Exception e) {
//                e.printStackTrace();
//            }
//            System.out.println("Finished"+f);
//        });
//        int bounds=0;int[] 
        File[] fs1 = new File("/data1/home/junxu/wheat/doubleAll/UniqPE").listFiles();
        fs1 = IOUtils.listFilesEndsWith(fs1, "UniqPE.sam");
        List<File> fList1 = Arrays.asList(fs1);
//        bounds = PArrayUtils.getSubsetsIndicesBySubsetSize(fs1.length,);
        fList1.stream().forEach(f -> {  
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s reverse ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat/rightchangewheat.gtf").append(" >> ");
//            sb.append(" "+this.gffFile).append(" >> ");
            sb.append(f.getName().replace("UniqPE.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File ("/data1/home/junxu/wheat/doubleAll/UniqPE/sh").getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();
                p.destroy();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
    public void TPMHTSeq(){
        String HTSeqFile="/Users/xujun/Desktop/TEP/TEPOut/countResult.txt";
        String gff="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
        GeneFeature gf = new GeneFeature(gff);
        RowTable t = new RowTable (HTSeqFile);
        double HNor11 [] = new double[t.getRowNumber()];double HNor12 [] = new double[t.getRowNumber()];double HNor13 [] = new double[t.getRowNumber()];double HNor14 [] = new double[t.getRowNumber()];
        double HTPM1 [] = new double [t.getRowNumber()];double HTPM2 [] = new double [t.getRowNumber()];double HTPM3 [] = new double [t.getRowNumber()];double HTPM4 [] = new double [t.getRowNumber()];
        double Hdepth1=0;double Hdepth2=0; double Hdepth3=0; double Hdepth4=0; 
        double HdepthN1=0;double HdepthN2=0; double HdepthN3=0; double HdepthN4=0;
        int[] HTSeq1=new int [t.getRowNumber()]; int[] HTSeq2=new int [t.getRowNumber()]; int [] HTSeq3=new int [t.getRowNumber()]; int [] HTSeq4=new int [t.getRowNumber()];
        List<String> geneName = new ArrayList <>();
        for(int i=0;i<t.getRowNumber();i++){
            geneName.add(t.getCellAsString(i, 0));
            int geneLength=gf.getGeneEnd(gf.getGeneIndex(t.getCellAsString(i, 0)))-gf.getGeneStart(gf.getGeneIndex(t.getCellAsString(i, 0)));
            int index=geneName.indexOf(t.getCell(i, 0));
            HTSeq1 [index]=t.getCellAsInteger(i, 17);HTSeq2 [index]=t.getCellAsInteger(i, 19);HTSeq3 [index]=t.getCellAsInteger(i, 18);HTSeq4 [index]=t.getCellAsInteger(i, 8);
            HNor11[index]=Math.log(t.getCellAsInteger(i, 17))-Math.log(geneLength);
            HNor12[index]=Math.log(t.getCellAsInteger(i, 19))-Math.log(geneLength);   
            HNor13[index]=Math.log(t.getCellAsInteger(i, 18))-Math.log(geneLength);   
            HNor14[index]=Math.log(t.getCellAsInteger(i, 8))-Math.log(geneLength); 
//            HNor11[index]=t.getCellAsInteger(i, 17)/geneLength;
//            HNor12[index]=t.getCellAsInteger(i, 19)/geneLength;
//            HNor13[index]=t.getCellAsInteger(i, 18)/geneLength;   
//            HNor14[index]=t.getCellAsInteger(i, 8)/geneLength; 
        }
        for(int i=0;i<HNor12.length;i++){
            Hdepth1+=Math.exp(HNor11[i]);
            Hdepth2+=Math.exp(HNor12[i]);
            Hdepth3+=Math.exp(HNor13[i]);
            Hdepth4+=Math.exp(HNor14[i]);
//            Hdepth1+=HNor11[i];
//            Hdepth2+=HNor12[i];
//            Hdepth3+=HNor13[i];
//            Hdepth4+=HNor14[i];
        }
        HdepthN1=Math.log(Hdepth1);HdepthN2=Math.log(Hdepth2);HdepthN3=Math.log(Hdepth3);HdepthN4=Math.log(Hdepth4);
        for(int i=0;i<HNor12.length;i++){
            HTPM1[i]=Math.exp(HNor11[i]-HdepthN1+Math.log(1e6));
            HTPM2[i]=Math.exp(HNor12[i]-HdepthN2+Math.log(1e6));
            HTPM3[i]=Math.exp(HNor13[i]-HdepthN3+Math.log(1e6));;
            HTPM4[i]=Math.exp(HNor14[i]-HdepthN4+Math.log(1e6));
//            HTPM1[i]=HNor11[i]/(Hdepth1/1000000);
//            HTPM2[i]=HNor12[i]/(Hdepth2/1000000);
//            HTPM3[i]=HNor13[i]/(Hdepth3/1000000);
//            HTPM4[i]=HNor14[i]/(Hdepth4/1000000);
        }
        try{
            BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/HTPM.txt");
            BufferedWriter bw1 = IOUtils.getTextWriter("/Users/xujun/HTPMPatern.txt");
//            bw.write("geneName"+"chr"+"s1116-1"+"\t"+"TPM"+"\t"+"PLATE-seq"+"\t"+"TPM"+"\t"+"Lexogen"+"\t"+"TPM"+"\t"+"SiPAS"+"\t"+"TPM");bw.newLine();
            bw.write("geneName"+"\t"+"chr"+"\t"+"s1116-1"+"\t"+"PLATE-seq"+"\t"+"Lexogen"+"\t"+"SiPAS");bw.newLine();
            bw1.write("geneName"+"\t"+"PLATE-seq");bw1.newLine();
            int count=0;int chr=0;
            for(int i=0;i<HNor12.length;i++){
                chr=gf.getGeneChromosome(gf.getGeneIndex(geneName.get(i)));
                if(HTSeq1[i]!=0||HTSeq2[i]!=0||HTSeq3[i]!=0){
                    bw.write(geneName.get(i)+"\t"+chr+"\t"+HTSeq1[i]+"\t"+HTSeq2[i]+"\t"+HTSeq3[i]+"\t"+HTSeq4[i]);bw.newLine();
                }
//                bw.write(geneName.get(i)+"\t"+chr+"\t"+HTSeq1[i]+"\t"+HTPM1[i]+"\t"+HTSeq2[i]+"\t"+HTPM2[i]+"\t"+HTSeq3[i]+"\t"+HTPM3[i]+"\t"+HTSeq4[i]+"\t"+HTPM4[i]);bw.newLine();
                if(HTPM2[i]>=1){
                    bw1.write(geneName.get(i)+"\t"+Math.round(HTPM2[i]));
                    bw1.newLine();
                }
                count++;
            }
            bw.flush();bw.close();
            bw1.flush();bw1.close();
        }
        catch (Exception ex) { 
            ex.printStackTrace();
        }
        
    }
    public void pattern(){
        String countFile="/Users/xujun/Desktop/TEP/TEPOut/geneCount/TEP_count.txt";
        String HTSeqFile="/Users/xujun/Desktop/HTPM.txt";
        String gff="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
        GeneFeature gf = new GeneFeature(gff);
        RowTable r = new RowTable (countFile);
        RowTable t = new RowTable (HTSeqFile);
        double Hcount2 [] = new double [r.getRowNumber()];
        double Hcount3 [] = new double [r.getRowNumber()];
        double Hcount4 [] = new double [r.getRowNumber()];
        List<String> geneName = new ArrayList <>();
        for (int i = 0; i < r.getRowNumber();i++){
            geneName.add(r.getCellAsString(i, 0));
        }
        for(int i=0;i<t.getRowNumber();i++){
            int index = geneName.indexOf(r.getCell(i, 0));
            if(t.getCellAsDouble(i, 1)<=1){
                Hcount2[index]=0;
            }else{
                Hcount2[index]=t.getCellAsDouble(i, 1);
            }
            if(t.getCellAsDouble(i, 2)<=1){
                Hcount3[index]=0;
            }else{
                Hcount3[index]=t.getCellAsDouble(i, 2);
            }
            if(t.getCellAsDouble(i, 3)<=1){
                Hcount4[index]=0;
            }else{
                Hcount4[index]=t.getCellAsDouble(i, 3);
            }
//            Hcount2[index]=t.getCellAsDouble(i, 1);
//            Hcount3[index]=t.getCellAsDouble(i, 2);
//            Hcount4[index]=t.getCellAsDouble(i, 3);
        }
        try{
            BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/HTPMPatern.txt");
            bw.write("\t"+"PLATE-seq"+"\t"+"Lexogen"+"\t"+"SiPAS");bw.newLine();
            for(int i=0;i<Hcount2.length;i++){
                bw.write(geneName.get(i)+"\t"+Hcount2[i]+"\t"+Hcount3[i]+"\t"+Hcount4[i]);
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch (Exception ex) { 
            ex.printStackTrace();
        }
        
        
    }
    public void CoDectededGene(){
        String countTable="/Users/xujun/Desktop/HTPM.txt";
        int count23=0;int count24=0;int count34=0;int countall=0;
        int count2=0;int count3=0;int count4=0;
        RowTable rt=new RowTable(countTable);
        for(int i=0;i<rt.getRowNumber();i++){
            if(rt.getCellAsDouble(i,0)>=1 && rt.getCellAsDouble(i,1)>=1){
                count23++;
            }
            if(rt.getCellAsDouble(i,0)>=1){
                count2++;
            }
            if(rt.getCellAsDouble(i,1)>=1){
                count3++;
            }
            if(rt.getCellAsDouble(i,2)>=1){
                count4++;
            }
            if(rt.getCellAsDouble(i,0)>=1 && rt.getCellAsDouble(i,2)>=1){
                count24++;
            }
            if(rt.getCellAsDouble(i,1)>=1 && rt.getCellAsDouble(i,2)>=1){
                count34++;
            }
            if(rt.getCellAsDouble(i,0)>=1 && rt.getCellAsDouble(i,1)>=1 && rt.getCellAsDouble(i,2)>=1){
                countall++;
            }
        }
        System.out.println(count23);
        System.out.println(count24);
        System.out.println(count34);
        System.out.println(countall);
        System.out.println("***************");
        System.out.println(count2);
        System.out.println(count3);
        System.out.println(count4);
    }
    public void countCompare(){
        String countTable="/Users/xujun/Desktop/TEP/TEPOut/countResult.txt";
        int count=0;
        RowTable rt=new RowTable(countTable);
        for(int i=0;i<rt.getRowNumber();i++){
            if(rt.getCellAsInteger(i,9)!=0 && rt.getCellAsInteger(i,20)!=0){
                count++;
            }
        }
        System.out.println(count);
    }
public void LibraryCompare(){
    String geneName=null;
    String geneName1=null;
    String gffFile="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
    xujun.analysis.rnaseq.GeneFeature gf=new xujun.analysis.rnaseq.GeneFeature(gffFile);
    try{
        BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/difGene.txt");
        BufferedWriter bw1 = IOUtils.getTextWriter("/Users/xujun/Desktop/allGene.txt");
        BufferedWriter bw2 = IOUtils.getTextWriter("/Users/xujun/Desktop/nonProteinGene.txt");
        String countFile="/Users/xujun/library.txt";
        RowTable<String> r =new RowTable<>(countFile);
        int countSame=0;int countDif=0;int countDifS=0;int countdif=0;int count=0;
        int countAll=0;
        for(int i=0;i<r.getRowNumber();i++){
            geneName1=r.getCellAsString(i, 3);
            int geneIndex1=gf.getGeneIndex(geneName1);            
            if(gf.getGeneBiotype(geneIndex1).equals("protein_coding") ){
                for(int j=0;j<gf.genes[geneIndex1].ts.size();j++){
                    if(!gf.genes[geneIndex1].ts.get(j).utr3List.isEmpty()){
                        break;
                    }else{
                        if(j<gf.genes[geneIndex1].ts.size()-1){
                            continue;
                        }else{
                            if(r.getCellAsInteger(i, 0)!=0|| r.getCellAsInteger(i, 1)!=0){
                                if(abs(r.getCellAsInteger(i, 0)-r.getCellAsInteger(i, 1))>=10){
                                    System.out.println(geneName1);
                                    count++;
                                }
                            }
                            countAll++;
                        } 
                    }
                }
//                bw2.write(geneName1+"\t"+gf.getGeneBiotype(geneIndex1)+"\t"+geneIndex1+"\t"+r.getCellAsInteger(i, 0)+"\t"+r.getCellAsInteger(i, 1));bw2.newLine();
//                if(r.getCellAsInteger(i, 0)!=r.getCellAsInteger(i, 1)){
//                    if(abs(r.getCellAsInteger(i, 0)-r.getCellAsInteger(i, 1))>=10){
//                        countdif++;
//                    }
//                    countSame++;
//                    bw2.write(geneName1+"\t"+gf.getGeneBiotype(geneIndex1)+"\t"+geneIndex1+"\t"+r.getCellAsInteger(i, 0)+"\t"+r.getCellAsInteger(i, 1));bw2.newLine();
//                }
            }
            if(abs(r.getCellAsInteger(i, 0)-r.getCellAsInteger(i, 1))>=10){
                countDif++;
                geneName=r.getCellAsString(i, 3);
                int geneIndex=gf.getGeneIndex(geneName);
                bw.write(geneName+"\t"+gf.getGeneBiotype(geneIndex)+"\t"+geneIndex+"\t"+r.getCellAsInteger(i, 0)+"\t"+r.getCellAsInteger(i, 1));bw.newLine();
            }else{
                countDifS++;
            }
            
            bw1.write(geneName1+"\t"+gf.getGeneBiotype(geneIndex1)+"\t"+geneIndex1+"\t"+r.getCellAsInteger(i, 0)+"\t"+r.getCellAsInteger(i, 1));bw1.newLine();
            
        }
        System.out.println(count);
        System.out.println(countAll);
//        System.out.println(countDifS);
        bw.flush();bw.close();
        bw1.flush();bw1.close();
        bw2.flush();bw2.close();
    }
    catch (Exception ex) {
        System.out.println(geneName1+"\t1234");  
        ex.printStackTrace();
               
    }
}
public void CompareTwoCountMethod(){
    String geneName=null;
    String gff="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
    GeneFeature gf = new GeneFeature(gff);
    try{
        String HTSeqFile="/Users/xujun/Desktop/TEP/TEPOut/countResult.txt";
        String Count="/Users/xujun/Desktop/TEP/TEPOut/geneCount/TEP_count.txt";
        RowTable<String> r =new RowTable<>(HTSeqFile);
        RowTable<String> t=new RowTable <> (Count);
        BufferedWriter bw = IOUtils.getTextWriter("/Users/xujun/Desktop/CompareTwoMethod.txt");
        BufferedWriter bw1 = IOUtils.getTextWriter("/Users/xujun/Desktop/CompareTwoMethodIn.txt");
        int[] TEPcount=new int[t.getRowNumber()];
        int [] HTSeqcount = new int [t.getRowNumber()];
        int[] TEPcountIn=new int[t.getRowNumber()];
        int [] HTSeqcountIn = new int [t.getRowNumber()];
        List<String> nameList=new ArrayList<>();
        List<String> incoList=new ArrayList<>();
        for(int i=0;i<r.getRowNumber();i++){
            if(gf.getGeneBiotype(gf.getGeneIndex(r.getCell(i, 0))).equals("protein_coding")){
                nameList.add(r.getCell(i, 0));
                
                int index=nameList.indexOf(r.getCell(i, 0));
                HTSeqcount[index]=r.getCellAsInteger(i, 1);
                int geneIndex=gf.getGeneIndex(r.getCell(i, 0));
                for(int j=0;j<gf.genes[geneIndex].ts.size();j++){
                    if(!gf.genes[geneIndex].ts.get(j).utr3List.isEmpty()){
                        break;
                    }else{
                        if(j<gf.genes[geneIndex].ts.size()-1){
                            continue;
                        }else{
                            incoList.add(r.getCell(i, 0));
                            int indexIn = incoList.indexOf(r.getCell(i, 0));
                            HTSeqcountIn[indexIn]=r.getCellAsInteger(i, 1);
                        } 
                    }
                }
            }    
        }
        for(int i=0;i<t.getRowNumber();i++){
            geneName=t.getCell(i, 0);
            if(nameList.contains(t.getCell(i, 0))){
                int index=nameList.indexOf(t.getCell(i, 0));
                TEPcount[index]=t.getCellAsInteger(i, 5);
    //            System.out.println(t.getCell(i, 0));
            }
            if(incoList.contains(t.getCell(i, 0))){
//                if(t.getCell(i, 0).equals("GRMZM5G881135")){
//                    System.out.print(r.getCell(i, 0));
//                }
                int indexIn1=incoList.indexOf(t.getCell(i, 0));
                TEPcountIn[indexIn1]=t.getCellAsInteger(i, 5);
            }
        }
        for(int i=0;i<nameList.size();i++){
//            System.out.println(countNumber[i]);
            bw.write(nameList.get(i)+"\t"+HTSeqcount[i]+"\t"+TEPcount[i]);
            bw.newLine();
        }
        for(int i=0;i<incoList.size();i++){
            bw1.write(incoList.get(i)+"\t"+HTSeqcountIn[i]+"\t"+TEPcountIn[i]);
            bw1.newLine();
        }
        bw.flush();bw.close();
        bw1.flush();bw1.close();
    }
    catch (Exception ex) {
        System.out.println(geneName+"\t1234");  
        ex.printStackTrace();
               
    }
}
    public void HTSeqCount1(){
        String HTSeqFile="/Users/xujun/Desktop/TEP/TEPOut/HTSeqLibrary/500PL_BC4_Plate4counts.txt";
        List<String> nameList=new ArrayList<>(); 
        String temp=null;String[] tem = null;
        int [] count =new int [49795];
        try{           
            BufferedReader br = IOUtils.getTextReader(HTSeqFile);
            while((temp = br.readLine()) != null){
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(tem[0].startsWith("gene")){
                    if(!nameList.contains(tem[0].split(":")[1])){
                        nameList.add(tem[0].split(":")[1]);
                    }
                    int index=nameList.indexOf(tem[0].split(":")[1]);
                    count[index]+=Integer.parseInt(tem[1]);
                }
            }
            for(int i=0;i<count.length;i++){
                System.out.println(count[i]);
            }
        }
        catch (Exception ex) {
            System.out.println(tem+"\t1234");  
            ex.printStackTrace();

        }   
    }
public void HTSeqCount(){//顺序是HTSeq里面的顺序

    List <String> nameList=new ArrayList<>();
    List<String> fileList=new ArrayList<>();
    String subFqDirS = new File ("/data1/home/junxu/wheat/doubleAlign-Plateseq/HTSeqCount").getAbsolutePath();
    File[] fs = new File(subFqDirS).listFiles();   
    fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
    List<File> fList = Arrays.asList(fs);
    int [][] count=new int[110790][fList.size()];
    fList.stream().forEach(f -> {
        String temp=null;String[] tem = null;
        fileList.add(f.getName());
        try{           
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            while((temp = br.readLine()) != null){
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(tem[0].startsWith("Traes")){
                    if(!nameList.contains(tem[0])){
                        nameList.add(tem[0]);
                    }
                    int index=nameList.indexOf(tem[0]);
                    count[index][fileList.indexOf(f.getName())]+=Integer.parseInt(tem[1]);
                }      
            }
            
        }
        catch (Exception ex) {
            System.out.println(tem[0]+"\t1234");  
            ex.printStackTrace();

        }
    });
    Collections.sort(fileList);
    String outputFileS = new File ("/data1/home/junxu/wheat/doubleAlign-Plateseq/HTSeqCount/double-Plateseq-countResult.txt").getAbsolutePath();
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
                BufferedReader br1 = IOUtils.getTextReader(infile1);
                BufferedReader br2 = IOUtils.getTextReader(infile2);
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
                BufferedReader br1 = IOUtils.getTextReader(infile1);                               
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
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
    public void testL33 () throws IOException {//
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/twice/clean_data";
        String outputDirS = "/Users/xujun/Desktop/RNA_seq/L3-1/clean_data";        
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
            String outfile2 = new File (outputDirS, "333.fq").getAbsolutePath();
            String seq = null;
            String seq1=null; 
            String temp=null;
            String phred=null;
            try {
                BufferedReader br1 = IOUtils.getTextReader(infile1);                              
                BufferedWriter bw = IOUtils.getTextWriter(outfile2);
                List w=new ArrayList();
                List q=new ArrayList();
                //String seq=null;
                while ((temp = br1.readLine()) != null){
                    int a=0;int phred1=0;                 
                    seq=br1.readLine(); 
                    br1.readLine();phred=br1.readLine();
                    for(int i=8;i<seq.length();i++){
                        if(seq.charAt(i)!='T'){
                            a++; 
                        }
                        if(a>2){
//                            if(i > 140){    
//                                seq1=null;
//                                w.add(142);
//                                q.add(0);   
//                            }else{
//                                seq1=seq.substring(i-2,142);
//                                w.add(i-1);
//                                for(int b=i-2;b<142;b++){
//                                    phred1+=(int)phred.charAt(b)-33;
//                                }   
//                                q.add(phred1/seq1.length());
//                                break;
//                            }
                            seq1=seq.substring(i-2);
                            w.add(i-7);
                            for(int b=i-2;b<seq.length();b++){
                                phred1+=(int)phred.charAt(b)-33;
                            }   
                            q.add(phred1/seq1.length());
                            break;
                        }   
                    }
                    if(a<=2){
                        seq1=null;
                        w.add(142);
                        q.add(0);
                    }
//                    bw.write(temp+ "\n");bw.write(seq1+ "\n");bw.write(br1.readLine() + "\n");
//                    bw.write(br1.readLine().substring(seq1.length()) + "\n");     
                }
                br1.close();
                bw.flush();bw.close();
                Collections.shuffle(w); 
                Collections.shuffle(q); 
                int randomSeriesLength = 10000;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength+1);
                List<Integer> randomq = q.subList(0, randomSeriesLength+1);
                System.out.println(p);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      System.out.println(randomq.get(k)+ " "); 
                      if(k==randomSeriesLength){
                           System.out.println("******************************************************** ");
                       } 
                   } else {
                       System.out.print(randomq.get(k)+ " ");
                           
                   }
                }
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
    public void testLexo ()  {//lexogen是不带barcpde的所以单独写了一个
        String barcodeFileS = "/Users/xujun/Desktop/barcodepool.txt";
        String inputDirS ="/Users/xujun/Desktop/RNA_seq/L3-1/clean_data/L3-1";
        String outputDirS = "/Users/xujun/Desktop/RNA_seq/L3-1/clean_data";        
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
            String phred=null;
            try {
                BufferedReader br1 = IOUtils.getTextReader(infile1);                              
                BufferedWriter bw = IOUtils.getTextWriter(outfile2);
                List w=new ArrayList();
                List q=new ArrayList();
                //String seq=null;
                while ((temp = br1.readLine()) != null){
                    int a=0;int phred1=0;                 
                    seq=br1.readLine(); 
                    br1.readLine();phred=br1.readLine();
                    for(int i=0;i<seq.length();i++){
                        if(seq.charAt(i)!='T'){
                            a++; 
                        }
                        if(a>2){
                            if(i > 140){    
                                seq1=null;
                                w.add(142);
                                q.add(0);   
                            }else{
                                seq1=seq.substring(i-2,142);
                                w.add(i-1);
                                for(int b=i-2;b<142;b++){
                                    phred1+=(int)phred.charAt(b)-33;
                                }   
                                q.add(phred1/seq1.length());
                                break;
                            }
                        }   
                    }
//                    bw.write(temp+ "\n");bw.write(seq1+ "\n");bw.write(br1.readLine() + "\n");
//                    bw.write(br1.readLine().substring(seq1.length()) + "\n");     
                }
                br1.close();
                bw.flush();bw.close();
                Collections.shuffle(w); 
                Collections.shuffle(q); 
                int randomSeriesLength = 1000;
                List<Integer> randomSeries = w.subList(0, randomSeriesLength+1);
                List<Integer> randomq = q.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                   if (k % 25 == 0){ 
                      System.out.println(randomq.get(k)+ " "); 
                      if(k==randomSeriesLength){
                           System.out.println("******************************************************** ");
                       } 
                   } else {
                       System.out.print(randomq.get(k)+ " ");
                           
                   }
                }
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
                BufferedReader br1 = IOUtils.getTextReader(inputFile1);   
                BufferedReader br2 = IOUtils.getTextReader(inputFile2);
                BufferedWriter bw = IOUtils.getTextWriter(outputFile);
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
                BufferedReader br1 = IOUtils.getTextReader(inputFile1);   
                BufferedWriter bw = IOUtils.getTextWriter(outputFile);
                
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
                BufferedReader br1 = IOUtils.getTextReader(inputFile1);   
                BufferedReader br2 = IOUtils.getTextReader(inputFile2);
                BufferedReader br3 = IOUtils.getTextReader(inputFile3);
                BufferedReader br4 = IOUtils.getTextReader(inputFile4);
                BufferedWriter bw = IOUtils.getTextWriter(outputFile);
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
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);            
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
             BufferedReader br1 = IOUtils.getTextReader(inputFile1);
             BufferedReader br2 = IOUtils.getTextReader(inputFile2);
             BufferedReader br3 = IOUtils.getTextReader(inputFile3);
             BufferedReader br4 = IOUtils.getTextReader(inputFile4);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);            
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
            BufferedReader br1 = IOUtils.getTextReader(inputfile1);
            BufferedReader br2 = IOUtils.getTextReader(inputfile2);
            BufferedWriter bw = IOUtils.getTextWriter(outputfile); 
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
             BufferedReader br1 = IOUtils.getTextReader(inputfile);
             BufferedWriter bw = IOUtils.getTextWriter(outputfile);            
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
    public static void main(String[] args) throws IOException, FileNotFoundException{
        new DemoSample();
//        new GBS();
    }
}




