/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

/**
 *
 * @author xujun
 */
import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import com.koloboke.collect.map.hash.HashIntObjMap;
import com.koloboke.collect.map.hash.HashIntObjMaps;
import pgl.infra.range.Range;
import pgl.infra.range.RangeValStr;
import pgl.infra.table.RowTable;
//import htsjdk.samtools.BAMFileReader;
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
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;
import rcaller.RCaller;
import rcaller.RCode;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author feilu
 */
public class WheatTEP {
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
    //The directory of gff file
    String rPath=null;
    String gffFile=null;
    GeneFeature.Gene[] genes;
    HashMap<String, RangeValStr> geneNameRangeMap = new HashMap<>();
    String[] geneNames = null;
    Integer[] chrs = null;
    String[][] geneNamesByChr = null;
    HashIntIntMap[] posGeneMaps = null;
    HashMap[] geneRangeMaps=null;
    
    int overhangLength = 150;
    
    int multiMapN = 10;
    
    float mismatchRate = (float)0.1;
    
    int minNMatch = 80;
    
    int[] barcodeLengths = null;
    
    String[] subDirS = {"subFastqs", "sams", "geneCount","countTable"};
    List<String> fqFileSList = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeTaxaMaps = null;
//    HashMap<Integer,Integer> 
    
    List<String>[] taxaLists = null;
    
    List<String> diffValumnMethodList = new ArrayList<String>();
    int geneNumber = 110790;
//    int geneNumber = 110882;

 
    public WheatTEP (String parameterFileS) {
//        this.getTPM("/Users/xujun/Desktop/TEP/TEPOut/sams/PL_BC1_Plate1Aligned.out.sam");
        this.parseParameters(parameterFileS);
//        this.dataValumn();
//        this.starAlignment();
//        this.starAlignmentPair();
//        this.mkPosGeneMap();
//        this.parseFq(); //Remove Ts?
//        this.mkIndexOfReference(); //one-time set, not requried for standard runs
//        this.starAlignmentPair();
//        this.mkGeneCountTableSE();
//        this.mkGeneCountTablePE();
//        this.HTSeqCountSingle();
        this.HTSeqCountDouble();
        this.MergeHTSeq();
//        this.countTPM();
        this.expGene();
    }
    public void expGene(){
//        RowTable rt = new RowTable(new File (this.outputDirS,subDirS[3]).getAbsolutePath()+"/countResult.txt");
//        RowTable rt = new RowTable("/data1/home/junxu/analysis0215/7/countResult.txt");
        RowTable rt = new RowTable("/Users/xujun/Desktop/eQTL/total/Homology/190425/countResult.txt");
        int isExp [][] = new int [this.geneNumber][rt.getColumnNumber()];
//        List<String> fileList=new ArrayList<>();
//        List<String> geneList=new ArrayList<>();
//        System.out.println(rt.getColumnNumber()+"/t"+rt.getRowNumber());
        for(int i=0; i< rt.getRowNumber();i++){
            for(int j=1; j < rt.getColumnNumber()-1;j++){
                if(rt.getCellAsInteger(i, j)>1){
                    isExp[i][j]=1;
                }else{
                    isExp[i][j]=0;
                }   
            }
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(new File("/Users/xujun/Desktop/eQTL/total/Homology/190425").getAbsolutePath()+"/expGene.txt");
            for(int i =0;i<isExp.length;i++){
                for(int j=0;j<isExp[i].length;j++){
                    bw.write(isExp[i][j]+"\t");
                }
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    public void countTPM(){
        RCaller caller = new RCaller();
        caller.setRscriptExecutable(this.rPath);
        RCode rCode = new RCode();
        rCode.addRCode(this.getTPM(new File (this.outputDirS,subDirS[3]).getAbsolutePath()+"/countResult.txt"));
        caller.setRCode(rCode);
        caller.runOnly();
    }
    public String getTPM(String f) {
        StringBuilder sb = new StringBuilder ();
        sb.append("x=read.csv(\""+f+"\",sep = \"\\t\")"+"\n");
        sb.append("countToTpm <- function(counts, effLen)"+"\n").append("{"+"\n");
        sb.append("\t"+"rate <- log(counts) - log(effLen)"+"\n");
        sb.append("\t"+"denom <- log(sum(exp(rate)))"+"\n");
        sb.append("\t"+"exp(rate - denom + log(1e6))"+"\n");
        sb.append("}"+"\n");
        sb.append("countTable<-data.frame(count=x[,6:101],length=x$End-x$Start)"+"\n");
        sb.append("for(i in 1:95){"+"\n");
        sb.append("\t"+"countTable[,i]<-with(countTable,countToTpm(countTable[,i],length))"+"\n");
        sb.append("}"+"\n");
        sb.append("write.table(countTable,file=\""+f.replace(".txt", "TPM.txt")+"\",append = F,quote = F,sep = \"\\t\",eol = \"\\n\",row.names = F,col.names = T)");
        String statement = sb.toString();
        System.out.println(statement);
        return statement;
    }
    public void dataValumn() {
        String subFqDirS = new File (this.outputDirS).getAbsolutePath();
//        String subFqDirS = "xujun/TEP/TEPOut/subFastqs";
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
//            if (fs[i].length() < 1100000000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fList.add(fs[i]);
        }
//        String command =null;
        int numCores = Runtime.getRuntime().availableProcessors();
        fList.stream().forEach(f -> {
                try { 
                    File dir = new File(new File (this.outputDirS,subDirS[0]).getAbsolutePath());
//                    for( int i = 1 ;i<=12;i++){
                    if (f.length() < 1810000000){
                        StringBuilder sb = new StringBuilder();
                        sb.append("cp "+f+" "+dir);
                        String command = sb.toString();
                        System.out.println(command); 
                        Process p=Runtime.getRuntime().exec(command);
                        p.waitFor();
                    }else{
                        StringBuilder sb = new StringBuilder();
//                        sb.append("/data1/home/junxu/wheat/seqtk/seqtk sample -s100 ").append(f+" ").append(1000000*i+" ").append(" >> ").append(" "+1000000*i+f.getName());
                        sb.append("/data1/home/junxu/wheat/seqtk/seqtk sample ").append(f+" ").append(5000000+" ").append(" >> ").append(" "+5000000+f.getName());
                        String command = sb.toString();
                        System.out.println(command);  
                        String []cmdarry ={"/bin/bash","-c",command};
                        Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                        p.waitFor();
                    }    
//                    } 
//                    for( int i = 1 ;i<=20;i++){
//                        StringBuilder sb = new StringBuilder();
////                        sb.append("/data1/home/junxu/wheat/seqtk/seqtk sample -s100 ").append(f+" ").append(1000000*i+" ").append(" >> ").append(" "+1000000*i+f.getName());
//                        sb.append("/data1/home/junxu/wheat/seqtk/seqtk sample ").append(f+" ").append(5000000+" ").append(" >> ").append(" "+5000000+f.getName().replace(".fq", "-"+i+".fq"));
//                        String command = sb.toString();
//                        System.out.println(command);  
//                        String []cmdarry ={"/bin/bash","-c",command};
//                        Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
//                        p.waitFor();
////                    }    
//                    } 
                }
                catch (Exception e) {
                    e.printStackTrace();
                }            
        });
    }
    public void HTSeqCountSingle(){
        String inputDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {  
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s no ");
            sb.append(f);
            sb.append(" "+"/data1/home/junxu/wheat/rightchangewheat.gtf").append(" >> ");
//            sb.append(" "+this.gffFile).append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File (this.outputDirS,subDirS[2]).getAbsolutePath());
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
    public void HTSeqCountDouble(){
        String inputDirS = new File (this.outputDirS,subDirS[1]).getAbsolutePath();
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {  
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s reverse ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat_v1.1_Lulab.gtf").append(" >> ");
//            sb.append(" "+this.gffFile).append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File (this.outputDirS,subDirS[2]).getAbsolutePath());
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
    public void MergeHTSeq(){
        HashMap<Integer,String> namePos = new HashMap();
        List<String> fileList=new ArrayList<>();
//        String subGeneCountDirS = new File (this.outputDirS,subDirS[2]).getAbsolutePath();
        String subGeneCountDirS = new File ("/data1/home/junxu/wheat/doubleAll/UniqPE/count").getAbsolutePath();
        File[] fs = new File(subGeneCountDirS).listFiles();   
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fList = Arrays.asList(fs);
        for(int i=0;i < fList.size();i++){
            fileList.add(fList.get(i).getName().replace("Count.txt", ""));
        }
        Collections.sort(fileList);
        int[][] geneCount = new int[this.geneNumber][]; // chr, geneByChr, taxon
        for (int i = 0; i < geneNumber; i++) {
            geneCount[i] = new int[fList.size()];   
        }
        fList.stream().forEach(f -> {
            String temp=null;String[] tem = null;
            int contG=0;
            try{           
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                int a=0;
                while((temp = br.readLine()) != null){
                    a++;
                    List<String> tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
    //                if(tem[0].startsWith("ERCC")){
                      if(tem[0].startsWith("TraesCS")){
//                        if(!nameList.contains(tem[0])){
//                            nameList.add(tem[0]);
//                        }
                        namePos.put(contG, tem[0]);
//                        int index=nameList.indexOf(tem[0]);
//                        geneCount[index][fileList.indexOf(f.getName().replace(".txt", ""))]=Integer.parseInt(tem[1]);
                        geneCount[contG][fileList.indexOf(f.getName().replace("Count.txt", ""))]=Integer.parseInt(tem[1]);
                    }  
                    contG++;
                    if(a>110790){
                        break;
                    }
                }
            }
            catch (Exception ex) {
                System.out.println(tem[0]+"\t1234");  
                ex.printStackTrace();
            }
        });
//        String outputFileS = new File (this.outputDirS,subDirS[3]).getAbsolutePath()+"/countResult.txt";
        String outputFileS = new File ("/data1/home/junxu/wheat/doubleAll/UniqPE/count/countResult.txt").getAbsolutePath();
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
//            sb.append("Gene\tChr\tStart\tEnd\tStrand\t");
            sb.append("Gene\t");
            for(int i=0;i<fileList.size();i++){            
                sb.append(fileList.get(i).replace(".txt", "")+"\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for(int i=0;i<geneCount.length;i++){
                sb = new StringBuilder();  
                for(int j=0;j<fileList.size();j++){
                    if(j==0){
//                        RangeValStr r = geneNameRangeMap.get(nameList.get(i));
//                        sb.append(nameList.get(i)+"\t");//.append(r.chr).append("\t");
                        sb.append(namePos.get(i)+"\t");
//                        sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str+"\t");
                    }
                    sb.append(geneCount[i][j]+"\t");           
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
    public void MergeHTSeqAndRPKM(){
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
        int[][] geneCount = new int[this.geneNumber][]; // chr, geneByChr, taxon
        double [][] TPM=new double[this.geneNumber][];
        int [] allReadsNumber = new int [fList.size()];
        for (int i = 0; i < geneNumber; i++) {
            geneCount[i] = new int[fList.size()];   
            TPM[i] = new double[fList.size()]; 
        }
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
        for(int i=0;i<fList.size();i++){
            for(int j=0;j<this.geneNumber;j++){
                RangeValStr v = geneNameRangeMap.get(nameList.get(j));
                TPM[j][i]=geneCount[j][i]*1000000/allReadsNumber[i]/(v.end-v.start+1)*1000;
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
                    sb.append(geneCount[i][j]+"\t"+TPM[i][j]+"\t");           
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
     private void mkGeneCountTable () {
//        String inputDirS = new File(this.outputDirS, this.subDirS[3]).getAbsolutePath();
        String inputDirS = new File("/data1/home/junxu/wheat/singleAlign-SiPASR2/sams").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace("Aligned.out.sam", ""));
            
        }
        Collections.sort(diffValumnMethodList);
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
        int [][][] MapQ=new int[chrs.length][][];
        int [][][] quality=new int[chrs.length][][];
        int [][][] mappedLength=new int[chrs.length][][];
        double [][][] TPM=new double[chrs.length][][];
        double [][][] effLen=new double [diffValumnMethodList.size()][chrs.length][];
        ArrayList<Integer> [][]readsLengthList = new ArrayList[diffValumnMethodList.size()][];
        for(int i=0;i<diffValumnMethodList.size();i++){
            readsLengthList[i]=new ArrayList[geneNames.length];
            for(int k=0;k<readsLengthList[i].length;k++){
                readsLengthList[i][k]=new ArrayList ();//在每种方法之下每一个基因下面都建立一个List用来储存不同的reads长度
            }
            for(int j=0;j<chrs.length;j++){
                effLen[i][j]=new double[geneNamesByChr[j].length];
            }
        }
        for (int i = 0; i < geneCount.length; i++) {
            geneCount[i] = new int[geneNamesByChr[i].length][];   
            TPM[i]= new double[geneNamesByChr[i].length][];
            MapQ[i] = new int [geneNamesByChr[i].length][];
            quality[i] = new int [geneNamesByChr[i].length][];
            mappedLength[i] = new int [geneNamesByChr[i].length][];
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[diffValumnMethodList.size()]; 
                TPM[i][j] = new double[diffValumnMethodList.size()];  
                MapQ[i][j] = new int[diffValumnMethodList.size()];
                quality[i][j] = new int[diffValumnMethodList.size()];
                mappedLength[i][j] = new int[diffValumnMethodList.size()];
            }
        }
        int [] expreNumber = new int [diffValumnMethodList.size()];
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                List<Range> geneRange=new ArrayList(); 
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
                int taxonIndex = Collections.binarySearch(diffValumnMethodList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;String Q=null;
                String cigar=null;int flag=0;
                while ((temp = br.readLine()) != null) {
                    int effLenPerReads =0; int mapq=0;
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));//应该是返回这个染色体所在的角标
                    if (chrIndex < 0) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];Q=tem[10];
                    mapq=Integer.valueOf(tem[4]);
                    flag=Integer.parseInt(tem[1]);
                    int cIndex = 0; int phred=0;
                    int cLength = 0;int skip=0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            if(c== 'M'){
                                effLenPerReads += Integer.parseInt(cigar.substring(cIndex, i));
                                for(int a=skip;a<skip+Integer.parseInt(cigar.substring(cIndex, i));a++){
                                    phred += (int)Q.charAt(a)-33;
                                }
                            }
                            if(c=='M'||c=='S'||c=='I'){
                                skip+=Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    if((flag & 16) ==16){//负链上面 好想PE里面所有的flag值&16都等于16或者是0 但是有一个数&16是90 再&上16的话就成了0 decode sam flag 
                        int index1=posGeneMaps[chrIndex*2+1].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2+1].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2+1].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                            mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                            MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                            quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        int index1=posGeneMaps[chrIndex*2].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                            mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                            MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                            quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }         
                }
            }
            catch (Exception e) {
                System.out.print(geneCount.length);
                e.printStackTrace();
            }
        });
        //standerdize
        double[][][]rate=new double[this.diffValumnMethodList.size()][chrs.length][];
        double [] denom=new double[this.diffValumnMethodList.size()];
        double [] allRate=new double [this.diffValumnMethodList.size()];
        double effLength=0;
        for(int i=0;i<this.diffValumnMethodList.size();i++){
            for(int j=0;j<chrs.length;j++){
                rate[i][j] = new double [geneNamesByChr[j].length];
                for(int k=0;k<geneNamesByChr[j].length;k++){        
                    if(geneCount[j][k][i]==0){}
                    else{
                        RangeValStr v = geneNameRangeMap.get(geneNamesByChr[j][k]);
                        effLength=v.end-v.start+1-effLen[i][j][k]/geneCount[j][k][i];
                        if(effLength<=0){
                            effLength = v.end-v.start+1;
                        }
                        rate[i][j][k]=Math.log(geneCount[j][k][i])-Math.log(effLength);
                        allRate[i]+=Math.exp(rate[i][j][k]);
                    }                    
                }
            }
            denom[i]=Math.log(allRate[i]);
        }
        for(int i=0;i<chrs.length;i++){
            for(int j=0;j<geneNamesByChr[i].length;j++){
                for(int k=0;k<this.diffValumnMethodList.size();k++){
                    if(geneCount[i][j][k]==0){TPM[i][j][k]=0;}
                    else{
                        TPM[i][j][k]=Math.exp(rate[k][i][j]-denom[k]+Math.log(1e6));
                        mappedLength[i][j][k]=mappedLength[i][j][k]/geneCount[i][j][k];
                        MapQ[i][j][k]=MapQ[i][j][k]/geneCount[i][j][k];
                        quality[i][j][k]=quality[i][j][k]/geneCount[i][j][k];
                    }                   
                }
            }
        } 
//        String outfileDirS = new File (this.outputDirS, this.subDirS[2]).getAbsolutePath();
        String outfileDirS = new File ("/data1/home/junxu/wheat/singleAlign-SiPASR2/geneCount").getAbsolutePath();
        String outputFileS = new File (outfileDirS, "TPM_diffValumnTEP_count.txt").getAbsolutePath();
//        String outputFileS1 = new File (outfileDirS, "readLemgth.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tStart\tEnd\tStrand\t");
            for (int i = 0; i < this.diffValumnMethodList.size(); i++) {
                sb.append(this.diffValumnMethodList.get(i)+"\t"+"readLength"+"\t"+"MapQ"+"\t"+"quality"+"\t");//.append(this.diffValumnMethodList.get(i)+":TPM");+"expre"+"\t"
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outputFileS1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append("\t");
                    sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str+"\t");
                    for (int k = 0; k < this.diffValumnMethodList.size(); k++) {
                        sb.append(TPM[i][j][k]+"\t"+mappedLength[i][j][k]+"\t"+MapQ[i][j][k]+"\t"+quality[i][j][k]+"\t");//.append("\t").append(geneCount[i][j][k])
//                        int expre =0;
//                        if(geneCount[i][j][k]==0){
//                            expre=0;
//                        }else{
//                            expre = 1;
//                            expreNumber[k]++;
//                        }
//                        sb.append(expre+"\t");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
//            for (int k = 0;  k<diffValumnMethodList.size(); k++) {
//                
//                   bw1.write(Integer.valueOf(expreNumber[k]));bw1.newLine();              
//                
//            }
            bw.flush();//bw1.flush();
            bw.close();//bw1.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
     private void mkGeneCountTableSE () {
//        String inputDirS = new File(this.outputDirS, this.subDirS[1]).getAbsolutePath();
        String inputDirS = new File("/data1/home/junxu/analysis0215/TruseqResult/result-all-withoutERCC/sams").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace("Aligned.out.sam", "")); 
        }
        Collections.sort(diffValumnMethodList);
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
        for (int i = 0; i < geneCount.length; i++) {
            geneCount[i] = new int[geneNamesByChr[i].length][];   
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[diffValumnMethodList.size()]; 
            }
        }
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                List<Range> geneRange=new ArrayList(); 
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
//                int taxonIndex = Collections.binarySearch(methodList, taxon);
                int taxonIndex = Collections.binarySearch(diffValumnMethodList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;String Q=null;
                String cigar=null;int flag=0;
                while ((temp = br.readLine()) != null) {
                    int effLenPerReads =0; int mapq=0;
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));//应该是返回这个染色体所在的角标
                    if (chrIndex < 0) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];Q=tem[10];
                    mapq=Integer.valueOf(tem[4]);
                    flag=Integer.parseInt(tem[1]);
                    int cIndex = 0; int phred=0;
                    int cLength = 0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    if(flag==16){//负链上面 
                        int index1=posGeneMaps[chrIndex*2+1].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2+1].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2+1].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        if(flag==0){
                            int index1=posGeneMaps[chrIndex*2].get(startPos);
                            if (index1 < 0) continue;
                            int index2=posGeneMaps[chrIndex*2].get(endPos);
                            if(index1 == index2){
                                geneIndex=index1;
                                geneRange=(List<Range>) geneRangeMaps[chrIndex*2].get(geneIndex);
                                if(geneRange.size()==1){
                                    geneCount[chrIndex][geneIndex][taxonIndex]++;
                                }else{
                                    for(int n=0;n<geneRange.size();n++){
                                        if(geneRange.get(n).isContain(chrIndex, startPos)){
                                            if(geneRange.get(n).isContain(chrIndex, endPos)){
                                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }         
                }
            }
            catch (Exception e) {
                System.out.print(geneCount.length);
                e.printStackTrace();
            }
        });
//        String outfileDirS = new File (this.outputDirS, this.subDirS[2]).getAbsolutePath();
        String outfileDirS = new File (this.outputDirS,subDirS[3]).getAbsolutePath();
        String outputFileS = new File (outfileDirS, "raw_count.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tStart\tEnd\tStrand\t");
            for (int i = 0; i < this.diffValumnMethodList.size(); i++) {
                sb.append(this.diffValumnMethodList.get(i)+"\t");//.append(this.diffValumnMethodList.get(i)+":TPM");+"expre"+"\t" "readLength"+"\t"+"MapQ"+"\t"+"quality"+"\t"
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outputFileS1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append("\t");
                    sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str+"\t");
                    for (int k = 0; k < this.diffValumnMethodList.size(); k++) {
                        sb.append(geneCount[i][j][k]+"\t");//.append("\t").append(geneCount[i][j][k]) +mappedLength[i][j][k]+"\t"+MapQ[i][j][k]+"\t"+quality[i][j][k]+"\t"
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();//bw1.flush();
            bw.close();//bw1.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
     private void mkGeneCountTablePE () {
//        String inputDirS = new File(this.outputDirS, this.subDirS[1]).getAbsolutePath();
        String inputDirS = new File("/data1/home/junxu/wheat/doubleAlign-SiPAS/sams").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace("Aligned.out.sam", "")); 
        }
        Collections.sort(diffValumnMethodList);
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
        for (int i = 0; i < geneCount.length; i++) {
            geneCount[i] = new int[geneNamesByChr[i].length][];   
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[diffValumnMethodList.size()]; 
            }
        }
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                List<Range> geneRange=new ArrayList(); 
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
//                int taxonIndex = Collections.binarySearch(methodList, taxon);
                int taxonIndex = Collections.binarySearch(diffValumnMethodList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;String Q=null;
                String cigar=null;int flag=0;
                while ((temp = br.readLine()) != null) {
                    int effLenPerReads =0; int mapq=0;
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));//应该是返回这个染色体所在的角标
                    if (chrIndex < 0) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];Q=tem[10];
                    mapq=Integer.valueOf(tem[4]);
                    flag=Integer.parseInt(tem[1]);
                    int cIndex = 0; int phred=0;
                    int cLength = 0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    if(flag==83 || flag==147){//负链上面 
                        int index1=posGeneMaps[chrIndex*2+1].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2+1].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2+1].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        if(flag==163|| flag==99 || flag==137){
                            int index1=posGeneMaps[chrIndex*2].get(startPos);
                            if (index1 < 0) continue;
                            int index2=posGeneMaps[chrIndex*2].get(endPos);
                            if(index1 == index2){
                                geneIndex=index1;
                                geneRange=(List<Range>) geneRangeMaps[chrIndex*2].get(geneIndex);
                                if(geneRange.size()==1){
                                    geneCount[chrIndex][geneIndex][taxonIndex]++;
                                }else{
                                    for(int n=0;n<geneRange.size();n++){
                                        if(geneRange.get(n).isContain(chrIndex, startPos)){
                                            if(geneRange.get(n).isContain(chrIndex, endPos)){
                                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }         
                }
            }
            catch (Exception e) {
                System.out.print(geneCount.length);
                e.printStackTrace();
            }
        });
//        String outfileDirS = new File (this.outputDirS, this.subDirS[2]).getAbsolutePath();
        String outfileDirS = new File (this.outputDirS,subDirS[3]).getAbsolutePath();
        String outputFileS = new File (outfileDirS, "raw_count.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tStart\tEnd\tStrand\t");
            for (int i = 0; i < this.diffValumnMethodList.size(); i++) {
                sb.append(this.diffValumnMethodList.get(i)+"\t");//.append(this.diffValumnMethodList.get(i)+":TPM");+"expre"+"\t" "readLength"+"\t"+"MapQ"+"\t"+"quality"+"\t"
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outputFileS1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append("\t");
                    sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str+"\t");
                    for (int k = 0; k < this.diffValumnMethodList.size(); k++) {
                        sb.append(geneCount[i][j][k]+"\t");//.append("\t").append(geneCount[i][j][k]) +mappedLength[i][j][k]+"\t"+MapQ[i][j][k]+"\t"+quality[i][j][k]+"\t"
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();//bw1.flush();
            bw.close();//bw1.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    private void mkGeneCountTableAndTPM () {
//        String inputDirS = new File(this.outputDirS, this.subDirS[1]).getAbsolutePath();
        String inputDirS = new File("/data1/home/junxu/wheat/TEPout/diffValumnSams/some").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace("Aligned.out.sam", "")); 
        }
        Collections.sort(diffValumnMethodList);
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
        int [][][] MapQ=new int[chrs.length][][];
        int [][][] quality=new int[chrs.length][][];
        int [][][] mappedLength=new int[chrs.length][][];
        double [][][] TPM=new double[chrs.length][][];
        double [][][] effLen=new double [diffValumnMethodList.size()][chrs.length][];
        ArrayList<Integer> [][]readsLengthList = new ArrayList[diffValumnMethodList.size()][];
        for(int i=0;i<diffValumnMethodList.size();i++){
            readsLengthList[i]=new ArrayList[geneNames.length];
            for(int k=0;k<readsLengthList[i].length;k++){
                readsLengthList[i][k]=new ArrayList ();//在每种方法之下每一个基因下面都建立一个List用来储存不同的reads长度
            }
            for(int j=0;j<chrs.length;j++){
                effLen[i][j]=new double[geneNamesByChr[j].length];
            }
        }
        for (int i = 0; i < geneCount.length; i++) {
            geneCount[i] = new int[geneNamesByChr[i].length][];   
            TPM[i]= new double[geneNamesByChr[i].length][];
            MapQ[i] = new int [geneNamesByChr[i].length][];
            quality[i] = new int [geneNamesByChr[i].length][];
            mappedLength[i] = new int [geneNamesByChr[i].length][];
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[diffValumnMethodList.size()]; 
                TPM[i][j] = new double[diffValumnMethodList.size()];  
                MapQ[i][j] = new int[diffValumnMethodList.size()];
                quality[i][j] = new int[diffValumnMethodList.size()];
                mappedLength[i][j] = new int[diffValumnMethodList.size()];
            }
        }
        int [] expreNumber = new int [diffValumnMethodList.size()];
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                List<Range> geneRange=new ArrayList(); 
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
//                int taxonIndex = Collections.binarySearch(methodList, taxon);
                int taxonIndex = Collections.binarySearch(diffValumnMethodList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;String Q=null;
                String cigar=null;int flag=0;
                while ((temp = br.readLine()) != null) {
                    int effLenPerReads =0; int mapq=0;
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));//应该是返回这个染色体所在的角标
                    if (chrIndex < 0) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];Q=tem[10];
                    mapq=Integer.valueOf(tem[4]);
                    flag=Integer.parseInt(tem[1]);
                    int cIndex = 0; int phred=0;
                    int cLength = 0;int skip=0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            if(c== 'M'){
                                effLenPerReads += Integer.parseInt(cigar.substring(cIndex, i));
                                for(int a=skip;a<skip+Integer.parseInt(cigar.substring(cIndex, i));a++){
                                    phred += (int)Q.charAt(a)-33;
                                }
                            }
                            if(c=='M'||c=='S'||c=='I'){
                                skip+=Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    if(flag==16){//负链上面 
                        int index1=posGeneMaps[chrIndex*2+1].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2+1].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2+1].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                            mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                            MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                            quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        if(flag==0){
                            int index1=posGeneMaps[chrIndex*2].get(startPos);
                            if (index1 < 0) continue;
                            int index2=posGeneMaps[chrIndex*2].get(endPos);
                            if(index1 == index2){
                                geneIndex=index1;
                                geneRange=(List<Range>) geneRangeMaps[chrIndex*2].get(geneIndex);
                                if(geneRange.size()==1){
                                    geneCount[chrIndex][geneIndex][taxonIndex]++;
                                    effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                    mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                    MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                    quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                                }else{
                                    for(int n=0;n<geneRange.size();n++){
                                        if(geneRange.get(n).isContain(chrIndex, startPos)){
                                            if(geneRange.get(n).isContain(chrIndex, endPos)){
                                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                                effLen[taxonIndex][chrIndex][geneIndex]+=effLenPerReads;
                                                mappedLength[chrIndex][geneIndex][taxonIndex]+=effLenPerReads;
                                                MapQ[chrIndex][geneIndex][taxonIndex]+=mapq;//MapQ值
                                                quality[chrIndex][geneIndex][taxonIndex]+=phred/effLenPerReads;//read的平均质量值
                                                break;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        
                    }         
                }
            }
            catch (Exception e) {
                System.out.print(geneCount.length);
                e.printStackTrace();
            }
        });
        //standerdize
        double[][][]rate=new double[this.diffValumnMethodList.size()][chrs.length][];
        double [] denom=new double[this.diffValumnMethodList.size()];
        double [] allRate=new double [this.diffValumnMethodList.size()];
        double effLength=0;
        for(int i=0;i<this.diffValumnMethodList.size();i++){
            for(int j=0;j<chrs.length;j++){
                rate[i][j] = new double [geneNamesByChr[j].length];
                for(int k=0;k<geneNamesByChr[j].length;k++){        
                    if(geneCount[j][k][i]==0){}
                    else{
                        RangeValStr v = geneNameRangeMap.get(geneNamesByChr[j][k]);
                        effLength=v.end-v.start+1-effLen[i][j][k]/geneCount[j][k][i];
                        if(effLength<=0){
                            effLength = v.end-v.start+1;
                        }
                        rate[i][j][k]=Math.log(geneCount[j][k][i])-Math.log(effLength);
                        allRate[i]+=Math.exp(rate[i][j][k]);
                    }                    
                }
            }
            denom[i]=Math.log(allRate[i]);
        }
        for(int i=0;i<chrs.length;i++){
            for(int j=0;j<geneNamesByChr[i].length;j++){
                for(int k=0;k<this.diffValumnMethodList.size();k++){
                    if(geneCount[i][j][k]==0){TPM[i][j][k]=0;}
                    else{
                        TPM[i][j][k]=Math.exp(rate[k][i][j]-denom[k]+Math.log(1e6));
                        mappedLength[i][j][k]=mappedLength[i][j][k]/geneCount[i][j][k];
                        MapQ[i][j][k]=MapQ[i][j][k]/geneCount[i][j][k];
                        quality[i][j][k]=quality[i][j][k]/geneCount[i][j][k];
                    }                   
                }
            }
        } 
//        String outfileDirS = new File (this.outputDirS, this.subDirS[2]).getAbsolutePath();
        String outfileDirS = new File (this.outputDirS,subDirS[2]).getAbsolutePath();
        String outputFileS = new File (outfileDirS, "rawandTPM-diffValumnTEP_count.txt").getAbsolutePath();
//        String outputFileS1 = new File (outfileDirS, "readLemgth.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tStart\tEnd\tStrand\t");
            for (int i = 0; i < this.diffValumnMethodList.size(); i++) {
                sb.append(this.diffValumnMethodList.get(i)+"\t"+"TPM"+"\t");//.append(this.diffValumnMethodList.get(i)+":TPM");+"expre"+"\t" "readLength"+"\t"+"MapQ"+"\t"+"quality"+"\t"
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
//            BufferedWriter bw1 = IOUtils.getTextWriter(outputFileS1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append("\t");
                    sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str+"\t");
                    for (int k = 0; k < this.diffValumnMethodList.size(); k++) {
                        sb.append(geneCount[i][j][k]+"\t"+TPM[i][j][k]+"\t");//.append("\t").append(geneCount[i][j][k]) +mappedLength[i][j][k]+"\t"+MapQ[i][j][k]+"\t"+quality[i][j][k]+"\t"
//                        int expre =0;
//                        if(geneCount[i][j][k]==0){
//                            expre=0;
//                        }else{
//                            expre = 1;
//                            expreNumber[k]++;
//                        }
//                        sb.append(expre+"\t");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
//            for (int k = 0;  k<diffValumnMethodList.size(); k++) {
//                
//                   bw1.write(Integer.valueOf(expreNumber[k]));bw1.newLine();              
//                
//            }
            bw.flush();//bw1.flush();
            bw.close();//bw1.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkPosGeneMap () {
        int chr = -1 ;
         String gffFile="/data1/home/junxu/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3";
//        String gffFile=this.gffFile;
        GeneFeature gf=new GeneFeature(gffFile); 
        String geneNameS=null;int gfIndex=0;
        try{
            BufferedReader br = IOUtils.getTextReader("/data1/home/junxu/wheat/rightchangewheat.gtf");
//            BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/wheat/rightchangewheat.gtf");
            String temp = null;
            Set<String> geneSet = new HashSet<String>();
            Set<Integer> chrSet = new HashSet<>();
            String[] tem = null;
            String geneName = null;
            HashMap<String, Integer> geneChrMap = new HashMap();
            HashMap<String, Integer> geneMinMap = new HashMap();
            HashMap<String, Integer> geneMaxMap = new HashMap();
            HashMap<String, Byte> geneStrandMap = new HashMap();
            while ((temp = br.readLine()) != null) {
                List<String> tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (!tem[2].startsWith("exon")) continue;
                String[] te = tem[8].split(";");
                geneName=te[1].split("\"")[1];
                if (!geneSet.contains(geneName)) {
                    geneMinMap.put(geneName, Integer.MAX_VALUE);
                    geneMaxMap.put(geneName, Integer.MIN_VALUE);
                }
                geneSet.add(geneName);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                if (geneMinMap.get(geneName) > min) geneMinMap.put(geneName, min);
                if (geneMaxMap.get(geneName) < max) geneMaxMap.put(geneName, max);
                chr = Integer.parseInt(tem[0]);            
                chrSet.add(chr);
                geneChrMap.put(geneName, chr);
                if (tem[6].startsWith("-")) geneStrandMap.put(geneName, (byte)0);
                else geneStrandMap.put(geneName, (byte)1);                              
            }
            geneNames = geneSet.toArray(new String[geneSet.size()]);
            Arrays.sort(geneNames);
            for (int i = 0; i < geneNames.length; i++) {
                RangeValStr r = new RangeValStr (geneChrMap.get(geneNames[i]), geneMinMap.get(geneNames[i]), geneMaxMap.get(geneNames[i]), 0, geneStrandMap.get(geneNames[i]));
                geneNameRangeMap.put(geneNames[i], r);
            }
            chrs = chrSet.toArray(new Integer[chrSet.size()]);
            Arrays.sort(chrs);
            geneNamesByChr = new String[chrs.length][];
            geneRangeMaps=new HashMap[chrs.length*2];
            posGeneMaps = new HashIntIntMap[chrs.length*2];
            List<String>[] geneListByChr = new List[chrs.length];
            for (int i = 0; i < chrs.length; i++) {
                geneListByChr[i] = new ArrayList<>();
            }
            for (int i = 0; i < geneNames.length; i++) {
                int cChr = geneNameRangeMap.get(geneNames[i]).chr;
                int index = Arrays.binarySearch(chrs, cChr);
                geneListByChr[index].add(geneNames[i]);
            }  
            for(int i=0;i<chrs.length*2;i++){
                geneRangeMaps[i]=new HashMap();
                posGeneMaps[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
            }
            for (int i = 0; i < chrs.length; i++) {
                String[] genes = geneListByChr[i].toArray(new String[geneListByChr[i].size()]);
                Arrays.sort(genes);
                geneNamesByChr[i] = genes;
                for (int j = 0; j < genes.length; j++) {
                    RangeValStr r = geneNameRangeMap.get(genes[j]);//通过基因的名字来get到基因的range
                    gfIndex=gf.getGeneIndex(genes[j]);
                    geneNameS=genes[j];
                        List<Range> Ranges=new ArrayList<>();
                        for(int a=0;a<gf.genes[gfIndex].ts.size();a++){//出来这个循环之后，这个基因上的每个位点都被遍历一次  
                            if(r.str==(byte)1){//1是正链 0是负链！！！
                                Ranges.add(new Range(i,gf.genes[gfIndex].ts.get(a).getTranscriptEnd()-500,gf.genes[gfIndex].ts.get(a).getTranscriptEnd()+150));
                            }else{
                                Ranges.add(new Range(i,gf.genes[gfIndex].ts.get(a).getTranscriptStart()-150,gf.genes[gfIndex].ts.get(a).getTranscriptStart()+500));
                            }
                        }
                        if(r.str==(byte)1){
                            geneRangeMaps[i*2].put(j,Merge(Ranges));//分正负链将基因的index和这个基因每个转录本的上游150下游500集合起来
                            for(int n=0;n<Merge(Ranges).size();n++){
                                for(int m=Merge(Ranges).get(n).start;m<Merge(Ranges).get(n).end;m++){
                                    posGeneMaps[i*2].put(m, j);//将merge起来的每个位点再与这个基因的index建立hashmap
                                }
                            }
                        }else{
                            geneRangeMaps[i*2+1].put(j,Merge(Ranges));
                            for(int n=0;n<Merge(Ranges).size();n++){
                                for(int m=Merge(Ranges).get(n).start;m<Merge(Ranges).get(n).end;m++){
                                    posGeneMaps[i*2+1].put(m, j);
                                }
                            }
                        }
//                        int a =Arrays.binarySearch(chrs, 1);
//                        System.out.println(a);
                   
                }
            }
        }
        catch (Exception e) {
            System.out.println(gfIndex);
            System.out.println(geneNameS+"\t"+gfIndex);
            e.printStackTrace();
            System.exit(1);    
        }
    }
    public List<Range> Merge(List<Range> Ranges) {
        List<Range> result = new ArrayList<Range>();
        if(Ranges==null||Ranges.size()==0)
            return result;
        Collections.sort(Ranges, new Comparator<Range>(){
            public int compare(Range i1, Range i2){
                if(i1.start!=i2.start)
                    return i1.start-i2.start;
                else
                    return i1.end-i2.end;
            }
        });
        Range pre = Ranges.get(0);
        for(int i=0; i<Ranges.size(); i++){
            Range curr = Ranges.get(i);
            if(curr.start>pre.end){
                result.add(pre);
                pre = curr;
            }else{
                Range merged = new Range(pre.chr,pre.start, Math.max(pre.end, curr.end));
                pre = merged;
            }
        }
        result.add(pre);
        return result;
    }
     private void starAlignmentPair ()  {
        long startTimePoint = System.nanoTime();
//        String subFqDirS = new File ("/data1/home/junxu/total/test-double").getAbsolutePath();
        String subFqDirS = new File (this.outputDirS,subDirS[0]).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
//        for(int i=0;i<fs.length;i++){
//            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
//            fList.add(fs[i]);
//        } 
        int numCores = Runtime.getRuntime().availableProcessors();
        nameSet.stream().forEach(f -> {  
            String infile1 = new File (subFqDirS, f+"_R1.fq").getAbsolutePath();
            String infile2 = new File (subFqDirS, f+"_R2.fq").getAbsolutePath();
            StringBuilder sb = new StringBuilder();
            sb.append(this.starPath).append(" --runThreadN ").append(numCores);
            sb.append(" --genomeDir ").append(this.referenceGenomeDirS);
            sb.append(" --genomeLoad LoadAndKeep");
            sb.append(" --readFilesIn ").append(infile1+" "+infile2);
//            sb.append(" --readFilesCommand zcat ");
            sb.append(" --outFileNamePrefix ").append(new File(new File(this.outputDirS,subDirS[1]).getAbsolutePath(), f)
                .getAbsolutePath()).append(" --outFilterMultimapNmax ").append(1);//this.multiMapN
            sb.append(" --outFilterMismatchNoverLmax ").append(this.mismatchRate)
                .append(" --outFilterIntronMotifs RemoveNoncanonicalUnannotated ");
            sb.append(" --outSAMtype SAM");
            sb.append(" --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin ").append(this.minNMatch);
            String command = sb.toString();
            System.out.println(command);
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    System.out.println(temp);
                }
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
//        StringBuilder time = new StringBuilder();
//        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
//        System.out.println(time.toString());  
    }
    private void starAlignment ()  {
        long startTimePoint = System.nanoTime();
//        String subFqDirS = new File (this.outputDirS, subDirS[0]).getAbsolutePath();
//        String subFqDirS = new File ("/data1/home/junxu/wheat/rnaseq20181204/Out/subFastqs/diffValumn").getAbsolutePath();
        String subFqDirS = new File (this.outputDirS,subDirS[0]).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fList.add(fs[i]);
        }
        int numCores = Runtime.getRuntime().availableProcessors();
//        StringBuilder 
        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append(this.starPath).append(" --runThreadN ").append(numCores);
            sb.append(" --genomeDir ").append(this.referenceGenomeDirS);
            sb.append(" --genomeLoad LoadAndKeep");
//            sb.append(" --sjdbGTFfile ").append(this.gffFile);
//            sb.append(" --sjdbGTFtagExonParentTranscript Parent");
            sb.append(" --readFilesIn ").append(f);
//            sb.append(" --readFilesCommand zcat ");
            sb.append(" --outFileNamePrefix ").append(new File(new File(this.outputDirS,subDirS[1]).getAbsolutePath(), f.getName().replaceFirst(".fq.gz", ""))
                .getAbsolutePath()).append(" --outFilterMultimapNmax ").append(this.multiMapN);
            sb.append(" --outFilterMismatchNoverLmax ").append(this.mismatchRate)
                .append(" --outFilterIntronMotifs RemoveNoncanonicalUnannotated ");
            sb.append(" --outSAMtype SAM");
            sb.append(" --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin ").append(this.minNMatch);
            String command = sb.toString();
            System.out.println(command);
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    System.out.println(temp);
                }
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
        
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());  
    }
    
    private void mkIndexOfReference () {//最终小麦产生index的代码
        int numCores = Runtime.getRuntime().availableProcessors();
        String referenceGenomeFileS = "/data1/home/junxu/wheat/abd_iwgscV1.fa";
//        String outputDirS = "/data1/home/junxu/wheat/starLib1";
        String outputDirS = "/data1/home/junxu/wheat/starLib1.1";
        try {
//            StringBuilder sb = new StringBuilder("/data1/programs/STAR-2.6.0c/bin/Linux_x86_64/STAR");
            StringBuilder sb = new StringBuilder("/data1/home/junxu/wheat/STAR-2.6.1c/bin/Linux_x86_64/STAR");
            sb.append(" --runThreadN ").append(numCores).append(" --runMode genomeGenerate --genomeDir ").append(outputDirS);
//            sb.append(" --sjdbGTFfile ").append("/data1/home/junxu/wheat/changewheat.gtf");//全局变量
            sb.append(" --sjdbGTFfile ").append("/data1/home/junxu/wheat_v1.1_Lulab.gtf");
//            sb.append(" --sjdbGTFtagExonParentTranscript Parent");
            sb.append(" --genomeFastaFiles ").append(referenceGenomeFileS);
            sb.append(" --sjdbOverhang ").append(140);
            sb.append(" --genomeChrBinNbits 17");  
            sb.append(" --genomeSAsparseD 2");
            sb.append(" --limitGenomeGenerateRAM 40000000000");//已经设成13位数了
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
        this.gffFile = pLineList.get(2);
        this.starPath = pLineList.get(3);
        this.outputDirS = pLineList.get(4);
        this.rPath=pLineList.get(5);
//        this.processTaxaAndBarcode();
        for (int i = 0; i < this.subDirS.length; i++) {
            new File(this.outputDirS, subDirS[i]).mkdir();
        }
    }
    
//    private void processTaxaAndBarcode () {
//        RowTable<String> t = new RowTable<>(this.sampleInformationFileS);
//        Set<String> fqSet = new HashSet<>();
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            fqSet.add(t.getCell(i, 3));
//        }
//        fqFileSList = new ArrayList<>(fqSet);
//        Collections.sort(fqFileSList);
//        barcodeLengths = new int[fqFileSList.size()];
//        barcodeLists = new ArrayList[fqFileSList.size()];
//        taxaLists = new ArrayList[fqFileSList.size()];
//        barcodeTaxaMaps = new HashMap[fqFileSList.size()];
//        int[] cnts = new int[fqFileSList.size()];
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
//            cnts[index]++;
//        }
//        for (int i = 0; i < cnts.length; i++) {
//            barcodeLists[i] = new ArrayList<>();
//            taxaLists[i] = new ArrayList<>();
//            barcodeTaxaMaps[i] = new HashMap<>();
//        }
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            int index = Collections.binarySearch(fqFileSList, t.getCell(i, 3));
//            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);
//            allTaxaList.add(taxon);
//            taxaLists[index].add(taxon);
//            barcodeLists[index].add(t.getCell(i, 1));
//            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
//            barcodeLengths[index] = t.getCell(i, 1).length();
//        }
//        Collections.sort(allTaxaList);
//    }
    
    public static void main(String args[]) { 
//        new WheatRNASeq20181107();
//        new WheatTEP(args[0]);
//        new DistinguishSample();
//        new DemoSample();
        new eQTL();
        
    }
    
}
