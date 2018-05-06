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
import format.range.RangeValStr;
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
import rcaller.RCaller;
import rcaller.RCode;
import utils.Benchmark;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author feilu
 */
public class ThreePrimeExpressionProfiler_new {
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
    String rPath=null;

    HashMap<String, RangeValStr> geneNameRangeMap = new HashMap<>();
    String[] geneNames = null;
    Integer[] chrs = null;
    String[][] geneNamesByChr = null;
    HashIntIntMap[] posGeneMaps = null;
    
    int overhangLength = 150;
    
    int multiMapN = 10;
    
    float mismatchRate = (float)0.1;
    
    int minNMatch = 80;
    
    String[] subDirS = {"subFastqs", "sams", "geneCount","readsLength"};
    List<String> fqFileSList = null;
    
    int[] barcodeLengths = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeTaxaMaps = null;
//    HashMap<Integer,Integer> 
    
    List<String>[] taxaLists = null;
    
    List<String> allTaxaList = new ArrayList<String>();
 
    public ThreePrimeExpressionProfiler_new (String parameterFileS) throws IOException {
        this.parseParameters(parameterFileS);
        mkPosGeneMap();
//        this.parseFq(); //Remove Ts?remove
        //////this.mkIndexOfReference(); //one-time set, not requried for standard runs
//        this.starAlignment(); 
        this.mkGeneCountTable();
    }
    private void mkGeneCountTable1 () {
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon 
        List<Integer>[] readsLength=new ArrayList[allTaxaList.size()];
        for(int i=0;i<allTaxaList.size();i++){
            readsLength[i]=new ArrayList();
        }
        for (int i = 0; i < geneCount.length; i++) { 
//            readsLength=new ArrayList [i][][];
            geneCount[i] = new int[geneNamesByChr[i].length][]; 
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[this.allTaxaList.size()];
            }
        }
        String inputDirS = new File(this.outputDirS, this.subDirS[1]).getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(new File(this.outputDirS,this.subDirS[3]),f.getName().replaceAll("Aligned.out.sam","")+"-1.2000.txt").getAbsolutePath());
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
                int taxonIndex = Collections.binarySearch(allTaxaList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;
                int[] count=new int [2000];
                String cigar=null;
                List w=new ArrayList();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));
                    if (chrIndex < 0||chrIndex>1) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];
                    int cIndex = 0;
                    int cLength = 0;
                    int cLength2=0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            if(c=='M'){
                                cLength2+=Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    int index1=posGeneMaps[chrIndex].get(startPos);
                    if (index1 < 0) continue;
                    int index2=posGeneMaps[chrIndex].get(endPos);
                    if(index1 == index2){                       
                        geneIndex=index1;
                        String name=geneNamesByChr[chrIndex][geneIndex];
                        if(!(w.contains(name))) w.add(name);
                        int indext=w.indexOf(name);
                        count[indext]++;
                        geneCount[chrIndex][geneIndex][taxonIndex]++;
                        readsLength[taxonIndex].add(cLength2);
                    }                  
                }              
                for(int i=0;i<allTaxaList.size();i++){
                    for(int j=0;j<readsLength[i].size();j++){
 //                       System.out.println(readsLength[i].get(j));
                       bw.write(readsLength[i].get(j)+"\n");
                    }                   
                }
            }
            catch (Exception e) {
                e.printStackTrace();
                System.out.println();
            }
        });
    }
    private void mkGeneCountTable () {
        String gffFile="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
        GeneFeature gf=new GeneFeature(gffFile);     
        int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
        double [][][] TPM=new double[chrs.length][][];
        double [][][] effLen=new double [allTaxaList.size()][chrs.length][];
        ArrayList<Integer> [][]readsLengthList = new ArrayList[allTaxaList.size()][];
        for(int i=0;i<allTaxaList.size();i++){
            readsLengthList[i]=new ArrayList[geneNames.length];
            for(int k=0;k<readsLengthList[i].length;k++){
                readsLengthList[i][k]=new ArrayList ();
            }
            for(int j=0;j<chrs.length;j++){
                effLen[i][j]=new double[geneNamesByChr[j].length];
            }
        }
        for (int i = 0; i < geneCount.length; i++) {
            geneCount[i] = new int[geneNamesByChr[i].length][];   
            TPM[i]= new double[geneNamesByChr[i].length][];
            for (int j = 0; j < geneCount[i].length; j++) {
                geneCount[i][j] = new int[this.allTaxaList.size()]; 
                TPM[i][j] = new double[this.allTaxaList.size()];              
            }
        }
        String inputDirS = new File(this.outputDirS, this.subDirS[1]).getAbsolutePath();
 //       String outputDirS = new File(this.outputDirS, this.subDirS[3]).getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);  
        int [][] length=new int [chrs.length][];
        int[][] count=new int[chrs.length][];
        int[][] position=new int[chrs.length][];
        String[][] name=new String[chrs.length][];
        String[][] type=new String[chrs.length][];
//        int[][][] q=new int[allTaxaList.size()][chrs.length][];
              
            for(int j=0;j<chrs.length;j++){
                int geneLength =0;
                length[j]=new int [geneNamesByChr[j].length];
                for(int k=0;k<geneNamesByChr[j].length;k++){                    
                    RangeValStr a=this.geneNameRangeMap.get(geneNamesByChr[j][k]);
                    length[j][k]=geneLength;
                    geneLength+=a.end-a.start+1;               
                }
                count[j]=new int [geneLength];    
                name[j]=new  String [geneLength];
                position[j]=new  int [geneLength];
                type[j]=new  String [geneLength];
            }
        fList.stream().forEach(f -> {
            int gfIndex=0;int CDSIndex=0;
            int UTR5Index=0;int UTR3Index=0;
            String genename=null;
            try {               
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
//                BufferedWriter bw = IOUtils.getTextWriter(new File(new File(this.outputDirS,this.subDirS[3]),f.getName().replaceAll("Aligned.out.sam","")+"-readsLength.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(new File(this.outputDirS,this.subDirS[3]),f.getName().replaceAll("Aligned.out.sam","")+"-onegene.txt").getAbsolutePath());
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
                int taxonIndex = Collections.binarySearch(allTaxaList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;
                List<String>[] w =new ArrayList[ allTaxaList.size()];
                for(int i=0;i<allTaxaList.size();i++){
                    w[i]=new ArrayList();
                }
                
                String cigar=null;    int cnt=0;                          
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));
                    if (chrIndex < 0) continue;
                    cnt++;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];
                    int index1=posGeneMaps[chrIndex].get(startPos);
                    if (index1 < 0) continue;                    
                    int count1=0;
                    int count2=0;
                    int cIndex = 0;
                    int cIndex1=0;
                    int cLength = 0;
                    int cLength2=0;
                    int startPos1=0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            if (c == 'M' || c == 'I' ) {
                                cLength2 += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;                           
                        }
                    }
                    endPos = startPos+cLength-1;                   
                    int index2=posGeneMaps[chrIndex].get(endPos);
                    if(index1 == index2){
                        geneIndex=index1;
                        genename=geneNamesByChr[chrIndex][geneIndex];
                        if(!(w[taxonIndex].contains(genename))) w[taxonIndex].add(genename);
                        int indext=w[taxonIndex].indexOf(genename);
                        geneCount[chrIndex][geneIndex][taxonIndex]++;
                        readsLengthList[taxonIndex][indext].add(cLength2);
                        effLen[taxonIndex][chrIndex][geneIndex]+=cLength2;
                        if(genename.equals("Zm00001d030935")){
                            
                        
                        RangeValStr a = this.geneNameRangeMap.get(geneNamesByChr[chrIndex][geneIndex]);
                        for (int i = 0; i < cigar.length(); i++) {
                            if (!Character.isDigit(cigar.charAt(i))) {
                                char c = cigar.charAt(i);
                                if(c == 'M' ) {
                                    startPos1=startPos+count1;
                                    for(int j=startPos1;j<startPos1+Integer.parseInt(cigar.substring(cIndex1, i));j++){
                                        count[chrIndex][length[chrIndex][geneIndex]+j-a.start]+=1;
   //                                     name[chrIndex][length[chrIndex][geneIndex]+j-a.start]=genename;
   //                                    position[chrIndex][length[chrIndex][geneIndex]+j-a.start]=j;
    //                                    gfIndex=gf.getGeneIndex(genename);
    //                                    for(int k=0;k<gf.genes[gfIndex].ts.size();k++){
    //                                        for(int o=0;o<gf.genes[gfIndex].ts.get(k).cdsList.size();o++){
    //                                            boolean b=gf.genes[gfIndex].ts.get(k).isWithinThisCDS(o, chrIndex, j);
    //                                            if(b){
    //                                                type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="CDS";
    //                                            }
    //                                        }
    //                                        for(int p=0;p<gf.genes[gfIndex].ts.get(k).utr3List.size();p++){
    //                                            boolean e=gf.genes[gfIndex].ts.get(k).isWithinThis3UTR(p, chrIndex, j);
    //                                            if(e){
    //                                                type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="3UTR";
    //                                            }
    //                                        }
    //                                        for(int v=0;v<gf.genes[gfIndex].ts.get(k).utr5List.size();v++){
    //                                            boolean g=gf.genes[gfIndex].ts.get(k).isWithinThis5UTR(v, chrIndex, j);
    //                                                if(g){
    //                                                    type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="5UTR";
    //                                                }
    //                                        }
                                            
 //                                           CDSIndex=gf.genes[gfIndex].ts.get(k).getCDSIndex(chrIndex, j);                                         
 //                                           if(CDSIndex<0){
 //                                               UTR3Index=gf.genes[gfIndex].ts.get(k).get3UTRIndex(chrIndex, j);
 //                                               if(UTR3Index<0){
 //                                                   UTR5Index=gf.genes[gfIndex].ts.get(k).get5UTRIndex(chrIndex, j);
 //                                                   if(UTR5Index<0){
 //                                                       type[chrIndex][length[chrIndex][geneIndex]+j-a.start]=null;
 //                                                   }else{
 //                                                       type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="5UTR";
 //                                                   }                                                   
 //                                               }else{
 //                                                   type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="3UTR";
 //                                              }
 //                                           }else{
 //                                               type[chrIndex][length[chrIndex][geneIndex]+j-a.start]="CDS";
 //                                           }                                           
 //                                       }                                                                                
                                    } 
                                    for(int j=count2;j<count2+Integer.parseInt(cigar.substring(cIndex1, i));j++){
//                                        q[taxonIndex][chrIndex][length[chrIndex][geneIndex]+j-a.start]+=tem[10].charAt(i)-33;
//                                        cnt++;
                                    }
                                    count1+=Integer.parseInt(cigar.substring(cIndex1,i));
                                    count2+=Integer.parseInt(cigar.substring(cIndex1,i));
                                }else{
                                    if(c=='D'||c=='N'){
                                        count1+=Integer.parseInt(cigar.substring(cIndex1,i));
                                    }else{ 
                                        count2+=Integer.parseInt(cigar.substring(cIndex1,i));
                                    }
                                }       
                                cIndex1 = i+1;                           
                            }
                        }
                        }
 //                           bws[taxonIndex].write(p/cont+"\n");
 //                           System.out.print(p/cont);System.out.println();
                        
                    }                  
                }
//                int number=0;
//                for(int i=0;i<chrs.length;i++){
//                    for(int j=0;j<count[i].length;j++){
//                            bw1.write(number+"\t"+name[i][j]+"\t"+position[i][j]+"\t"+type[i][j]+"\t"+count[i][j]);
//                            bw1.newLine();                   
//                        number++;
//                    }
//                }
                if(taxonIndex==10){
                    for(int i=0;i<count[taxonIndex].length;i++){
                       // System.out.print(i+1+"\t"+count[taxonIndex][i]+"\n");
                        bw1.write(i+1+"\t"+count[taxonIndex][i]);bw1.write("\n");
                    }                   
                }
                br.close();
                bw1.flush();bw1.close();

            }
            catch (Exception e) {
                System.out.println(genename+"\t"+gfIndex);
                e.printStackTrace();
            }
        });
        //standerdize
/*        double[][][]rate=new double[allTaxaList.size()][chrs.length][];
        double [] denom=new double[allTaxaList.size()];
        double [] allRate=new double [allTaxaList.size()];
        double effLength=0;
        for(int i=0;i<allTaxaList.size();i++){
            for(int j=0;j<chrs.length;j++){
                rate[i][j] = new double [geneNamesByChr[j].length];
                for(int k=0;k<geneNamesByChr[j].length;k++){        
                    if(geneCount[j][k][i]==0){}
                    else{
                        RangeValStr v = this.geneNameRangeMap.get(geneNamesByChr[j][k]);
                        effLength=v.end-v.start+1-effLen[i][j][k]/geneCount[j][k][i]+1;
                        rate[i][j][k]=Math.log(geneCount[j][k][i])-Math.log(effLength);
                        allRate[i]+=Math.exp(rate[i][j][k]);
                    }                    
                }
            }
            denom[i]=Math.log(allRate[i]);
        }
        for(int i=0;i<chrs.length;i++){
            for(int j=0;j<geneNamesByChr[i].length;j++){
                for(int k=0;k<allTaxaList.size();k++){
                    if(geneCount[i][j][k]==0){TPM[i][j][k]=0;}
                    else{
                        TPM[i][j][k]=Math.exp(rate[k][i][j]-denom[k]+Math.log(1e6));
                    }                   
                }
            }
        }        
        String outfileDirS = new File (this.outputDirS, this.subDirS[2]).getAbsolutePath();
        String outputFileS = new File (outfileDirS, "TEP_count.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tgeneIndexByChr\tStart\tEnd\tStrand");
            for (int i = 0; i < this.allTaxaList.size(); i++) {
                sb.append("\t").append(allTaxaList.get(i)).append("\t").append(allTaxaList.get(i)+":TPM");
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = this.geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append(sb).append("\t");
                    sb.append(r.start).append("\t").append(r.end).append("\t").append(r.str);
                    for (int k = 0; k < this.allTaxaList.size(); k++) {
                        sb.append("\t").append(geneCount[i][j][k]).append("\t").append(TPM[i][j][k]);
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.out.println();
        }*/
    }
    private void findOneGene(){
        String inputFile="/Users/xujun/Desktop/TEP/TEPOut/readsLength/PL_BC4_Plate4-taxon.txt";
        String outputFile="/Users/xujun/Desktop/TEP/TEPOut/readsLength/oneGenePatern.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw = IOUtils.getTextWriter(outputFile);
            while(br.readLine()!=null){
                
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    private void mkMap(){
        try{
            BufferedReader br = IOUtils.getTextReader(geneAnnotationFileS);
            String temp = null;
            Set<String> geneSet = new HashSet<String>();
            Set<Integer> chrSet = new HashSet<>();
            String[] tem = null;
            String  geneName = null;
            byte strand=0;
            HashMap<String, Integer> geneChrMap = new HashMap();
            HashMap<String, Integer> geneMinMap = new HashMap();
            HashMap<String, Integer> geneMaxMap = new HashMap();
            HashMap<String, Byte> geneStrandMap = new HashMap();
            while ((temp = br.readLine())!=null) {
                char s = temp.charAt(0);
                if ((int)s < 48 || (int)s > 57) continue;
                List<String> tList= FStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if(tem[2].startsWith("CDS")){
                    
                }                 
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
            posGeneMaps = new HashIntIntMap[chrs.length];
            List<String>[] geneListByChr = new List[chrs.length];
            for (int i = 0; i < chrs.length; i++) {
                geneListByChr[i] = new ArrayList<>();
            }
            for (int i = 0; i < geneNames.length; i++) {
                int cChr = geneNameRangeMap.get(geneNames[i]).chr;
                int index = Arrays.binarySearch(chrs, cChr);
                geneListByChr[index].add(geneNames[i]);
            }
            for (int i = 0; i < chrs.length; i++) {
                String[] genes = geneListByChr[i].toArray(new String[geneListByChr[i].size()]);
                Arrays.sort(genes);
                geneNamesByChr[i] = genes;
                posGeneMaps[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
                for (int j = 0; j < genes.length; j++) {
                    RangeValStr r = geneNameRangeMap.get(genes[j]);
                    for (int k = r.getRangeStart(); k < r.getRangeEnd(); k++) {
                        posGeneMaps[i].put(k, j);
                    }
                }
            }
            for(int i=0;i<chrs.length;i++){
                
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    private void mkPosGeneMap () {
        try {
            BufferedReader br = IOUtils.getTextReader(geneAnnotationFileS);
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
                geneName=te[1].split(":")[1].substring(0,te[1].split(":")[1].length()-1);
                if (!geneSet.contains(geneName)) {
                    geneMinMap.put(geneName, Integer.MAX_VALUE);
                    geneMaxMap.put(geneName, Integer.MIN_VALUE);
                }
                geneSet.add(geneName);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                if (geneMinMap.get(geneName) > min) geneMinMap.put(geneName, min);
                if (geneMaxMap.get(geneName) < max) geneMaxMap.put(geneName, max);
                int chr = Integer.parseInt(tem[0]);
                chrSet.add(chr);
                geneChrMap.put(geneName, chr);
                if (tem[6].startsWith("-")) geneStrandMap.put(geneName, (byte)1);
                else geneStrandMap.put(geneName, (byte)0);                              
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
            posGeneMaps = new HashIntIntMap[chrs.length];
            List<String>[] geneListByChr = new List[chrs.length];
            for (int i = 0; i < chrs.length; i++) {
                geneListByChr[i] = new ArrayList<>();
            }
            for (int i = 0; i < geneNames.length; i++) {
                int cChr = geneNameRangeMap.get(geneNames[i]).chr;
                int index = Arrays.binarySearch(chrs, cChr);
                geneListByChr[index].add(geneNames[i]);
            }
            for (int i = 0; i < chrs.length; i++) {
                String[] genes = geneListByChr[i].toArray(new String[geneListByChr[i].size()]);
                Arrays.sort(genes);
                geneNamesByChr[i] = genes;
                posGeneMaps[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
                for (int j = 0; j < genes.length; j++) {
                    RangeValStr r = geneNameRangeMap.get(genes[j]);
                    for (int k = r.getRangeStart(); k < r.getRangeEnd(); k++) {
                        posGeneMaps[i].put(k, j);
                    }
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void starAlignment () throws IOException {
        long startTimePoint = System.nanoTime();
        String subFqDirS = new File (this.outputDirS, subDirS[0]).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq");
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fList.add(fs[i]);
        }
        int numCores = Runtime.getRuntime().availableProcessors();
        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append(this.starPath).append(" --runThreadN ").append(numCores).append(" --genomeDir ").append(this.referenceGenomeDirS);
            sb.append(" --sjdbGTFfile ").append(this.geneAnnotationFileS);
            sb.append(" --readFilesIn ").append(f);
            sb.append(" --outFileNamePrefix ").append(new File(new File(this.outputDirS, subDirS[1]).getAbsolutePath(), f.getName().replaceFirst(".fq", ""))
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
    
    private void mkIndexOfReference () {
        int numCores = Runtime.getRuntime().availableProcessors();
        String referenceGenomeFileS = "/Users/feilu/Documents/database/maize/reference/AGPv4/maizeAGPv4.fa";
        String outputDirS = "/Users/feilu/Documents/database/maize/reference/starLib";
        try {
            StringBuilder sb = new StringBuilder("/Users/feilu/Software/STAR-2.5.4b/bin/MacOSX_x86_64/STAR");
            sb.append(" --runThreadN ").append(numCores).append(" --runMode genomeGenerate --genomeDir ").append(outputDirS);
            sb.append(" --sjdbGTFfile ").append(geneAnnotationFileS);
            sb.append(" --genomeFastaFiles ").append(referenceGenomeFileS);
            sb.append(" --sjdbOverhang ").append(140);
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
        this.geneAnnotationFileS = pLineList.get(2);
        this.starPath = pLineList.get(3);
        this.outputDirS = pLineList.get(4);
        this.rPath = pLineList.get(5);
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
            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);
            allTaxaList.add(taxon);
            taxaLists[index].add(taxon);
            barcodeLists[index].add(t.getCell(i, 1));
            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
            barcodeLengths[index] = t.getCell(i, 1).length();
        }
        Collections.sort(allTaxaList);
    }
    
    public static void main(String args[]) throws IOException {
        new ThreePrimeExpressionProfiler_new(args[0]);
    }
    
}
