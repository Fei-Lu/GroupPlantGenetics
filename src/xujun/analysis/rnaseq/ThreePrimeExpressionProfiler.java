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
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import xuebo.analysis.annotation.FStringUtils;

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

    HashMap<String, RangeValStr> geneNameRangeMap = new HashMap<>();
    String[] geneNames = null;
    Integer[] chrs = null;
    String[][] geneNamesByChr = null;
    HashIntIntMap[] posGeneMaps = null;
    HashIntIntMap[] pos3UTRMaps = null;
    HashMap[] geneRangeMaps=null;
    
    int overhangLength = 150;
    
    int multiMapN = 10;
    
    float mismatchRate = (float)0.1;
    
    int minNMatch = 80;
    
    int[] barcodeLengths = null;
    
    String[] subDirS = {"subFastqs", "sams", "geneCount"};
    List<String> fqFileSList = null;
    
    List<String>[] barcodeLists = null;
    
    HashMap<String, String>[] barcodeTaxaMaps = null;
//    HashMap<Integer,Integer> 
    
    List<String>[] taxaLists = null;
    
    List<String> allTaxaList = new ArrayList<String>();
 
    public ThreePrimeExpressionProfiler (String parameterFileS) {
        this.parseParameters(parameterFileS);
//        this.mkPosGeneMap();
        this.parseFq(); //Remove Ts?
        //////this.mkIndexOfReference(); //one-time set, not requried for standard runs
//        this.starAlignment(); 
        this.mkGeneCountTable();
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
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                List<Range> geneRange=new ArrayList();
                int effLenPerReads =0;
                String temp = null;
                int chrIndex = -1;
                int geneIndex = -1;
                String taxon = f.getName().replaceFirst("Aligned.out.sam", "");
                int taxonIndex = Collections.binarySearch(allTaxaList, taxon);
                String[] tem=null;
                int startPos = -1;
                int endPos = -1;
                String cigar=null;int flag=0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("@")) continue;
                    List<String> tList= FStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    chrIndex=Arrays.binarySearch(chrs, Integer.parseInt(tem[2]));//应该是返回这个染色体所在的角标
                    if (chrIndex < 0) continue;
                    startPos=Integer.parseInt(tem[3]);                   
                    cigar=tem[5];
                    flag=Integer.parseInt(tem[1]);
                    int cIndex = 0;
                    int cLength = 0;
                    for (int i = 0; i < cigar.length(); i++) {
                        if (!Character.isDigit(cigar.charAt(i))) {
                            char c = cigar.charAt(i);
                            if (c == 'M' || c == 'D' || c == 'N') {
                                cLength += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            if(c== 'M'){
                                effLenPerReads += Integer.parseInt(cigar.substring(cIndex, i));
                            }
                            cIndex = i+1;
                        }
                    }
                    endPos = startPos+cLength-1;
                    if(flag==16){
                        int index1=posGeneMaps[chrIndex*2].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                effLen[taxonIndex][chrIndex][geneIndex]=effLenPerReads;
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            effLen[taxonIndex][chrIndex][geneIndex]=effLenPerReads;
                                            break;
                                        }
                                    }
                                }
                            }
                        }
                    }else{
                        int index1=posGeneMaps[chrIndex*2+1].get(startPos);
                        if (index1 < 0) continue;
                        int index2=posGeneMaps[chrIndex*2+1].get(endPos);
                        if(index1 == index2){
                            geneIndex=index1;
                            geneRange=(List<Range>) geneRangeMaps[chrIndex*2+1].get(geneIndex);
                            if(geneRange.size()==1){
                                geneCount[chrIndex][geneIndex][taxonIndex]++;
                                effLen[taxonIndex][chrIndex][geneIndex]=effLenPerReads;
                            }else{
                                for(int n=0;n<geneRange.size();n++){
                                    if(geneRange.get(n).isContain(chrIndex, startPos)){
                                        if(geneRange.get(n).isContain(chrIndex, endPos)){
                                            geneCount[chrIndex][geneIndex][taxonIndex]++;
                                            effLen[taxonIndex][chrIndex][geneIndex]=effLenPerReads;
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
        double[][][]rate=new double[allTaxaList.size()][chrs.length][];
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
                        effLength=v.end-v.start+1-effLen[i][j][k]/geneCount[j][k][i];
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
        String outputFileS1 = new File (outfileDirS, "TEP_count-new.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\tChr\tStart\tEnd\tStrand");
            for (int i = 0; i < this.allTaxaList.size(); i++) {
                sb.append("\t").append(allTaxaList.get(i));
            }
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(outputFileS1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < chrs.length; i++) {
                for (int j = 0; j < this.geneNamesByChr[i].length; j++) {
                    sb = new StringBuilder();
                    RangeValStr r = this.geneNameRangeMap.get(geneNamesByChr[i][j]);
                    sb.append(geneNamesByChr[i][j]).append("\t").append(r.chr).append("\t");
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
        }
    }
    
    private void mkPosGeneMap () {
         String gffFile="/Users/xujun/Desktop/TEP/Zea_mays.AGPv4.38.modified.gff3";
//        String gffFile="";
        GeneFeature gf=new GeneFeature(gffFile); 
        String geneNameS=null;int gfIndex=0;
        try{
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
            for(int i=0;i<geneRangeMaps.length;i++){
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
                    if(gf.getGeneBiotype(gfIndex).equals("protein_coding")){
                        List<Range> Ranges=new ArrayList<>();
                        for(int a=0;a<gf.genes[gfIndex].ts.size();a++){//出来这个循环之后，这个基因上的每个位点都被遍历一次  
                            if(r.str==(byte)1){//如果这个基因位于负链上
                                Ranges.add(new Range(i,gf.genes[gfIndex].ts.get(a).getTranscriptStart()-150,gf.genes[gfIndex].ts.get(a).getTranscriptStart()+500));
                            }else{
                                Ranges.add(new Range(i,gf.genes[gfIndex].ts.get(a).getTranscriptEnd()-500,gf.genes[gfIndex].ts.get(a).getTranscriptEnd()+150));
                            }
                        }
                        if(r.str==(byte)1){
                            geneRangeMaps[i*2+1].put(j,Merge(Ranges));
                            for(int n=0;n<Merge(Ranges).size();n++){
                                for(int m=Merge(Ranges).get(n).start;m<Merge(Ranges).get(n).end;m++){
                                    posGeneMaps[i*2+1].put(m, j);
                                }
                            }
                        }else{
                            geneRangeMaps[i*2].put(j,Merge(Ranges));
                            for(int n=0;n<Merge(Ranges).size();n++){
                                for(int m=Merge(Ranges).get(n).start;m<Merge(Ranges).get(n).end;m++){
                                    posGeneMaps[i*2].put(m, j);
                                }
                            }
                        }
                    } 
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            System.out.println(geneNameS+"\t"+gfIndex);
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
        barcodeLengths = new int[fqFileSList.size()];//不同的样本可能带有不同长度的barcode
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
            String taxon = t.getCell(i, 0) + "_"+ t.getCell(i, 2);//这个连接起来的taxon就是我们要找的index信息
            allTaxaList.add(taxon);
            taxaLists[index].add(taxon);
            barcodeLists[index].add(t.getCell(i, 1));
            barcodeTaxaMaps[index].put(t.getCell(i, 1), taxon);
            barcodeLengths[index] = t.getCell(i, 1).length();
        }
        Collections.sort(allTaxaList);
    }
    
    public static void main(String args[]) {
        new ThreePrimeExpressionProfiler(args[0]);
    }
    
}
