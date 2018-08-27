/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import com.google.common.collect.Table;
import format.genomeAnnotation.GeneFeature;
import format.range.Range;
import format.table.RowTable;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import graphcis.r.DensityPlot;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class VariantSummary {
    
    public VariantSummary(){
        //this.countSite();
        //this.subSetTest();
        //this.densityTest_deprecated();
        //this.densityTest();
        //this.density();
        //this.filterHmp321Info();
        //this.filterHmp321Info_bysiftTrans();
        //this.summarizeTranscript();
        
        this.summarizeTranscript2();
    }
    
   
    
    private void summarizeTranscript2(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
//        String infileDirS = "/data1/home/aoyue/maizeLoad/001_variantSummary/001_hmp321Info_filter";
//        String geneFeatureFileS = "/data1/publicData/maize/gene/Zea_mays.AGPv4.38.pgf";
//        String outfileS = "/data1/home/aoyue/maizeLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
        double gerpCut = 0;
        File[] fs = new File(infileDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        byte[][] snpAnc = new byte[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        //下面这一段将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里 
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashMap<String, Integer> geneCDSLengthMap = new HashMap();
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的起始位点*/
        List<String> genesList = new ArrayList<>();
        //String[] genes = new String[gf.getGeneNumber()];
        int cntchr11and12 = 0;
        int cntchr1to10 = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrIndex = gf.getGeneChromosome(i)-1;
            if (chrIndex >9) {
                cntchr11and12++;
                continue;
            }
            cntchr1to10++;
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String geneName = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            //genes[i] = geneName;
            genesList.add(geneName);
            List<Range> cdsList = gf.getCDSList(i, longTransIndex); /*得到基因的最长转录本的CDSList*/
            int cnt = 0;
         
            
        /*对于每一个基因的编码序列，我们进行一个for循环，对于编码序列的起始和终止位置，我们又得到一个循环*/    
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
                    否则，如果还位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k);
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList); /*最终将posGeneMap绘图完成*/
                    }
                    cnt++; /*每一个CDS位点加完geneName，cnt就加1，最终cnt是所有cdslist相加的和*/
                }
            }
            geneCDSLengthMap.put(geneName, cnt);
        }
        
        System.out.println(cntchr11and12 + "genes are not used");
        System.out.println(cntchr1to10 + "genes are used");
        
        String[] genes = genesList.toArray(new String[genesList.size()]);
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length];
        List<File> hmpList = Arrays.asList(fs);
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                /*这一段主要是将有基因的位点列出来，及将转录组的位点列出来*/
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt%1000000 == 0) System.out.println("Hmp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ###hmpInfo Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains("<") || l.get(3).contains(",")) continue;
                    int pos = Integer.valueOf(l.get(1));
                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));
                        snpCount[index]++; //该位点有几个基因，就有几次snpCount.
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(l.get(4).getBytes()[0]);  //返回该位点的ancestral碱基的二进制码
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray();
            snps[chrIndex] = snpList.toArray();
            snpAnc[chrIndex] = snpAncList.toArray();
        });
        
        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length];
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];
        
        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];
        
        TIntArrayList[] delPosList = new TIntArrayList[chrNum];
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            delPosList[chrIndex] = new TIntArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Sift\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ### SIFT Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(16).startsWith("NA")) continue;
                    if (l.get(18).startsWith("NA")) continue;
                    String gene = l.get(18);
                    int geneIndex = Arrays.binarySearch(genes, gene); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) continue;
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                    if (index < 0) continue;
                    if (snps[chrIndex][index] != l.get(3).getBytes()[0]) continue;
                    byte ref = l.get(2).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                        derivedState = 1; //mean b73 carries derived allele
                    }
                    else if (snpAnc[chrIndex][index] == ref) {
                        derivedState = 0;
                    }
                    
                    if (derivedState == -1) noAncCount[geneIndex]++;
                    
                    String type = null;
                    if (l.get(16).equals("NA")) {
                        
                    }
                    else {
                        if (l.get(16).equals("SYNONYMOUS")) {
                            type = "Syn";
                            synCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73SynCount[geneIndex]++;
                            }
                        }
                        else {
                            type = "Non";
                            nonCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73NonCount[geneIndex]++;
                            }
                            if (l.get(17).startsWith("NA")) {
                                    naCount[geneIndex]++;
                            }
                            else{
                                if (Double.valueOf(l.get(17)) < 0.05) {
                                    delCount[geneIndex]++;
                                    delPosList[chrIndex].add(pos);
                                    if (derivedState == 1) {
                                        b73DelCount[geneIndex]++;
                                    }
                                }
                            }
                        }
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        int[][] delPos = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delPos[i] = delPosList[i].toArray();
            Arrays.sort(delPos[i]);
        }
        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpTree = new double[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpTree = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];
        
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) { //gerp文件没有表头
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Gerp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ### Gerp Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    if(l.get(14).equals("NA")) continue;
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        if (gene == null) continue;
                        int geneIndex = Arrays.binarySearch(genes, gene);
                        
                        double treeValue = Double.valueOf(l.get(14));
                        double scoreValue = Double.valueOf(l.get(15));
                        if (treeValue == 0) continue;
                        gerpAlignCount[geneIndex]++;
                        gerpTree[geneIndex]+=treeValue;
                        gerpScore[geneIndex]+=scoreValue;
                        int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                        if (index < 0) continue;

                        byte derivedState = 0; //mean ancestral allele is not defined or ancestral allele is alt
                        if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                            derivedState = 1; //mean b73 carries derived allele
                        }
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpTree[geneIndex]+=treeValue;
                        snpGerpScore[geneIndex]+=scoreValue;
                        index = Arrays.binarySearch(delPos[chrIndex], pos);
                        if (index < 0) continue;
                        if (scoreValue <= gerpCut) continue;
                        delHGCount[geneIndex]++;
                        if (derivedState == 1) {
                            b73DelHGCount[geneIndex]++;
                        }
                    }
                    
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
     
        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            header = header +"\tNumAmbigousAnc\tB73NumberOfSyn\tB73PercentageSyn\tB73NumberOfNon\tB73PercentageNon\tB73NumberOfDeleterious\tB73PercentageDeleterious\tB73NumberOfHGDeleterious\tB73PercentageHGDeleterious";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < genes.length; i++) {
                //if(gf.getGeneChromosome(i) > 10) continue;
                StringBuilder sb =new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double)snpCount[i]/cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) ifSiftAligned = 0;
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double)synCount[i]/cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double)nonCount[i]/cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) ratio = Double.NaN;
                else ratio = (double)nonCount[i]/synCount[i];
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double)delCount[i]/cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double)delHGCount[i]/cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) ifGerpAligned = 0;
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double)gerpAlignCount[i]/cdsLength).append("\t");
                if (gerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN).append("\t");
                else sb.append((double)gerpTree[i]/gerpAlignCount[i]).append("\t").append((double)gerpScore[i]/gerpAlignCount[i]).append("\t");
                if (snpGerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN);
                else sb.append((double)snpGerpTree[i]/snpGerpAlignCount[i]).append("\t").append((double)snpGerpScore[i]/snpGerpAlignCount[i]);
                
                double cdsL = (double)(snpCount[i]-noAncCount[i])/snpCount[i]*cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(b73SynCount[i]).append("\t").append((double)b73SynCount[i]/cdsL).append("\t");
                sb.append(b73NonCount[i]).append("\t").append((double)b73NonCount[i]/cdsL).append("\t");
              
                sb.append(b73DelCount[i]).append("\t").append((double)b73DelCount[i]/cdsL).append("\t");
                sb.append(b73DelHGCount[i]).append("\t").append((double)b73DelHGCount[i]/cdsL);
                
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
   
    private void summarizeTranscript(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
//        String infileDirS = "/data1/home/aoyue/maizeLoad/001_variantSummary/001_hmp321Info_filter";
//        String geneFeatureFileS = "/data1/publicData/maize/gene/Zea_mays.AGPv4.38.pgf";
//        String outfileS = "/data1/home/aoyue/maizeLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
        double gerpCut = 0;
        File[] fs = new File(infileDirS).listFiles();
        int chrNum = fs.length;
        HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum];
        int[][] snpPos = new int[chrNum][];
        byte[][] snps = new byte[chrNum][];
        byte[][] snpAnc = new byte[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            posGeneMap[i] = new HashMap();
        }
        //下面这一段将posGeneMap建立完整，使每个位点对应哪些基因名字，都装进这个map里 
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        HashMap<String, Integer> geneCDSLengthMap = new HashMap();
        /*将所有基因的名字进行for循环输入到数组genes中，对应于每一个基因，我们通过getTranscriptName得到转录本的名字，通过getCDSList方法得到编码序列的起始位点*/
        String[] genes = new String[gf.getGeneNumber()];
        int cntchr11and12 = 0;
        int cntchr1to10 = 0;
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String geneName = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            genes[i] = geneName;
            List<Range> cdsList = gf.getCDSList(i, longTransIndex); /*得到基因的最长转录本的CDSList*/
            int cnt = 0;
            int chrIndex = gf.getGeneChromosome(i)-1;
            if (chrIndex >9) {
                cntchr11and12++;
                continue;
            }
            cntchr1to10++;
        /*对于每一个基因的编码序列，我们进行一个for循环，对于编码序列的起始和终止位置，我们又得到一个循环*/    
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
                    否则，如果还位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k);
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList); /*最终将posGeneMap绘图完成*/
                    }
                    cnt++; /*每一个CDS位点加完geneName，cnt就加1，最终cnt是所有cdslist相加的和*/
                }
            }
            geneCDSLengthMap.put(geneName, cnt);
        }
        
        System.out.println(cntchr11and12 + "genes are not used");
        System.out.println(cntchr1to10 + "genes are used");
        
        Arrays.sort(genes); //genes指所有基因对应最长转录本的名字的组合，是一个数组。
        int[] snpCount = new int[genes.length];
        List<File> hmpList = Arrays.asList(fs);
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            TIntArrayList snpPosList = new TIntArrayList();
            TByteArrayList snpList = new TByteArrayList();
            TByteArrayList snpAncList = new TByteArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                /*这一段主要是将有基因的位点列出来，及将转录组的位点列出来*/
                while ((temp = br.readLine()) != null) {
                    cnt++; //检测程序运行的情况
                    if (cnt%1000000 == 0) System.out.println("Hmp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ###hmpInfo Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains("<") || l.get(3).contains(",")) continue;
                    int pos = Integer.valueOf(l.get(1));
                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));
                        snpCount[index]++; //该位点有几个基因，就有几次snpCount.
                    }
                    /*将该位置加入snpPosList，表明该位置是有基因的，没有基因的位置就不必加到snpPosList中去*/
                    snpPosList.add(pos);
                    snpList.add(l.get(3).getBytes()[0]); //返回该位点的alt碱基的二进制码
                    snpAncList.add(l.get(4).getBytes()[0]);  //返回该位点的ancestral碱基的二进制码
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            snpPos[chrIndex] = snpPosList.toArray();
            snps[chrIndex] = snpList.toArray();
            snpAnc[chrIndex] = snpAncList.toArray();
        });
        
        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length];
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];
        
        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];
        
        TIntArrayList[] delPosList = new TIntArrayList[chrNum];
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            delPosList[chrIndex] = new TIntArrayList();
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Sift\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ### SIFT Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(16).startsWith("NA")) continue;
                    if (l.get(18).startsWith("NA")) continue;
                    String gene = l.get(18);
                    int geneIndex = Arrays.binarySearch(genes, gene); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) continue;
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                    if (index < 0) continue;
                    if (snps[chrIndex][index] != l.get(3).getBytes()[0]) continue;
                    byte ref = l.get(2).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                        derivedState = 1; //mean b73 carries derived allele
                    }
                    else if (snpAnc[chrIndex][index] == ref) {
                        derivedState = 0;
                    }
                    
                    if (derivedState == -1) noAncCount[geneIndex]++;
                    
                    String type = null;
                    if (l.get(16).equals("NA")) {
                        
                    }
                    else {
                        if (l.get(16).equals("SYNONYMOUS")) {
                            type = "Syn";
                            synCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73SynCount[geneIndex]++;
                            }
                        }
                        else {
                            type = "Non";
                            nonCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73NonCount[geneIndex]++;
                            }
                            if (l.get(17).startsWith("NA")) {
                                    naCount[geneIndex]++;
                            }
                            else{
                                if (Double.valueOf(l.get(17)) < 0.05) {
                                    delCount[geneIndex]++;
                                    delPosList[chrIndex].add(pos);
                                    if (derivedState == 1) {
                                        b73DelCount[geneIndex]++;
                                    }
                                }
                            }
                        }
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        
        int[][] delPos = new int[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delPos[i] = delPosList[i].toArray();
            Arrays.sort(delPos[i]);
        }
        int[] gerpAlignCount = new int[genes.length];
        int[] snpGerpAlignCount = new int[genes.length];
        double[] gerpTree = new double[genes.length];
        double[] gerpScore = new double[genes.length];
        double[] snpGerpTree = new double[genes.length];
        double[] snpGerpScore = new double[genes.length];
        
        hmpList.stream().forEach(f -> {
            int chrIndex = Integer.valueOf(f.getName().split("_chr")[1].replaceFirst("_AGPv4_AnnoDB.txt", ""))-1;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                String temp = br.readLine();
                int cnt = 0;
                while ((temp = br.readLine()) != null) { //gerp文件没有表头
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println("Gerp\tchr"+String.valueOf(chrIndex+1)+"\t"+String.valueOf(cnt) + " ### Gerp Process");
                    List<String> l = PStringUtils.fastSplit(temp);
                    int pos = Integer.valueOf(l.get(1));
                    if(l.get(14).equals("NA")) continue;
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        if (gene == null) continue;
                        int geneIndex = Arrays.binarySearch(genes, gene);
                        
                        double treeValue = Double.valueOf(l.get(14));
                        double scoreValue = Double.valueOf(l.get(15));
                        if (treeValue == 0) continue;
                        gerpAlignCount[geneIndex]++;
                        gerpTree[geneIndex]+=treeValue;
                        gerpScore[geneIndex]+=scoreValue;
                        int index = Arrays.binarySearch(snpPos[chrIndex], pos);
                        if (index < 0) continue;

                        byte derivedState = 0; //mean ancestral allele is not defined or ancestral allele is alt
                        if (snpAnc[chrIndex][index] == snps[chrIndex][index]) {
                            derivedState = 1; //mean b73 carries derived allele
                        }
                        snpGerpAlignCount[geneIndex]++;
                        snpGerpTree[geneIndex]+=treeValue;
                        snpGerpScore[geneIndex]+=scoreValue;
                        index = Arrays.binarySearch(delPos[chrIndex], pos);
                        if (index < 0) continue;
                        if (scoreValue <= gerpCut) continue;
                        delHGCount[geneIndex]++;
                        if (derivedState == 1) {
                            b73DelHGCount[geneIndex]++;
                        }
                    }
                    
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
     
        gf.sortGeneByStartPosition();
        try {
            String header = "Transcript\tCDSLength\tSNPNumber\tSNPPercentage\tIfSiftAligned\tNumberOfSyn\tPercentageSyn\tNumberOfNon\tPercentageNon\tNonVsSynRatio\tNumberOfDeleterious\tPercentageDeleterious\tNumberOfHGDeleterious\tPercentageHGDeleterious\tIfGerpAligned\tGerpAlignedCount\tPercentageGerpAlignedCount\tMeanGerpTreeLength\tMeanGerpScore\tSNPMeanGerpTreeLength\tSNPMeanGerpScore";
            header = header +"\tNumAmbigousAnc\tB73NumberOfSyn\tB73PercentageSyn\tB73NumberOfNon\tB73PercentageNon\tB73NumberOfDeleterious\tB73PercentageDeleterious\tB73NumberOfHGDeleterious\tB73PercentageHGDeleterious";
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < (genes.length - cntchr11and12); i++) {
                if(gf.getGeneChromosome(i) > 10) continue;
                StringBuilder sb =new StringBuilder(genes[i]);
                int cdsLength = geneCDSLengthMap.get(genes[i]);
                sb.append("\t").append(cdsLength).append("\t").append(snpCount[i]).append("\t").append((double)snpCount[i]/cdsLength).append("\t");
                int ifSiftAligned = 1;
                if (naCount[i] == nonCount[i]) ifSiftAligned = 0;
                sb.append(ifSiftAligned).append("\t").append(synCount[i]).append("\t").append((double)synCount[i]/cdsLength).append("\t");
                sb.append(nonCount[i]).append("\t").append((double)nonCount[i]/cdsLength).append("\t");
                double ratio = 0;
                if (synCount[i] == 0) ratio = Double.NaN;
                else ratio = (double)nonCount[i]/synCount[i];
                sb.append(ratio).append("\t").append(delCount[i]).append("\t").append((double)delCount[i]/cdsLength).append("\t");
                sb.append(delHGCount[i]).append("\t").append((double)delHGCount[i]/cdsLength).append("\t");
                int ifGerpAligned = 1;
                if (gerpAlignCount[i] == 0) ifGerpAligned = 0;
                sb.append(ifGerpAligned).append("\t").append(gerpAlignCount[i]).append("\t").append((double)gerpAlignCount[i]/cdsLength).append("\t");
                if (gerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN).append("\t");
                else sb.append((double)gerpTree[i]/gerpAlignCount[i]).append("\t").append((double)gerpScore[i]/gerpAlignCount[i]).append("\t");
                if (snpGerpAlignCount[i] == 0) sb.append(Double.NaN).append("\t").append(Double.NaN);
                else sb.append((double)snpGerpTree[i]/snpGerpAlignCount[i]).append("\t").append((double)snpGerpScore[i]/snpGerpAlignCount[i]);
                
                double cdsL = (double)(snpCount[i]-noAncCount[i])/snpCount[i]*cdsLength;
                sb.append("\t").append(noAncCount[i]).append("\t").append(b73SynCount[i]).append("\t").append((double)b73SynCount[i]/cdsL).append("\t");
                sb.append(b73NonCount[i]).append("\t").append((double)b73NonCount[i]/cdsL).append("\t");
              
                sb.append(b73DelCount[i]).append("\t").append((double)b73DelCount[i]/cdsL).append("\t");
                sb.append(b73DelHGCount[i]).append("\t").append((double)b73DelHGCount[i]/cdsL);
                
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
    
    
    private void filterHmp321Info_bysiftTrans(){
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
        String outDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/002_hmp321SiftTrans_filter";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "hmp");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, "hmp321Info_filterbySift_" + f.getName().split("_filter_")[1].replaceFirst("_AGPv4_AnnoDB.txt", ".txt")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            try{
                String header = br.readLine();
                bw.write(header);
                bw.newLine();
                String temp = null;
                List<String> l = null;
                int cnt = 0;
                while((temp = br.readLine()) != null){
                    l = PStringUtils.fastSplit(temp);
                    String trans = l.get(17);
                    if(trans.equals("NA"))continue;
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt)+"\t"+ f.getName().split("_filter_")[1].replaceFirst("_AGPv4_AnnoDB.txt", "") + " trans sites");
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
                
            }
        });
        
        
    }
    
    private void filterHmp321Info () {
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_new";
        String outDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/testtt";
        int minMinorDepth = 3;
        double maxIndiDepth = 7;
        int siteDepthCut = 5000;
        double siteHeterozygousCut = 0.15;
        File[] fs = new File(inDirS).listFiles();
        List<File> fList = Arrays.asList(fs);
        
        fList.stream().forEach(f -> {
            String outfileS = new File(outDirS, f.getName().replaceFirst("hmp321Info", "hmp321Info_filter")).getAbsolutePath();
            String temp = null;
            String checkNA = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write(br.readLine()); //把表头读进去
                bw.newLine();
                //String temp = null;
                int cnt = 0;
                int cntfilter_number =0;
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    if((checkNA = l.get(5)).equals("NA")){
                        cnt++; 
                    }
                    else{
                        cntfilter_number++;
                        if (l.get(3).contains(",")) continue;
                        if (Integer.valueOf(l.get(10)) < minMinorDepth) continue;  // 第10列是MaxMinorDepth Remove sites with maximum minor allele count 
                        if (Integer.valueOf(l.get(11)) > siteDepthCut) continue;
                        double siteHeterozygous = Double.valueOf(l.get(13))/Double.valueOf(l.get(12)); //12列是位点的次数，13列是位点的杂合子数
                        if (siteHeterozygous > siteHeterozygousCut) continue;
                        double indiDepth = Double.valueOf(l.get(11))/Double.valueOf(l.get(12)); // 第11列是该位点的测序深度。12列是该位点有几个taxa被测到。
                        if (indiDepth>maxIndiDepth) continue;
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName().split("_")[1] + ": " + cnt + "sites are NA about basic annotation.");
                System.out.println(f.getName().split("_")[1] + ": " + cntfilter_number + "sites are filter under this standard.");
                
            }
            catch (Exception e) {
                e.printStackTrace();
                System.out.println(temp);
            }
        });
    }
    
    
    public void density(){
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_new";
        String outDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/SNPSiteDensity";
        File[] fs = IOUtils.listFilesStartsWith(new File(inDirS).listFiles(), "hmp");
        List<File> fsList = Arrays.asList(fs);
        fsList.stream().forEach(f -> {
            String outfileS = new File(outDirS,f.getName().replaceFirst("_AnnoDB.txt", ".density.pdf")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
        try{
            String header = br.readLine();
            String temp = null;
            List<String> l = null;
            List<Double> SiteDepthList = new ArrayList();
            List<Double> SiteCountList = new ArrayList();
            int cnt = 0;
            while((temp = br.readLine()) != null){
                double r = Math.random();
                if (r > 0.003) continue;
                l = PStringUtils.fastSplit(temp);
                if((temp = l.get(11)).equals("NA")){
                    cnt++;
                }
                else{
                    SiteDepthList.add(Double.parseDouble(l.get(11)));
                    SiteCountList.add(Double.parseDouble(l.get(12)));
                }
            }
            System.out.println(f.getName().split("_")[1] + ": " + cnt + " is NA");
            br.close();
            double[] depth = new double[SiteDepthList.size()];
            double[] siteCount = new double[SiteCountList.size()];
            for(int i =0; i< depth.length;i++){
                depth[i] = SiteDepthList.get(i);
                siteCount[i] = SiteCountList.get(i);
            }
            double[] d = new double[depth.length];
            for (int i = 0; i < d.length; i++) {
                d[i] = depth[i]/siteCount[i]; /*平均一个样在1个位点测了几次*/
            }
            DensityPlot dd = new DensityPlot(d);
            dd.setTitle("Kernel Density of Coverage of " + f.getName().split("_")[1]);
            dd.setXLim(0, 20);
            dd.setYLim(0, 0.4);
            dd.setXLab("Coverage");
            dd.setYLab("Density");
            dd.saveGraph(outfileS);
            dd.showGraph();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        });
        
        
    }
    public void densityTest(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/testfilter/test.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/testfilter/d.pdf";
        BufferedReader br = IOUtils.getTextReader(infileS);
        try{
            String header = br.readLine();
            String temp = null;
            List<String> l = null;
            List<Double> SiteDepthList = new ArrayList();
            List<Double> SiteCountList = new ArrayList();
            int cnt = 0;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                if((temp = l.get(11)).equals("NA")){
                    cnt++;
                }
                else{
                    SiteDepthList.add(Double.parseDouble(l.get(11)));
                    SiteCountList.add(Double.parseDouble(l.get(12)));
                }
            }
            System.out.println(cnt + " is NA");
            br.close();
            double[] depth = new double[SiteDepthList.size()];
            double[] siteCount = new double[SiteCountList.size()];
            for(int i =0; i< depth.length;i++){
                depth[i] = SiteDepthList.get(i);
                siteCount[i] = SiteCountList.get(i);
            }
            double[] d = new double[depth.length];
            for (int i = 0; i < d.length; i++) {
                d[i] = depth[i]/siteCount[i]; /*平均一个样在1个位点测了几次*/
            }
            DensityPlot dd = new DensityPlot(d);
            dd.setTitle("Kernel Density of Coverage Per chromosome");
            dd.setXLab("Coverage");
            dd.setYLab("Density");
            dd.saveGraph(outfileS);
            dd.showGraph();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void densityTest_deprecated() {
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/testfilter/test_ok.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/testfilter/test_removeNAok.txt";
        RowTable<String> t = new RowTable<> (infileS);
//        List<String> l = t.getHeader();
//        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//        try{
//            StringBuilder sb = new StringBuilder();
//            for(int i =0; i< l.size(); i++){
//                sb.append(l.get(i)).append("\t");
//            }
//            sb.deleteCharAt(sb.length()-1);
//            bw.write(sb.toString());
//            bw.newLine();
//            bw.flush();
//            bw.close();
//        }
//        catch(Exception e){
//            e.printStackTrace();
//            System.exit(1);
//        }
//        RowTable<String> t1 = new RowTable<> (outfileS);
        int cnt =0;
        for(int i=0; i< t.getRowNumber();i++){
            if(String.valueOf(t.getCell(i, 11)).equals("NA")){
                cnt++;
                System.out.println(i + " is NA");
                t.removeRow(i);             
            }
        }
        System.out.println(cnt + " is NA");
        t.writeTextTable(outfileS, IOFileFormat.Text);
        RowTable<String> t1 = new RowTable<> (outfileS);
        double[] depth = t1.getColumnAsDoubleArray(11); /*获取table第11列的SiteDepth值，组成一个数组。1210份taxa测序5X，该位点总共测的次数*/
        double[] siteCount = t1.getColumnAsDoubleArray(12);/*获取table第11列的SiteCount值，组成一个数组。该位点有几个taxa测到*/
        double[] d = new double[depth.length];
        for (int i = 0; i < d.length; i++) {
            d[i] = depth[i]/siteCount[i]; /*平均一个样在1个位点测了几次*/
        }
        DensityPlot dd = new DensityPlot(d);
        dd.showGraph();
    }
    
    private void subSetTest () {
/*Chr	Pos	Ref	Alt	Ancestral(Gerp)	Major	Minor	MinorAlleleFrequency	MinorPresence	MinorTotalDepth	MaxMinorDepth	SiteDepth	SiteCount	HetCount
10	228730	A	T	NA	A	T	0.013307985	21	64	16	3591	789	1*/
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_newnew/hmp321Info_chr010_AGPv4_AnnoDB.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/testfilter/test_ok.txt";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                double r = Math.random();
                if (r > 0.002) continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(3).contains(",")) continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void countSite () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_new";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if(fs[i].isHidden()){
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles(); //将文件删除后，重新将文件列表打印出来，此时，fs不包含隐藏文件。
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = -1;
            try {
                BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    cnt++;
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(String.valueOf(cnt)+"\t"+fs[i].getName());
            sum+=cnt;
        }
        System.out.println(String.valueOf(sum));
    }
}
