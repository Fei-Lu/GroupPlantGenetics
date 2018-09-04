/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import com.google.common.collect.Table;
import format.genomeAnnotation.GeneFeature;
import format.range.Range;
import format.range.Ranges;
import format.table.RowTable;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TDoubleArrayList;
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
import static utils.IOFileFormat.Text;
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
        //this.summarizeTranscript2();
       //this.classifySNPs();
       //this.test();
       this.mkBarplotOfSNPs();
    }
    
  
    
    private void mkBarplotOfSNPs () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class";
        String countFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/classCount.txt";
        String mafDistrubutionFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/mafSFS.txt";
        String dafDistrubutionFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/dafSFS.txt";
        //int sampleSize = 10000;
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if(fs[i].isHidden()){
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles(); //将文件删除后，重新将文件列表打印出来，此时，fs不包含隐藏文件。
        /*建立一个边界数组bound，大小是100；bound[i] = 1/100*i;即，均分为100等分！
        建立一个二维数组mafFrequency 和 dafFrequency， 长度为class的种类长，100宽；
        建立一个count 和 dafCount 数组，长度为class分类长；
        建立一个 mafList 和 dafList 数组，长度为class分类长；
        
        进入for循环，对class文件一一遍历，以读表格的形式进入文件
        Chr	Pos	MinorAllele	MAF	DerivedAllele	DAF
        1	92716	A	0.0077619664	NA	NA
        1	93774	G	8.130081E-4	A	0.9991869919
        1	122123	C	0.0012254902	C	0.0012254902
        做以下几件事情：1，数行数，看每个分类的个数；每个位点一定有maf的值，但不一定有daf的值，因为DA allele不一定存在，若不存在，则DAF值为空。
        2，将maf的值加入mafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则mafFrequency[i][index]++;
        如果DAF不以N开头，则将daf的值加入 dafList;在bound数组里搜索，如果未搜到，则index = - index -2; 否则dafFrequency[i][index]++; dafCount数组加一
        在此循环内，进入进入另一个for循环j，做统计：
        计算出maf 和daf 在 index 1-100 bin范围内,各个位点所占的百分比。即 index=1, 所有位点count= rownumber， frenquency = index /count；
        */
        int size = 100;
        double[] bound = new double[size];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double)1/size*i;
        }
        double[][] mafFrequency = new double[fs.length][size];
        double[][] dafFrequency = new double[fs.length][size];
        int[] count = new int[fs.length];
        int[] dafCount = new int[fs.length];
        TDoubleArrayList[] mafList = new TDoubleArrayList[fs.length];
        TDoubleArrayList[] dafList = new TDoubleArrayList[fs.length];
        for (int i = 0; i < fs.length; i++) {
            mafList[i] = new TDoubleArrayList(); //这里为何要再new一下？
            dafList[i] = new TDoubleArrayList(); //这里为何要再new一下？
            String infileS = fs[i].getAbsolutePath(); 
            RowTable<String> t = new RowTable<>(infileS);
            count[i] = t.getRowNumber();
            for (int j = 0; j < t.getRowNumber(); j++) {
                double value = t.getCellAsDouble(j, 3);
                mafList[i].add(value);
                int index = Arrays.binarySearch(bound, value);
                if (index < 0) index = -index -2;
                //例如 0.112021856在bound搜索结果为-13，则此时index为11，及0.1-0.2范围内。好神奇！！又如0.112394，index也是11.
                //又如0.21652在bound搜索结果中为-23,这样index=21， 这样就将maf的值按照1-100分布开来。
                mafFrequency[i][index]++; 
                if (!t.getCell(j, 5).startsWith("N")) {
                    value = t.getCellAsDouble(j, 5);
                    dafList[i].add(value);
                    index = Arrays.binarySearch(bound, value);
                    if (index < 0) index = -index -2;
                    dafFrequency[i][index]++;
                    dafCount[i]++;
                }
            }
            for (int j = 0; j < mafFrequency[i].length; j++) {
                mafFrequency[i][j] = mafFrequency[i][j]/count[i];
                dafFrequency[i][j] = dafFrequency[i][j]/dafCount[i];
            }
        }
        /*打表格，输入表头，去掉最后一个\t键；表头为4个文件的文件名；
        第二行输入每个分类的conut数；
        */
        try {
            BufferedWriter bw =IOUtils.getTextWriter(countFileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb.append(fs[i].getName().replaceFirst(".txt", "")).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length-1; i++) {
                bw.write(String.valueOf(count[i])+"\t");
            }
            bw.write(String.valueOf(count[count.length-1]));
            bw.newLine();
            bw.flush();
            bw.close();
            bw =IOUtils.getTextWriter(mafDistrubutionFileS);
            bw.write("MAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t"+fs[i].getName().replaceFirst(".txt", ""));
            }
            bw.newLine();
            /*double[][] mafFrequency = new double[fs.length][size] 文件长度为4，size为100*/
            for (int i = 0; i < mafFrequency[0].length; i++) { //i小于第一个文件的长度100
                bw.write(String.valueOf(bound[i]));
                for (int j = 0; j < mafFrequency.length; j++) { //j小于400，将1-100的频率写出来
                    bw.write("\t"+mafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw =IOUtils.getTextWriter(dafDistrubutionFileS); 
            bw.write("DAF");
            for (int i = 0; i < fs.length; i++) {
                bw.write("\t"+fs[i].getName().replaceFirst(".txt", ""));
            }
            bw.newLine();
            for (int i = 0; i < dafFrequency[0].length; i++) {
                bw.write(String.valueOf(bound[i]));
                for (int j = 0; j < dafFrequency.length; j++) {
                    bw.write("\t"+dafFrequency[j][i]);
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    
    public void test(){
        int size = 100;
        double[] bound = new double[size];
        for (int i = 1; i < bound.length; i++) {
            bound[i] = (double)1/size*i;
        }
        double b = 0.112021856;
        double c = 0.21652;
        double d = 0.112394;
        int index = Arrays.binarySearch(bound, b);
        int index1 = Arrays.binarySearch(bound, c);
        int index2 = Arrays.binarySearch(bound, d);
        int a =3;
    }
    
      /**
     * Step 1: filtered out low quality gene models. Non/syn ratio less than 2.5
     * Step 2: classify SNPs into different classes based on only high confidence gene models. Deleterious mutations (SIFT less than 0.05, GERP greater than 2)
     */
    
    
    private void classifySNPs () {
        String geneFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String geneSummaryFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/003_transcriptSummary/transcriptSummary.txt";
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/002_hmp321SiftTrans_filter";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class";
        String transcriptFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String[] types = {"Synonymous", "Non_Synonymous_Tolerent", "Non_Synonymous_Deleterious", "Non_Synonymous_Deleterious_High_GERP"};
        
        double gerpCut = 0;
        double nonSynRatioCut = 2.5;//filter incorrect gene model
        double nonSynGeneCut = 0.02;
        RowTable t = new RowTable (geneSummaryFileS);
        boolean[] ifOut = new boolean[t.getRowNumber()];
        ArrayList<String> tranNameList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (!t.getCell(i, 4).equals("1")) continue; //有sift信息
            if (!t.getCell(i, 14).equals("1")) continue; //filter incorrect gene model by picking up gerp aligned genes
            double ratio = t.getCellAsDouble(i, 9); //非同义突变除以同义突变
            double nonSyn = t.getCellAsDouble(i, 8);
            if (Double.isNaN(ratio)) { //说明比值是一个NaN，此基因中没有同义突变，只有非同义突变， 若非同义突变大于0.02就删去不要
                if (nonSyn > nonSynGeneCut) continue;
                tranNameList.add(t.getCellAsString(i, 1));
                ifOut[i] = true;
            }
            else {
                if (ratio > nonSynRatioCut) continue; /*非同义突变除以同义突变 大于2.5，则删除不要*/
                tranNameList.add(t.getCellAsString(i, 0));
                ifOut[i] = true; //如果是对的，就将true赋值给ifOut
            }
        }
        t.writeTextTable(transcriptFileS, Text, ifOut);
        String[] tranName = tranNameList.toArray(new String[tranNameList.size()]);
        System.out.println(String.valueOf(tranName.length) + "genes kept");
        Arrays.sort(tranName);
        /*对gf文件进行过滤，挑出高质量的gene model，如果gf文件中的基因在tranName中搜到，代表是高质量的基因，后续继续得到该基因的cdsList,
            根据该基因的cdsList，我们知道了该基因的区域，后续判断中，如果输入文件是在该区域，则进行sift等信息的统计*/
        GeneFeature gf = new GeneFeature(geneFileS);
        gf.sortGeneByStartPosition();
        ArrayList<Range> cdsList = new ArrayList();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String query = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            //String query = gf.getTranscriptName(i, 0);
            int index = Arrays.binarySearch(tranName, query);
            if (index < 0) continue;
            List<Range> cdsl = gf.getCDSList(i, longTransIndex); /*获得每个基因的cdslist*/
            for (int j = 0; j < cdsl.size(); j++) {
                cdsList.add(cdsl.get(j)); /*cdsList是所有基因的cdslist之和*/
            }
        }
        
        /*把每个基因的cdslist放进cdsList,并按照起始位置排序*/
        Rangescp rs = new Rangescp (cdsList, "cds");
        
        rs.sortByStartPosition();
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if(fs[i].isHidden()){
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles(); //将文件删除后，重新将文件列表打印出来，此时，fs不包含隐藏文件。
        Arrays.sort(fs);
        try {
            BufferedWriter[] bw = new BufferedWriter[types.length];
            for (int i = 0; i < types.length; i++) {
                String outfileS = new File (outfileDirS, types[i]+".txt").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileS);
                bw[i].write("Chr\tPos\tMinorAllele\tMAF\tDerivedAllele\tDAF");
                bw[i].newLine();
            }
            /*对每个文件（每条染色体）进行derived allele 计算， */
            for (int i = 0; i < fs.length; i++) {
                RowTable<String> t1 = new RowTable<> (fs[i].getAbsolutePath());
                int chr = t1.getCellAsInteger(0, 0);
                for (int j = 0; j < t1.getRowNumber(); j++) {
                    int pos = t1.getCellAsInteger(j, 1);
                    if (!rs.isInRanges(chr, pos)) continue; /*判断rs是否是ranges，不是的话，就跳过*/
                    double maf = t1.getCellAsDouble(j, 7);
                    String da = "NA"; /* da = derived allele */
                    String daf = "NA"; /* daf = derived allele frequence */
                    if (t1.getCell(j, 4).length() == 1) { /*是因为有些位点的ancestral allele是NA值，我们将这些位点排除掉*/
                        /*如果ancestral allele存在,且等于major，则da等于第6列的minor, daf 就等于maf
                          如果ancestral allele存在,且等于minor，则da等于第5列的major, daf 就等于 1-maf
                        意思是，如果祖先基因和major相等，derived allele 就等于minor; 否则就等于 major*/
                        String an = t1.getCell(j, 4);
                        if (an.equals(t1.getCell(j, 5))) {
                            da = t1.getCell(j, 6);
                            daf = t1.getCell(j, 7);
                        }
                        else if (an.equals(t1.getCell(j, 6))) {
                            da = t1.getCell(j, 5);
                            daf = String.valueOf(1-Double.valueOf(t1.getCell(j, 7)));
                        }
                    }
                    /*
                    Chr	Pos	MinorAllele	MAF	DerivedAllele	DAF
                    1	111527	A	0.0039893617	G	0.9960106383
                    1	111542	T	0.0012987013	C	0.9987012987
                    */
                    if (t1.getCell(j, 16).equals("SYNONYMOUS")) {
                        if (t1.getCell(j, 17).startsWith("N")) continue;
                        double sift = t1.getCellAsDouble(j, 17);
                        bw[0].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.getCell(j, 6)+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                        bw[0].newLine();
                    }
                    if (t1.getCell(j, 16).equals("NONSYNONYMOUS")) {
                        if (t1.getCell(j, 17).startsWith("N")) continue;
                        double sift = t1.getCellAsDouble(j, 17);
                        if (sift < 0.05) {
                            if (t1.getCell(j, 15).startsWith("N")) continue;
                            if (t1.getCellAsDouble(j, 15) > gerpCut) {
                                bw[3].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.getCell(j, 6)+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                                bw[3].newLine();
                            }
                            bw[2].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.getCell(j, 6)+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                            bw[2].newLine();
                        }
                        else {
                            bw[1].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.getCell(j, 6)+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
                            bw[1].newLine();
                        }
                        
                    }
//                    if (t1.content[j][14].equals("StopLoss")) {
//                        bw[2].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
//                        bw[2].newLine();
//                    }
//                    if (t1.content[j][14].equals("StopGain")) {
//                        bw[3].write(String.valueOf(chr)+"\t"+String.valueOf(pos)+"\t"+t1.content[j][6]+"\t"+String.valueOf(maf)+"\t"+da+"\t"+daf);
//                        bw[3].newLine();
//                    }
                }
            }
            for (int i = 0; i < types.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }

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
            /*这个地方是先过滤数据，将定位在11号12号染色体上的基因过滤掉，并且跳出循环*/
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
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos); // 
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
