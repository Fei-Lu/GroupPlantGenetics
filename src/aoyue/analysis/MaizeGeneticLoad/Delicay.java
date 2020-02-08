/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import com.itextpdf.text.pdf.parser.Path;
import pgl.infra.dna.FastaByte;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import static java.util.Collections.list;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;


/**
 *
 * @author Aoyue
 */
public class Delicay {
    public Delicay() {
        //this.mkDir();
        //this.mkDirs();
        //this.mkHapPosAllele();
        //创建文件并写入字符串
        //this.mkfile1();
        //建立一个具有表格属性的文本
        /*
        调用外部命令
        */
       /*
        将/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt中
        的K16HL0019	Tzi 8	Tzi_8 第1列和第3列信息提取出来，带有表头 测序样本编号	材料名称	Taxa，
        将最后生成的K16HL0019	Tzi_8 文本存放在stanTaxaMap.txt 中。
        */
        //this.stanTaxaMap();
        /*
        将/Users/Aoyue/Documents/mds/1210matrix.txt 的第一列taxa信息提取出来，存放到 1210taxa.txt中
        */
        //this.taxa1210();    
        //this.getTaxaMds();
//        this.rm851mds();
//        this.mds500bp();
        /*
        对TASSEL5产生的genotype_Summary中的TaxaSummary文件进行分组，加标签,851SM和1210SM
        */
//        this.addgroup();
//        this.get1210taxaName();
//        this.mkgroup();
//        this.mkdesTable();
//输出851个 0.0，将FR218和Lan766-4-2这2行的矩阵结果NaN，用 0.0 表示。
        //this.dealwithNaN();
        //删除抽样的VCF中FR218和Lan766-4-2这2列的基因型信息。
        //this.rm2taxa();
        //this.rm2taxa2();
        ///从849taxaIBS矩阵中先挑出100个taxa,然后作树看下情况。
        //this.sampleMatrix();
        //this.get432taxavcf();
        //this.selectVCFByTaxa();
        
        //this.getGBS_reextractDNAID();
        //this.getMissPlate();
        //this.addGeneBiotype();
        //this.parseGff3();
        //this.mkSIFTxlsSimple();
        //this.getVariantType();
        //this.countVariant();
        //this.mergeCountVariantFile();  
        //this.mkSIFTvcfSimple();
        //this.countVariant2();
        //this.mkV4CentromerePos();
        //this.testStringSpilt();
        //this.testmathrandom();
        //this.checkCategory();
        /**
         * 在文件内部过滤数据
         */
        //this.filterGroup();
        /**
         * 根据外部条件过滤数据
         */
        //this.filterGroupUseList();
        /**
         * 查找某列数据的种类及个数,含表头
         */
        this.countCaseInGroup();
        /**
         * 对数据进行排序，打印出来。
         */
        //this.sort();
        //this.testbreak();
        
        
        
        
              
         
    }
    
    public void testbreak(){
        String a = "009";
        int i = Integer.valueOf(a);
        System.out.println(i);

        
    
        
        
    }
    
    
    
            
    
    public void sort(){
        String infileS = "/Users/Aoyue/Documents/国家分布134.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(0);
        Collections.sort(l);
        System.out.println(l.size());
        for(String a : l){
             System.out.println(a);
         }  
        
            /**
             * 遍历hashmap的方法
             */
//            Iterator iter = hm.entrySet().iterator(); 
//            while(iter.hasNext()){
//                Map.Entry entry = (Map.Entry) iter.next(); 
//                Object key = entry.getKey(); 
//                Object val = entry.getValue();
//                if(val.equals("")){
//                    cnt++;
//                    //System.out.println(key);
//                }
//            }
    }
    
    public void countCaseInGroup(){
        //String infileS = "/Users/Aoyue/Documents/国家分布134.txt";
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/005_vcf/018_annoDB/007_checkDelCount/delSNP/chrA.subgenome.delSNP.changeChrPos.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(4);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
        public void filterGroupUseList(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage/phosphoproteome_sortbyGene_addAverage.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_sortbyGene_addAverage_removedupli.txt";
        String tfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_removelist_Dupli40_addP01.txt";
        RowTable<String> t = new RowTable<>(tfileS);
        List<String> l = t.getColumn(0);
        String[] rl = l.toArray(new String[l.size()]);
        Arrays.sort(rl);
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String name = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(rl, name);
                if (index < 0) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void filterGroup(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage/phosphoproteome_sortbyGene_addAverage.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_sortbyGene_addAverage_removedupli.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String group = PStringUtils.fastSplit(temp).get(3);
                if (group.equals("Teo")) {
                    bw.write(temp);
                    bw.newLine();
                }

            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void checkCategory(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/005_loadVSB73Distance/004_reccesiveDeleterious_merge_notraitMixedUnknow.txt";
        RowTable<String> t = new RowTable<>(infileS);
        Set<String> DataSet = new HashSet<>();
        for(int i=0; i < t.getRowNumber(); i++){
            DataSet.add(t.getCellAsString(i, 3));
        }
        String[] datasets = DataSet.toArray(new String[DataSet.size()]);
        System.out.println( DataSet);
        for(int i=0; i< datasets.length; i++){
            try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            Set<String> categorySet = new HashSet<>();
            while((temp = br.readLine()) != null){
                String Dataset = PStringUtils.fastSplit(temp).get(2);
                String category = PStringUtils.fastSplit(temp).get(4);
                if(Dataset.equals(datasets[i])){
                    categorySet.add(category);
                }
            }
            System.out.println(datasets[i] + categorySet);
            br.close();
            
            }
            catch(Exception e){
                //System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            }
        }
    }    
    
    public void testmathrandom(){
        int snpNum = 5875644;
        int size = 10000;
        int[] indices = new int[size]; 
            for (int i = 0; i < size; i++) {
                double a = Math.random();
                indices[i] = (int)(a*snpNum);
            }
            Arrays.sort(indices);
        
    }
    
    public void testStringSpilt(){
        String infileS = "/Users/Aoyue/Documents/SIFTCalculate/01SIFTresult/allresult/hmp321_agpv4_chr10_SIFTannotations.txt";
        String[] spilt = infileS.split("_");
        for(int i = 0; i < spilt.length; i++){
            System.out.println(spilt[i]);
        }
        int chrIndex = Integer.parseInt(infileS.split("_")[2].replaceFirst("chr", ""));
        System.out.println(chrIndex);
        //int chrIndex = Integer.parseInt(f.getName().split("_")[2].replaceFirst("chr", ""))-1;
    }
    
    public void mkV4CentromerePos(){
        this.getchrLength();
        this.countLengthCentromere();
    }

    public void countLengthCentromere(){
        String infile1S = "/Users/Aoyue/Documents/Data/referenceGenome/position/chrlength.txt";
        String infile2S = "/Users/Aoyue/Documents/Data/referenceGenome/position/ChrLenCentPosi_agpV3.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/referenceGenome/position/lengthCompare.txt";
        RowTable t4 = new RowTable(infile1S);
        RowTable t3 = new RowTable(infile2S);
        List<Integer> l = new ArrayList<Integer>();
        for(int i = 0; i< t3.getRowNumber(); i++){
            //int value = l.add(Integer.parseInt(String.valueOf(t4.getCell(i, 3))));
            int value = Integer.parseInt(String.valueOf(t3.getCell(i, 3))) - Integer.parseInt(String.valueOf(t3.getCell(i, 2)));
            int j = i +1;
            System.out.println(value);
            System.out.println("###This is the chr" + j + "centromere length ###");
            
        }
    }
    
    public void getchrLength(){
        String infileS = "/Users/Aoyue/Documents/Data/referenceGenome/AGPv4/maizeAGPv4.fa";
        String outfileS = "/Users/Aoyue/Documents/Data/referenceGenome/position/chrlength.txt";
        FastaByte ref = new FastaByte(infileS);
        int[] chr = new int [10];
        int cnt = 0;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chromosome" + "\t" + "Length(V4)" + "\t" + "CentromereS" + "\t" + "CentromereE");
            bw.newLine();
            
            for(int i=0; i < chr.length; i++){
                int j = i+1;
                chr[i] = ref.getSeqLength(j);
                StringBuilder sb = new StringBuilder(); /*一定要在里面建立sb，否则每次循环写入的东西不会清楚，又重复写到下一行*/
                sb.append(j).append("\t").append(chr[i]).append("\t");
                bw.write(sb.toString());
                bw.newLine(); 
                cnt = cnt + chr[i];
                System.out.println(chr[i]);
                System.out.println("###This is the " + j + "chr length ###");
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("This is the code result " + ref.getTotalSeqLength());
        System.out.println("This is calculated result " + ref.getTotalSeqLength());
        /*This is the code result 2135098301
            This is calculated result 2135098301*/
    }
    
    public void countVariant2(){
        String infileDirS = "/Volumes/LuLab3T_30/allresult/hmp321_agpv4_chr7_SIFTpredictions.vcf";
        String outfileDirS = "/Users/Aoyue/Documents/SIFTCalculate/countVariant/hmp321_agpv4_chr7_SIFTpredictions.txt";
        BufferedReader br = IOUtils.getTextReader(infileDirS);
        BufferedWriter bw = IOUtils.getTextWriter(outfileDirS);
        String temp = null;
        int cntSynonymous = 0;
        int cntNonsynonymous = 0;
        int cntStartlost = 0;
        int cntDeleterious = 0;
        int cntTolerent = 0;
        try{
            List<String> l = null;
                List<String> infoList = null;
                List<String> siftInfoList = null;
                String siftInfo = null;
                while ((temp = br.readLine()) != null) {
                    if(temp.startsWith("#")) continue;
                    l = PStringUtils.fastSplit(temp);
                    String info = l.get(7);
                    infoList = PStringUtils.fastSplit(info, ";");
                    /*将List转化为数组，如果数组等于中的元素以SIFTINFO开头，将此元素转化siftInfoList为集合*/
                    //String siInfo = infoList.
                   String[] infoArray = infoList.toArray(new String[infoList.size()]);                   
                   for (int i = 0; i < infoArray.length; i++){
                       if(infoArray[i].startsWith("SIFTINFO")){
                           siftInfo = infoArray[i];
                           siftInfoList = PStringUtils.fastSplit(siftInfo, "|");
                           /*siftInfoList一共包含13个元素，我们只关心第6列的变异类型和第13列的SIFTScore值*/
                           //if(siftInfoList.get(5).equals("NA")) continue; 
                           if(siftInfoList.get(5).equals("SYNONYMOUS")) cntSynonymous++; 
                           if(siftInfoList.get(5).equals("NONSYNONYMOUS")) {
                               cntNonsynonymous++;
                               if(siftInfoList.get(8).equals("NA")) continue;
                               if(Float.parseFloat(siftInfoList.get(8))<0.05){
                                   cntDeleterious++;                                
                               }
                               else{
                                   cntTolerent++;                                  
                               }
                           }
                           if(siftInfoList.get(5).equals("cntStartlost")) cntStartlost++;                       
                       }                 
                   }  
                }
                
                StringBuilder sb = new StringBuilder();
                sb.append("CHROM").append("\t").append("NONSynonymous").append("\t").append("SIFT_deleterious_mutation")
                        .append("\t").append("Non-synonymous_tolerent_mutation")
                        .append("\t").append("Neutral_mutation").append("\t").append("Startlost_mutation").append("\n");
               
                sb.append("hmp321_agpv4_chr7").append("\t").append(cntNonsynonymous).append("\t")
                        .append(cntDeleterious).append("\t").append(cntTolerent)
                        .append("\t").append(cntSynonymous).append("\t").append(cntStartlost);
                bw.write(sb.toString());
                bw.newLine();
                
                System.out.println("###### " + "hmp321_agpv4_chr7" + " ######result:");
                System.out.println("The number of Neutral mutation (Synonymous) is: " + cntSynonymous);
                System.out.println("The number of NONSynonymous is: " + cntNonsynonymous);
                System.out.println("The number of SIFT_deleterious_mutation is: " + cntDeleterious);
                System.out.println("The number of Non-synonymous_tolerent_mutation is: " + cntTolerent);
                System.out.println("The number of Startlost is: " + cntStartlost);
                
                bw.flush();
                bw.close();
                br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }               
    }
    
    public void mkSIFTvcfSimple(){
        String infileDirS = "/data1/home/aoyue/sift/output/allresult";
        String outfileDirS = "/data1/home/aoyue/sift/output/result_SIFTvcfSimple";
        new File (outfileDirS).mkdir();
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", "_simple.vcf")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    if(temp.startsWith("##")) continue;
                    temp = temp.substring(0, 420);
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t").append(l.get(7));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt%100000 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            } 
        });       
    }
    
     public void mergeCountVariantFile(){
        String infileDirS = "/Users/Aoyue/Documents/SIFTCalculate/ffffffff";
        String outfileS = "/Users/Aoyue/Documents/SIFTCalculate/allSIFTpredictions.txt";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        try{                       
            StringBuilder sb = new StringBuilder();
            sb.append("CHROM").append("\t").append("NONSynonymous").append("\t").append("SIFT_deleterious_mutation")
                        .append("\t").append("Non-synonymous_tolerent_mutation")
                        .append("\t").append("Neutral_mutation").append("\t").append("Startlost_mutation").append("\n");
            for(int i = 0; i < fs.length; i++){
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                String subFileS = fs[i].getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(subFileS);
                String temp = null;
                while((temp = br.readLine()) != null){
                    if(temp.startsWith("CHROM")) continue;
                    sb.append(temp).append("\n");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();  
            }
        }
        catch(Exception e){
            System.out.println();
            e.printStackTrace();
            System.exit(1);              
        }
        
        RowTable t = new RowTable("/Users/Aoyue/Documents/SIFTCalculate/allSIFTpredictions.txt");
        
        
        
//        try{
//            String total = "/Users/Aoyue/Documents/SIFTCalculate/totalSIFTpredictions.txt";
//            BufferedReader br = IOUtils.getTextReader(outfileS);
//            BufferedWriter bw = IOUtils.getTextWriter(total);
//            int cntSynonymous = 0;
//            int cntNonsynonymous = 0;
//            int cntStartlost = 0;
//            int cntDeleterious = 0;
//            int cntTolerent = 0;
//            String temp;
//            List<String> l = null;         
//            while((temp = br.readLine()) != null){
//                if(temp.startsWith("CHROM")) continue;
//                l = PStringUtils.fastSplit(temp);
//                cntNonsynonymous = cntNonsynonymous + Integer.parseInt(l.get(1));
//                cntDeleterious = cntDeleterious + Integer.parseInt(l.get(2));
//                cntTolerent = cntTolerent + Integer.parseInt(l.get(3));
//                cntSynonymous = cntSynonymous + Integer.parseInt(l.get(4));
//                cntStartlost = cntStartlost + Integer.parseInt(l.get(5));                            
//            }
//            StringBuilder sb = new StringBuilder();
//            sb.append("Total").append("\t").append(String.valueOf(cntNonsynonymous)).append("\t").
//                    append(String.valueOf(cntDeleterious)).append("\t").append(String.valueOf(cntTolerent)).append("\t")
//                    .append(String.valueOf(cntSynonymous)).append("\t").append(String.valueOf(cntStartlost)).append("\n");
//            bw.write(sb.toString());
//            bw.newLine();
//            bw.flush();
//            bw.close();
//            br.close();
//        }
//        catch(Exception e){
//            System.out.println();
//            e.printStackTrace();
//            System.exit(1); 
//            
//        }
    }
    
    
    public void countVariant () {
        String infileDirS = "/Users/Aoyue/Documents/SIFTCalculate/result_SIFTvcfSimple";
        String outfileDirS = "/Users/Aoyue/Documents/SIFTCalculate/02Plot_SIFT/countVariant";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS,f.getName().replaceFirst("_SIFTpredictions.vcf", "_SIFTpredictions.txt")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cntSynonymous = 0;
            int cntNonsynonymous = 0;
            int cntStartlost = 0;
            int cntDeleterious = 0;
            int cntTolerent = 0;
            try{
                List<String> l = null;
                List<String> infoList = null;
                List<String> siftInfoList = null;
                String siftInfo = null;
                while ((temp = br.readLine()) != null) {
                    if(temp.startsWith("#")) continue;
                    l = PStringUtils.fastSplit(temp);
                    String info = l.get(7);
                    infoList = PStringUtils.fastSplit(info, ";");
                    /*将List转化为数组，如果数组等于中的元素以SIFTINFO开头，将此元素转化siftInfoList为集合*/
                   String[] infoArray = infoList.toArray(new String[infoList.size()]);
                   
                   for (int i = 0; i < infoArray.length; i++){
                       if(infoArray[i].startsWith("SIFTINFO")){
                           siftInfo = infoArray[i];
                           siftInfoList = PStringUtils.fastSplit(siftInfo, "|");
                           /*siftInfoList一共包含13个元素，我们只关心第6列的变异类型和第13列的SIFTScore值*/
                           //if(siftInfoList.get(5).equals("NA")) continue; 
                           if(siftInfoList.get(5).equals("SYNONYMOUS")) cntSynonymous++; 
                           if(siftInfoList.get(5).equals("NONSYNONYMOUS")) {
                               cntNonsynonymous++;
                               if(siftInfoList.get(8).equals("NA")) continue;
                               if(Float.parseFloat(siftInfoList.get(8))<0.05){
                                   cntDeleterious++;                                
                               }
                               else{
                                   cntTolerent++;                                  
                               }
                           }
                           if(siftInfoList.get(5).equals("cntStartlost")) cntStartlost++;                       
                       }                 
                   }  
                }
                
                StringBuilder sb = new StringBuilder();
                sb.append("CHROM").append("\t").append("NONSynonymous").append("\t").append("SIFT_deleterious_mutation")
                        .append("\t").append("Non-synonymous_tolerent_mutation")
                        .append("\t").append("Neutral_mutation").append("\t").append("Startlost_mutation").append("\n");
               
                sb.append(f.getName().replaceFirst("_SIFTpredictions.vcf", "")).append("\t").append(cntNonsynonymous).append("\t")
                        .append(cntDeleterious).append("\t").append(cntTolerent)
                        .append("\t").append(cntSynonymous).append("\t").append(cntStartlost);
                bw.write(sb.toString());
                bw.newLine();
                
                System.out.println("###### " + f.getName().replaceFirst("_SIFTannotations_simple.txt", "") + " ######result:");
                System.out.println("The number of Neutral mutation (Synonymous) is: " + cntSynonymous);
                System.out.println("The number of NONSynonymous is: " + cntNonsynonymous);
                System.out.println("The number of SIFT_deleterious_mutation is: " + cntDeleterious);
                System.out.println("The number of Non-synonymous_tolerent_mutation is: " + cntTolerent);
                System.out.println("The number of Startlost is: " + cntStartlost);
                
                bw.flush();
                bw.close();
                br.close();
       
            }
           
            catch(Exception e){
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);                
            }        
            //System.out.println("Non-synonymous: " + cntNonsynonymous);
           // System.out.println("SIFT deleterious mutation (Non-synonymous, SIFT < 0.05)" + "is" + cntNonsynonymous);          
        });
    }
        
    public void getVariantType(){
        String infileDirS = "/Users/Aoyue/Documents/SIFTCalculate/allsimpleresult";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder();
                List<String> l = null;
                Set<String> variantType = new HashSet();
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(5));
                    cnt++;
                    variantType.add(l.get(5)); 
                }
                br.close();
                int countSNP = cnt-1;
                System.out.println(String.valueOf(countSNP) + " SNPs output from " + f.getAbsolutePath());
                System.out.println(variantType);
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            } 
        });      
    }
    
    public void mkSIFTxlsSimple () {
        String infileDirS = "/Users/Aoyue/Documents/SIFTCalculate/allresult";
        String outfileDirS = "/Users/Aoyue/Documents/SIFTCalculate/allsimpleresult";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_simple.txt")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder();
                List<String> l = null;
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp);
                    sb = new StringBuilder(l.get(0));
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").append(l.get(4))
                            .append("\t").append(l.get(8)).append("\t").append(l.get(12)).append("\t").append(l.get(13)).append("\t").append(l.get(16));
                    bw.write(sb.toString());
                    bw.newLine();
                    if (cnt%100000 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath());
                System.out.println(temp);
                e.printStackTrace();
                System.exit(1);
            } 
        });
    }
        
    public void parseGff3(){
        String InputFileS = "/Users/Aoyue/Documents/SIFTCalculate/gene_annotation/Zea_mays.AGPv4.38.gtf";
        String[] OutFileS = new String[10];
        //String OutmtFileS = "";
        //"/Users/Aoyue/Documents/Zea_mays.AGPv4.38.gtf";
        for(int i = 0; i < 10; i++){
            int j = i+1;
            OutFileS[i] = "/Users/Aoyue/Documents/Zea_mays.AGPv4.38_" + j + ".gtf";
            System.out.println(OutFileS[i] + "\n");        
        }
        try{
            for(int i = 0; i < 10; i++){
                BufferedReader br = IOUtils.getTextReader(InputFileS);
                BufferedWriter bw = IOUtils.getTextWriter(OutFileS[i]);
                String tem = null;
                List<String> l = null;
                while((tem = br.readLine()) != null){
                     l = PStringUtils.fastSplit(tem);
                     int chr = i + 1;
                     String j = Integer.toString(chr);
                     if (j.equals(l.get(0))){
                         StringBuilder sb = new StringBuilder();
                         sb.append(tem);
                         bw.write(sb.toString());
                         bw.newLine();                         
                     }
                }
                bw.flush();
                bw.close();
                br.close();                
            }
            BufferedReader br = IOUtils.getTextReader(InputFileS);
            //BufferedWriter bw = IOUtils.getTextWriter(OutmtFileS);
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);            
        }        
    }
    
    public void addGeneBiotype(){
        String InputFileS = "/Users/Aoyue/Documents/SIFTCalculate/gene_annotation/Zea_mays.AGPv4.40.gene.gtf";
        String OutputFileS = "/Users/Aoyue/Documents/SIFTCalculate/gene_annotation/Zea_mays.AGPv4.40.gtf";
        try{
            BufferedReader br = IOUtils.getTextReader(InputFileS); 
            BufferedWriter bw = IOUtils.getTextWriter(OutputFileS);
            String temp = null;
            while ((temp = br.readLine()) != null){         
                    //String plate = null;
                    //plate = PStringUtils.fastSplit(temp).get(2);
                    //if(plate.equalsIgnoreCase("exon") | plate.equalsIgnoreCase("CDS") | plate.equalsIgnoreCase("stop_codon") | plate.equalsIgnoreCase("start_codon")){
                    List<String> l = null;
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3)).append("\t").
                            append(l.get(4)).append("\t").append(l.get(5)).append("\t").
                                    append(l.get(6)).append("\t").append(l.get(7)).append("\t").append(l.get(8)).append(" ")
                                            .append("gene_biotype").append(" ").append("\"protein_coding\";");
                    bw.write(sb.toString());
                    bw.newLine();         
            }
            bw.flush();
            bw.close();
            br.close();                 
            }           
            catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        }              
    
    public void getMissPlate(){
        String OriginalInfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/RerxtractDNA_Plate.txt";
        String outfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/Re_extractDNA_ID.txt";
        RowTable<String> oT = new RowTable<>(OriginalInfileS);
        RowTable<String> dT = new RowTable<>(outfileS);
        Set<String> plateSet = new HashSet (oT.getColumn(2));
        Set<String> IDSet = new HashSet (dT.getColumn(1));
        try{
            BufferedReader br = IOUtils.getTextReader(OriginalInfileS);           
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("S")) continue;             
                else{
                    
                    String tem = temp;
                    String plate = null;
                    plate = PStringUtils.fastSplit(tem).get(2);
                    if(IDSet.add(plate)){
                        System.out.println (plate); 
                    }            
                }    
            }                
                 br.close();
            }           
            catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        }                          
    
    public void getGBS_reextractDNAID(){
        String PlateinfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/Re2extractDNA_Plate.txt";
        String OriginalInfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/GBS_ID2.txt";
        String outfileS = "/Users/Aoyue/Documents/IGDB/4-PHD/wheat/GBS取样/实验记录/Re2-2extractDNA_ID.txt";
        RowTable<String> t = new RowTable<>(PlateinfileS);
        List<String> plateList = t.getColumn(2);
        Collections.sort(plateList);
        try {
            BufferedReader br = IOUtils.getTextReader(OriginalInfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("SerialNumber\tPlate\tID");
            bw.newLine();
            String temp = br.readLine();
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (Collections.binarySearch(plateList, tem[1]) < 0) continue;
                cnt++;
                bw.write(String.valueOf(cnt)+"\t"+temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    public void selectVCFByTaxa () {
        String inputFileS = "/Users/Aoyue/Documents/tree/backup/chr10000.vcf";
        String taxaFileS = "/Users/Aoyue/Documents/tree/432taxa/432taxa.tab.txt";
        String outfileS = "/Users/Aoyue/Documents/tree/432taxa/chr10000_432taxa.vcf";
        RowTable<String> t = new RowTable<>(taxaFileS);
        List<String> taxaList = t.getColumn(0);
        Collections.sort(taxaList);
        TIntArrayList indexList = new TIntArrayList();
        int[] indices = null;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    bw.write(temp);
                    bw.newLine();
                }
                else if (temp.startsWith("#C")) {
                    String[] tem = temp.split("\t");
                    for (int i = 0; i < 9; i++) {
                        indexList.add(i);
                    }
                    for (int i = 9; i < tem.length; i++) {
                        if (Collections.binarySearch(taxaList, tem[i]) < 0) continue;
                        indexList.add(i);
                    }
                    indices = indexList.toArray();
                    sb = new StringBuilder();
                    for (int i = 0; i < indices.length; i++) {
                        sb.append(tem[indices[i]]).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                else {
                    String[] tem = temp.split("\t");
                    sb = new StringBuilder();
                    for (int i = 0; i < indices.length; i++) {
                        sb.append(tem[indices[i]]).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void get432taxavcf(){
        String infileS = "/Users/Aoyue/Documents/tree/432taxa/432taxa.txt";
        String infile2S = "/Users/Aoyue/Documents/tree/1wansite/chr10000_rm2taxa.vcf";
        String outfileS = "/Users/Aoyue/Documents/tree/432taxa/chr10000_432taxa_column.txt";
        String outfile2S = "/Users/Aoyue/Documents/tree/432taxa/chr10000_432taxabyyue.vcf";
        int[] indices = null;
        try{
            BufferedReader br = IOUtils.getTextReader(infileS); 
            BufferedReader br2 = IOUtils.getTextReader(infile2S);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            String line = null;
            String lines = null;
            int k = 0;            
            Set<String> taxaset = new HashSet<>();
            while ((line = br.readLine()) != null){
                taxaset.add(line);
                k++;             
            }
            System.out.println ("检验 " + k + " 个taxa name是否放入set中"); 
           
//            Iterator<String> itr = taxaset.iterator();
//            while(itr.hasNext()){
//                System.out.println (itr.next());    
//            }
            List<String> l = null;
            List<String> colume = null;
            while((lines =br2.readLine()) != null){
                if (lines.startsWith("##")) continue; 
                if (lines.startsWith("#CHROM")){
                    l = PStringUtils.fastSplit(lines);
                    System.out.println(l.size() + " 查看每行文本中的列数");
                    int n = 0;
                    for (int i = 0; i < l.size(); i++){
                        if (taxaset.add(l.get(i))) continue;
                        else {
                            n++;
                            sb.append(i).append("\n"); 
                            
                        }                       
                    }
                    bw.write(sb.toString());
                    System.out.println(n + " 统计每行文本中432列数信息是否提取出来");
                }         
            }
            bw.flush();
            bw.close();
            br.close();
            br2.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);           
        }
        try{
            BufferedReader br = IOUtils.getTextReader(outfileS);
            BufferedReader br2 = IOUtils.getTextReader(infile2S);
            BufferedWriter bw = IOUtils.getTextWriter(outfile2S);
            StringBuilder sb = new StringBuilder();
            String line = null;
            String lines = null;
            int k = 0; 
            TIntArrayList columnlist = new TIntArrayList();          
            columnlist.add(0);
            columnlist.add(1);
            columnlist.add(2);
            columnlist.add(3);
            columnlist.add(4);
            columnlist.add(5);
            columnlist.add(6);
            columnlist.add(7);
            columnlist.add(8);          
            while ((line = br.readLine()) != null){
                columnlist.add(Integer.parseInt(line));
                k++;             
            }
            System.out.println ("k = " + k + " 把432列数放入一个list");
            indices = columnlist.toArray();
//            for (int i = 0; i < columnset.size(); i++){
//                System.out.print(columnset.get(i) + " ");
//                
//            }            
            List<String> l = null;           
            while((lines =br2.readLine()) != null){
                if (lines.startsWith("##")) {
                    bw.write(lines);
                    bw.newLine();
                }
                //else {
                    //if (lines.startsWith("#CHROM")){
                    
                    l = PStringUtils.fastSplit(lines); 
                    //int n = 0;
                    sb = new StringBuilder();
                    for (int i = 0; i < indices.length; i++){
                        //n++ ;
                         
                        //int a = Integer.parseInt(columnlist.get(i));                      
                        
                        sb.append(l.get(indices[i])).append("\t");                            
                    } 
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                    //System.out.println(n);                    
                //}                                    
                //}
                
            }           
            bw.flush();
            bw.close();
            br.close();
            br2.close();           
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);            
        }       
    }
    
    ///没写完
    public void sampleMatrix(){
        String infileS = "/Users/Aoyue/Documents/matrix.txt";
        String outfileS = "/Users/Aoyue/Documents/matrix100.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> l1 = null;
            List<String> l2 = null;
            //List<Integer> removeIndexList = new ArrayList<Integer>();
            for (int i = 0; i < 101; i++){
                if (temp.startsWith("849")) continue; 
                
                                   
                    l2 = PStringUtils.fastSplit(temp);
                 
//                    for (int i = 0; i < removeIndexList.size(); i++) {
//                        l2.remove((int)removeIndexList.get(i));
//                    }      
                    StringBuilder sb = new StringBuilder();
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();                                                    
            } 
            bw.flush();
            bw.close();
            br.close();
            }      
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }   
        
    }
    ///没写完
    public void rm2taxa2(){
        String infileS = "/Users/Aoyue/Documents/bcftools_viewTaxa/tabletest.vcf";
        String outfileS = "/Users/Aoyue/Documents/bcftools_viewTaxa/out.vcf";
        RowTable t = new RowTable("infileS");
        t.removeColumn(269);
        t.removeColumn(456);
        
        
        
        
        
    }
    
    public void rm2taxa(){
        String infileS = "/Users/Aoyue/Documents/tree/backup/chr10000.vcf";
        String outfileS = "/Users/Aoyue/Documents/tree/backup/chr10000_rm2taxa.vcf";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            List<String> l1 = null;
            List<String> l2 = null;
            //List<Integer> removeIndexList = new ArrayList<Integer>();
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("##")) continue; 
                if (temp.startsWith("#CHROM")){
                    l1 = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    for (int i =0 ; i < l1.size(); i++){
                        if (l1.get(i).equals("FR218") || l1.get(i).equals("Lan766-4-2") ){
                            System.out.println(i + "\t" + "this is the wrong data column");
                            l1.remove(i);
                            //removeIndexList.add(i);
                            i--;
                        }
                        else {
                            sb.append(l1.get(i)).append("\t");
                        }                    
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();  
                    System.out.println(l1.size() + " " + "biaotou");
                 
                }
                else{                   
                    l2 = PStringUtils.fastSplit(temp);
                    l2.remove(269);
                    l2.remove(456);
//                    for (int i = 0; i < removeIndexList.size(); i++) {
//                        l2.remove((int)removeIndexList.get(i));
//                    }      
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < l2.size(); i++){    
                        sb.append(l2.get(i)).append("\t");                                                
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();                                       
                } 
                
            } 
            bw.flush();
            bw.close();
            br.close();
            }      
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                                
    }
    

    public void dealwithNaN(){
        String infileS = "/Users/Aoyue/Documents/nan.txt";
        String outfileS = "/Users/Aoyue/Documents/tree/chr5000matrix_rmNaN.txt";
//        RowTable t = new RowTable(infileS);
//        for (int i = 0 ; i < t.getRowNumber(); i++){
//            t.getRow(i);            
//            System.out.println(t.getRow(i));            
//        }
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
               
                    String tem = temp;
                    List<String> l = null;
                    l = PStringUtils.fastSplit(tem);
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 851; i++){
                        sb.append(0.0).append("\t");               
                    }
                    bw.append(sb.toString());
                    bw.newLine();               
                       } 
            bw.flush();
            bw.close();
            br.close();
            }      
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                            
    }
    
    public void mkdesTable(){
        String infileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_10000TaxaSummary_grouped.txt";
        String outfileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_heter.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
               
                    String tem = temp;
                    List<String> l = null;
                    l = PStringUtils.fastSplit(tem);
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(1)).append("\t").append(l.get(4)).append("\t").append(l.get(6)).append("\t").append(l.get(9));
                    bw.append(sb.toString());
                    bw.newLine();               
                       } 
            bw.flush();
            bw.close();
            br.close();
            }      
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                          
    }
    
    public void mkgroup(){
        String taxaS = "/Users/Aoyue/Documents/genotype_summary/1210taxaName.txt";
        String infileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_10000TaxaSummary_addgroup.txt";
        String outfileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_10000TaxaSummary_grouped.txt";
        Set<String> taxaSet = new HashSet();  
        try{
            BufferedReader br = IOUtils.getTextReader(taxaS);                 
            String temp = null;
            while ((temp = br.readLine()) != null){
                String tem = temp;
                taxaSet.add(tem);                          
            }          
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("Taxa")) {
                    String tem = temp;
                    StringBuilder sb = new StringBuilder();
                    sb.append(tem);
                    bw.append(sb.toString());
                    bw.newLine();    
                }
                else{
                    String tem = temp;
                    String taxa = null;
                    taxa = PStringUtils.fastSplit(tem).get(1);
                    if(taxaSet.add(taxa)){
                        StringBuilder sb = new StringBuilder();
                        sb.append(tem);
                        bw.append(sb.toString());
                        bw.newLine();         
                    }
                    else{
                        List<String> l = null;
                        l = PStringUtils.fastSplit(tem);
                        StringBuilder sb = new StringBuilder();
                        sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).
                                append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t").
                                append(l.get(5)).append("\t").append(l.get(6)).append("\t").append(l.get(7)).
                                append("\t").append(l.get(8)).append("\t").append(l.get(9).replaceFirst("851SM", "1210SM"));
                        bw.write(sb.toString());
                        bw.newLine();                    
                    }             
                }        
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                     
    }
    
    public void get1210taxaName(){
        String infileS = "/Users/Aoyue/Documents/genotype_summary/chrlibrary_100001TaxaSummary.txt";
        String outfileS = "/Users/Aoyue/Documents/genotype_summary/1210taxaName.txt";
        List<String> taxaList = new ArrayList <>();
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("Taxa")) continue;
                String tem = temp;
                List<String> l = null;
                l = PStringUtils.fastSplit(tem);
                StringBuilder sb = new StringBuilder(l.get(1));
                bw.write(sb.toString());
                bw.newLine();        
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }               
    }
    
    public void addgroup(){
        String infileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_10000TaxaSummary.txt";
        String outfileS = "/Users/Aoyue/Documents/genotype_summary/chrmerge_10000TaxaSummary_addgroup.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                String tem = temp;
                if (temp.startsWith("Taxa")) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(tem).append("\t").append("Group");
                    bw.write(sb.toString());
                    bw.newLine();     
                }
                else{
                    StringBuilder sb = new StringBuilder();
                    sb.append(tem).append("\t").append("851SM");
                    bw.write(sb.toString());
                    bw.newLine();                    
                }                       
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                       
    }
    
    public void mds500bp(){
        String infileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/434sample.txt";
        String infile2S = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/stanTaxaMap.txt";
        String infile3S = "/Users/Aoyue/Documents/mds/851orimds.txt";
        String outfileS = "/Users/Aoyue/Documents/mds/434taxa.txt";
        String outfile2S = "/Users/Aoyue/Documents/mds/434mds.txt";
        HashMap<String, String> taxaMap = new HashMap<>();    
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedReader br2 = IOUtils.getTextReader(infile2S);
        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
        BufferedWriter bw2 = IOUtils.getTextWriter(outfile2S);
        try{          
            String temp = null;
            String line = null;
            List<String> l = null;
            while ((temp = br2.readLine()) != null){
                l = PStringUtils.fastSplit(temp); ///将此行的temp放进l list中
                taxaMap.put(l.get(0), l.get(1)); /// 将第1列测序ID当作key，第2列taxa信息当作value.                                                     
            }
            while ((line = br.readLine()) != null){ ///把434sample的测序ID读进去
                StringBuilder sb = new StringBuilder();
                sb.append(taxaMap.get(line)); ///得到434对应的taxa信息
                bw.write(sb.toString());
                bw.newLine();
            }          
            bw.flush();
            bw.close();
            br2.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }    
        Set<String> infile = new HashSet();  
        try{
            BufferedReader br3 = IOUtils.getTextReader(outfileS); ///将434taxa读进去                
            String temp = null;
            while ((temp = br3.readLine()) != null){ ///将434taxa放入hashset中
                infile.add(temp);                          
            } 
            String line2 = null;
            BufferedReader br4 = IOUtils.getTextReader(infile3S); 
            while ((line2 = br4.readLine()) != null){
                String taxa = null;
                taxa = PStringUtils.fastSplit(line2).get(0); ///取每一行的taxa信息
                if(infile.add(taxa)) continue;///
                else{
                    StringBuilder sb = new StringBuilder();
                    sb.append(line2);
                    bw2.write(sb.toString());
                    bw2.newLine();                    
                }                
            }           
            bw2.flush();
            bw2.close();
            br3.close();
            br4.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }             
    }
    
    public void rm851mds(){
        String infileS = "/Users/Aoyue/Documents/mds/851mds.txt";
        String outfileS = "/Users/Aoyue/Documents/mds/851oriMds.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                String tem = temp;
                StringBuilder sb = new StringBuilder();
                sb.append(PStringUtils.fastSplit(tem).get(0).replaceFirst("2:", "")).append("\t").append(PStringUtils.fastSplit(tem).get(1)).append("\t")
                    .append(PStringUtils.fastSplit(tem).get(2)).append("\t").append(PStringUtils.fastSplit(tem).get(3));
                bw.write(sb.toString());
                bw.newLine();                                        
            }          
            br.close();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }                      
    }
    
    public void getTaxaMds(){
        String orifileS = "/Users/Aoyue/Documents/mds/2061matrix.txt";
        String infileS = "/Users/Aoyue/Documents/mds/1210taxa.txt";
        String outfileS = "/Users/Aoyue/Documents/mds/1210mds.txt";
        String outfile2S = "/Users/Aoyue/Documents/mds/851mds.txt";
        Set<String> infile = new HashSet();  
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);                 
            String temp = null;
            while ((temp = br.readLine()) != null){
                String tem = temp;
                infile.add(tem);                          
            }          
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }              
        try{
            BufferedReader br = IOUtils.getTextReader(orifileS);     
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfile2S);
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("\t")) continue;
                String tem = temp;
                String taxa = null;
                taxa = PStringUtils.fastSplit(tem).get(0);
                if(infile.add(taxa)) {
                   //System.out.println(tem);
                   StringBuilder sb2 = new StringBuilder();
                   sb2.append(tem);
                   bw2.write(sb2.toString());
                   bw2.newLine();                    
                }
                else{
                    StringBuilder sb = new StringBuilder();
                    sb.append(tem);
                    bw.write(sb.toString());
                    bw.newLine();                    
                }                
            }
            bw2.flush();
            bw2.close();
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }              
    }
    
    public void taxa1210(){
        String infileS = "/Users/Aoyue/Documents/mds/1210matrix.txt";
        String outfileS = "/Users/Aoyue/Documents/mds/1210taxa.txt";
        RowTable <String> t = new RowTable<>(infileS);
        List<String> taxaList = new ArrayList <>();
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("\t")) continue;
                String tem = temp;
                List<String> l = null;
                l = PStringUtils.fastSplit(tem);
                StringBuilder sb = new StringBuilder(l.get(0));
                bw.write(sb.toString());
                bw.newLine();        
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }             
    }
    
    public void stanTaxaMap(){
        String infileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/stanTaxaMap.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                if (temp.startsWith("测")) continue;
                String tem = temp;
                List<String> l = null;
                l = PStringUtils.fastSplit(tem);
                StringBuilder sb = new StringBuilder(l.get(0)).append("\t").append(l.get(2));
                bw.write(sb.toString());
                bw.newLine();        
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }       
    }
    
    
    public void mkfile1(){
        String input = "this is the first time that I have write a txt using java !";
        byte[] sourceByte = input.getBytes();
        if (null != sourceByte){
            try {
                File f = new File ("/Users/Aoyue/Documents/afile/a.txt");
                if (!f.exists()){
                    File dir = new File (f.getParent());
                    dir.mkdirs();
                    f.createNewFile();       
                }
                FileOutputStream outStream = new FileOutputStream(f);
                outStream.write(sourceByte);
                outStream.close();              
            }
            catch (Exception e){
                e.printStackTrace();        
            }           
        }     
    }
    
    public void mkHapPosAllele () {
        String infileDirS = "/Users/Aoyue/Desktop/";
        String outfileDirS = "/Users/Aoyue/Desktop/out/"; //该路径必须提前设置好
        //public String[] list() Returns an array of strings naming the files and directories in the directory denoted by this abstract pathname.
        //返回由此抽象路径名所表示的目录中的文件和目录的名称所组成字符串数组。
        //public File[] listFiles() Returns an array of abstract pathnames denoting the files in the directory denoted by this abstract pathname.
        //返回一个抽象路径名数组，这些路径名表示此抽象路径名所表示目录中的文件。
        File fsall = new File(infileDirS);
        File[] fs = fsall.listFiles();
        //String [] fs1 = new File (infileDirS).list(); //String类型
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs); //public static <T> List<T> asList(T... a) {} 把以vcf.gz结尾的文件放入fsList中
        fsList.parallelStream().forEach(f -> { //forEach流是Stream接口里的方法，对单个文件vcf.gz进行操作。
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath()); //把vcf.gz解压，将文件读进缓冲字符流中去
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".posAllele.txt.gz")).getAbsolutePath();//输出文件名设置
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //将文件写出来
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder(); //new一个新的StringBuilder开始写文件（一个可变的字符序列）
                sb.append("Chr\tPos\tRef\tAlt"); //开始建立表头 chr Pos Alt
                bw.write(sb.toString()); //Writer类 是抽象类，用于写出字符流。
                bw.newLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) { //把当vcf.gz按行读不为空时，
                    if (temp.startsWith("#")) continue; //public boolean startsWith(String prefix) 以#开头时，继续读下一行
                    temp = temp.substring(0, 40); // 取每行的前0-40 字符串
                    l = PStringUtils.fastSplit(temp); //l最开始为一个空的list, 后将temp以 \t 分割后，放入l中。
                    sb = new StringBuilder(l.get(0)); //取集合List的第0列，即染色体号。上文已讲StringBuilder是一个可变字符序列
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)); // 将文件第0 1 3 4列的信息提取出来
                    bw.write(sb.toString()); 
                    bw.newLine();
                    if (cnt%2 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush(); //将字符输出流冲出来
                bw.close(); //将字符输出流关闭
                br.close(); //将字符输入流关闭
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath()); //如果文件有异常，就将文件路径打印出来
                System.out.println(temp); 
                e.printStackTrace();
                System.exit(1);
            }
            
        });
    }
    
    public static File[] listFilesEndsWith (File[] all, String endStr) {
        ArrayList<File> al = new ArrayList(); //建立一个新的集合类 al， 并且规定该类是 File类型.
        //将以endStr结尾的文件添加到集合al 中
        for (int i = 0; i < all.length; i++) { 
            if (all[i].getName().endsWith(endStr)) al.add(all[i]);
        }
        return al.toArray(new File[al.size()]); //size()返回集合元素的数量; new File[3]; toArray返回T[] a（即返回集合中的元素）
    }
    
    public void mkDirs(){
        String dirdemo = "/Users/Aoyue/Documents/SIFTCalculate/outputbatches";
        for(int i=1; i< 11; i++ ){
            String dir = dirdemo.replaceFirst("/Users/Aoyue/Documents/SIFTCalculate/outputbatches", "/Users/Aoyue/Documents/SIFTCalculate/outputbatches" + i);         
            File f = new File (dir);
            if (!f.exists()) {  
            f.mkdir();  
            }
            System.out.println("Perform the end "+ dir);          
        }
        
    }
    
    public void mkDir (){
        String dir = "/Users/Aoyue/Documents/SIFTCalculate";
        File f = new File (dir); 
        if (!f.exists()) {  
            f.mkdir();  
        }  
        System.out.println("Perform the end "+ dir);  
    }  
    public static void main (String[] args) {        
        new Delicay(); 
        System.out.println("I can do it !");    
        
        
    }
          
} 
     
