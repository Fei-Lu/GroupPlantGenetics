/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.graph.r.DensityPlot;
import pgl.graph.r.Histogram;
import pgl.graph.r.ScatterPlot;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

import static pgl.infra.utils.IOFileFormat.Text;

/**
 *
 * @author Aoyue
 */
public class VariantSummary {
    
    public VariantSummary(){
        this.countSite();
        this.subSetTest();
        //this.densityTest_deprecated();
        //this.densityTest();
        //this.density();
        this.filterHmp321Info();
        this.filterHmp321Info_bysiftTrans();
        this.filterHmp321Info_bysiftTrans_useDataBase();
        //this.testarraylength();
//        this.summarizeTranscript_deprecated();
        this.summarizeTranscript2();
       this.classifySNPs();
//       this.binarySearchtest();
       this.mkBarplotOfSNPs();
//       this.mkHmp321MafPlot();
//       this.mkHmp321MafPlot_useR();
 //      this.SiftGerp_Correlation();
      this.countDeleteriousHmp321();
//       this.checkTaxaNameSamewithFeiGroup();
       this.mkDepthOfHmp321();
       this.mkDepthSummary();
       this.countDeleteriousHmp32HighDepth();
//       this.mkBurdenFile_deprecated();
//       this.mkDistanceToB73();
//       this.mergeRecDelHighAndDsitanceToB73();
//       this.mkMDSplot();
//       this.countVarient();
//       
       
    }
    
    /**
     * 模板
     */
    
    public void demo(){
        String infileS = "";
        String outfileS = "";
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tGroup");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
    }
    
    public void countVarient(){
        String infileDirS = "/Volumes/LuLab3T_30/hmp321_agp4";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if(fs[i].isHidden()){
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        }
        fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf.gz");
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if(temp.startsWith("#"))continue;
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
    
    public void mkMDSplot(){
        //this.addGroupInfotoMds();
        this.mdsHmp321HighDepthFiltered();
        
    }
    
    public void mdsHmp321HighDepthFiltered(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/008_mds/003_mds_addGroup.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/008_mds/005_mds_addGroup_highDepth_Filtered.txt";
        String taxaFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/003_IBD/005_reccesiveDeleterious_merge_noTraitMixedUnknowteo.txt";
        RowTable<String> t = new RowTable<>(taxaFileS);
        List<String> taxaList = new ArrayList<>();
        taxaList = t.getColumn(0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa, key);
                if(index < 0) continue;
                bw.write(temp);
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
        
    }
    private void addGroupInfotoMds(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/008_mds/001_mds.txt";
        String taxaGroupFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/002_hmp321TaxaGroup/hmp321_taxaGroup.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/008_mds/003_mds_addGroup.txt";
        RowTable<String> t = new RowTable<>(taxaGroupFileS);
        HashMap<String, String> hm = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++){
            hm.put(t.getCellAsString(i, 0),t.getCellAsString(i, 1));
        }
        String temp = null;
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tGroup");bw.newLine();
            //String temp = null;
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                bw.write(temp);bw.write("\t");bw.write(hm.get(key));
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
        
        RowTable<String> tt = new RowTable<>(outfileS);
        Set<String> group = new HashSet<>();
        for (int i=0; i<tt.getRowNumber(); i++){
            group.add(tt.getCellAsString(i, 4));
        }
        System.out.println("1210份玉米种质的分类" + group);
        //[ss, Tripsacum, popcorn, sweet corn, mixed, sweet, Landrace, unknown, nss, ts, Teo, teo]  12种，太乱七八糟额！
       
    }
    
    private void mergeRecDelHighAndDsitanceToB73(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321_highDepth.txt";
        String distanceFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/003_IBD/002_1210taxaDistanceToB73.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/003_IBD/003_reccesiveDeleterious_merge.txt";
        RowTable<String> t = new RowTable<>(distanceFileS);
        HashMap<String,String> taxaDistance = new HashMap<>();
        for(int i= 0; i< t.getRowNumber(); i++){
            taxaDistance.put(t.getCellAsString(i, 1),t.getCellAsString(i, 2));
        }
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.write("\tIBS_Distance");bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String key = PStringUtils.fastSplit(temp).get(0);
                bw.write(temp);bw.write("\t");bw.write(taxaDistance.get(key));
                bw.newLine();
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /** 
     * Convert the matrix table to only oneVSone table
     */
    
    private void mkDistanceToB73(){
        String matrixFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/003_IBD/001B73matrix_sub.txt";
        String taxaListFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/002_hmp321TaxaGroup/TaxaList.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/003_IBD/002_1210taxaDistanceToB73.txt";
        RowTable<String> t = new RowTable<>(taxaListFileS);
        List<String> taxaList = new ArrayList<>();
        taxaList = t.getColumn(0);
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        try{
            BufferedReader br = IOUtils.getTextReader(matrixFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxon1\tTaxon2\tIBS_Distance");
            bw.newLine();
            String temp = br.readLine();
            List<String> matrixList = PStringUtils.fastSplit(temp);
            for(int i = 0; i < taxaList.size(); i++){
                bw.write("B73\t");
                bw.write(taxaList.get(i));bw.write("\t");
                int j = i+1;
                bw.write(matrixList.get(j));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * [ss, Tripsacum, popcorn, mixed, sweet, Landrace, unknown, ts, nss, Teo, teo] GroupSet 11个 （Teo 和teo重复） 多了 tripsacum(摩擦禾属)  mixed  unknow 
     */
    private void mkBurdenFile_deprecated(){
        
        String infileS = "/Users/Aoyue/Documents/additiveDeleterious_hmp321_highDepth.txt";
        String outfileS = "/Users/Aoyue/Documents/add_hmp321_burden.txt";
        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
        Set<String> group = new HashSet<>();
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            
            String header = br.readLine();
            bw.write(header);bw.write("\tGroupID");
            bw.newLine();
            String temp = null;
            List<String> l = null;
            while((temp=br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                String groupName = l.get(3);
                group.add( groupName);
                if(groupName.equals("unknow")) continue;
                if(groupName.equals("mixed")) continue;
                if(groupName.equals("unknow")) continue;
                
                if(groupName.equals("ss")){
                    StringBuilder sb = new StringBuilder();
                    sb.append(temp).append("\t").append("1");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                
            }
            bw.flush();bw.close();br.close();
            System.out.println(group);
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
            
        }
    }
    
    /**
     * 
     */
    private void countDeleteriousHmp32HighDepth () {
        String taxaSummaryFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/001_hmp321Depth/taxaDepth.summary.txt";
        String addInfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/additiveDeleterious_hmp321.txt";
        String addOutfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/additiveDeleterious_hmp321_highDepth.txt";
        String recInfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321.txt";
        String recOutfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321_highDepth.txt";
        String taxaGroupFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/004_taxaSummary/source/hmp321_taxaGroup_version3.txt";
        //1.54=3x, 2.75 = 5x
        double depthCut = 1.54;
        RowTable t = new RowTable (taxaSummaryFileS);
        ArrayList<String> taxaList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsDouble(i, 2) < depthCut) continue;
            taxaList.add(t.getCellAsString(i, 0));
        }
        String[] taxa = taxaList.toArray(new String[taxaList.size()]);
        Arrays.sort(taxa);
        t = new RowTable (taxaGroupFileS);
        HashMap<String, String> taxaGroupMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            taxaGroupMap.put(t.getCellAsString(i, 0), t.getCellAsString(i, 1));
        }
        t = new RowTable (addInfileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(addOutfileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth\tGroup\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 1))/Double.valueOf(t.getCellAsDouble(i, 2));
                bw.write(taxa[index]+"\t"+t.getCellAsString(i, 1)+"\t"+t.getCellAsString(i, 2)+"\t"+taxaGroupMap.get(taxa[index])+"\t"+String.valueOf(ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        t = new RowTable (recInfileS);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(recOutfileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth\tGroup\tRatio");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                int index = Arrays.binarySearch(taxa, t.getCellAsString(i, 0));
                if (index < 0) continue;
                double ratio = Double.valueOf(t.getCellAsDouble(i, 1))/Double.valueOf(t.getCellAsDouble(i, 2));
                bw.write(taxa[index]+"\t"+t.getCellAsString(i, 1)+"\t"+t.getCellAsString(i, 2)+"\t"+taxaGroupMap.get(taxa[index])+"\t"+String.valueOf(ratio));
                bw.newLine();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkDepthSummary () {
        //思想：对每一个taxa做统计，每读进一个taxa，就统计所有深度的和，再除以表格的行数，也即是抽查的位点数，最终得出总得深度除以位点数，得出平均每个位点的深度。
        String taxaDepthDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/001_hmp321Depth/taxa";
        String taxaSummaryFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/001_hmp321Depth/taxaDepth.summary.txt";
        File[] fs = new File (taxaDepthDirS).listFiles();
        Arrays.sort(fs);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(taxaSummaryFileS);
            bw.write("Taxa\tID\tMeanDepth");
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                RowTable t = new RowTable (fs[i].getAbsolutePath());
                String taxaName = String.valueOf(t.getHeader().get(0)).replaceFirst("_siteDepth", "");
                double value = 0;
                for (int j = 0; j < t.getRowNumber(); j++) {
                    value+=t.getCellAsDouble(j, 0);
                }
                bw.write(taxaName+"\t"+String.valueOf(i+1)+"\t"+String.valueOf((double)value/t.getRowNumber()));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void mkDepthOfHmp321 () {
        //思想：随机抽取1万个位点，进行每个taxa的每个位点的depth统计
        /**hmpInfoFileS file
         * Chr	Pos	Ref	Alt	Ancestral(Gerp)	Major	Minor	MajorAlleleFrequency	MinorAlleleFrequency	SiteDepth	HetCount
            10	228730	A	T	NA	A	T	0.986692	0.013307985	3591	1
            10	228743	T	<DEL>	NA	T	<DEL>	0.9964455	0.0035545023	4072	0
            10	228744	A	<DEL>	NA	A	<DEL>	0.9964413	0.0035587188	4054	0
         */
        String hmpFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10.vcf.gz";
        String hmpInfoFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter/hmp321Info_filter_chr010_AGPv4_AnnoDB.txt";
        String taxaDepthDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/001_hmp321Depth/taxa/";
        int snpNum = 0;
        //int size = 20000; //这里我增加了抽样的位点数
        int size = 50000;
        try {
            BufferedReader br = IOUtils.getTextReader(hmpInfoFileS);
            String temp = br.readLine();
            int cnt = 0; 
            while ((temp = br.readLine()) != null) {
                cnt++;
            }
            snpNum = cnt;
            /*建立一个10000个 indices的整型数组， 赋值为indice[0]= */
            int[] indices = new int[size]; 
            for (int i = 0; i < size; i++) {
                indices[i] = (int)(Math.random()*snpNum); //随机挑选10000个数目
            }
            Arrays.sort(indices);
            br = IOUtils.getTextGzipReader(hmpFileS);
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l = PStringUtils.fastSplit(temp, "\t");
            String[] taxa = new String[l.size()-9];
            for (int i = 0; i < taxa.length; i++) {
                taxa[i] = l.get(i+9);
            }
            TIntArrayList[] depthList = new TIntArrayList[taxa.length]; //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
            for (int i = 0; i < taxa.length; i++) depthList[i] = new TIntArrayList(); //对list进行初始化
            cnt = 0;
            while ((temp = br.readLine()) != null) {
                cnt++; // 对snp开始计数
                if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines");
                int idx = Arrays.binarySearch(indices, cnt-1);
                if (idx < 0) continue;
                l = PStringUtils.fastSplit(temp, "\t");
                for (int i = 0; i < taxa.length; i++) {
                    String genoS = l.get(i+9);
                    if (genoS.startsWith(".")) {
                        depthList[i].add(0);
                        continue;
                    }
                    List<String> ll = PStringUtils.fastSplit(genoS, ":");
                    List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                    int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1));
                    depthList[i].add(depth);
                }
            }
            for (int i = 0; i < taxa.length; i++) {
                String outfileS = new File (taxaDepthDirS, PStringUtils.getNDigitNumber(5, i+1)+"depth.txt").getAbsolutePath();
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    bw.write(taxa[i]+"_siteDepth");
                    bw.newLine();
                    int[] depth = depthList[i].toArray();
                    for (int j = 0; j < depth.length; j++) {
                        bw.write(String.valueOf(depth[j]));
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void checkTaxaNameSamewithFeiGroup(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/002_hmp321TaxaGroup/hmp321_taxaGroup.txt";
        String infileS2 = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/002_hmp321TaxaGroup/TaxaList.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/002_hmp321TaxaGroup/CompareTaxaList.txt";
        RowTable<String> t = new RowTable<>(infileS2);
        List<String> taxaListYue = t.getColumn(0);
        String[] taxaYue = taxaListYue.toArray(new String[taxaListYue.size()]);
        Arrays.sort(taxaYue);
        t = new RowTable<>(infileS);
        HashMap taxaVSGroup = new HashMap();
        List<String> taxaListFei = t.getColumn(0);
        String[] taxaFei = taxaListFei.toArray(new String[taxaListFei.size()]);
        Arrays.sort(taxaFei);
        int cnt = 0;
        for(int i = 0; i< t.getRowNumber(); i++){
            String testTaxa = t.getCell(i, 0);
            String group = t.getCell(i, 1);
            taxaVSGroup.put(testTaxa, group);
            int index = Arrays.binarySearch(taxaYue, testTaxa);
            if (index <0){
                System.out.println(testTaxa + "\t"+i +" is not the same with ME");
                System.out.println(t.getRow(i));
            }
            
        }
       
        
        /**
         * 08-64	1 is not the same with Fei
            282set_DE-2	1020 is not the same with Fei
            成功构建 (总时间: 0 秒)
            * 
            * 80-64	68 is not the same with ME
            [80-64, unknown]
            282set_DE	192 is not the same with ME
            [282set_DE, unknown]
            1210
         */
        
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa_Yue" + "\t" + "Taxa_Fei" + "\t" + "TaxaGroup");
            bw.newLine();
            for(int i = 0; i < taxaYue.length;i++){
                
                bw.write(taxaYue[i] + "\t" + taxaFei[i] + "\t" + String.valueOf(taxaVSGroup.get(taxaFei[i])));
                bw.newLine();
                if (taxaYue[i].equals(taxaFei[i])){
                    cnt++;
                }
                else{
                    System.out.println(taxaYue[i] + "\t" + taxaFei[i] + "\t" + String.valueOf(i+1) + "th row are not the same");
                    
                }
            }
            bw.flush();
            bw.close();
             System.out.println(String.valueOf(cnt));
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    /**
     * 大方向是 每个位点pos 的ref 和 derived allele是多变的，所以在确定显隐性的时候每个位点情况都是不同。
     * 对于10条染色体的每个位点，都是这样每个位点每个位点进行统计的，所有位点的统计结果都汇总到taxa中去，
     * 这里以taxa为单位，可以计算出 dele 的位点数，也可以算出，在这些为点数中，有多少个位点是 dele。 
     */
      private void countDeleteriousHmp321 () {
        String hmpDirS = "/Volumes/Lulab3T_14/hmp321_agp4";
        String infoFileS = "/Users/Aoyue/Documents/Data/referenceGenome/position/ChrLenCentPosi_agpV4.txt";
        String deleFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class/Non_Synonymous_Deleterious_High_GERP.txt";
        String addCountFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/additiveDeleterious_hmp321.txt";
        String recCountFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/reccesiveDeleterious_hmp321.txt";
        int minDepth = 2;//inclusive
        /**
         *将info读进表格
         * 定义一个整型数组 chrLength[] 分别为{1，2，3，4，5，6，7，8，9，10}
         * 定义一个整型集合chrList 分别为[1,2,3,4,5,6,7,8,9,10]
         * 
         * 将deleFileS读进表格
         * 定义一个整型集合类的数组 posList[10]，并进行初始化集合， 分别为 posList[0] posList[1]...
         * 定义一个字符型集合类的数组 charList[10]，并进行初始化集合， 分别为charList[0] charList[1]...
         * 如果DerivedAllele个数(A, T)大于1，跳过
         * 将表格pos信息读进对应的posList[i];将DerivedAllele信息读进对应的charList[i]
         * 将集合类的数组posList[10]和charList[10] 转化为二维数组 delePos[10][]  deleChar[10][]
         * 
         */
        RowTable t = new RowTable (infoFileS);
        int chrNum = t.getRowNumber();
        int[] chrLength = new int[chrNum];
        ArrayList<Integer> chrList = new ArrayList();
        for (int i = 0; i < chrNum; i++) {
            chrLength[i] = t.getCellAsInteger(i, 1);
            chrList.add(i+1);
        }
        /**
         * Chr	Pos	MinorAllele	MAF	DerivedAllele	DAF
            1	93755	T	0.0031847134	T	0.0031847134
            1	122211	T	0.0018773467	T	0.0018773467
         */
        t = new RowTable (deleFileS);
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        TCharArrayList[] charList = new TCharArrayList[chrNum];
        for (int i = 0; i < chrNum; i++) { //集合类数组，要初始化每一个集合
            posList[i] = new TIntArrayList();
            charList[i] = new TCharArrayList();
        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            int index = t.getCellAsInteger(i, 0)-1; 
            if (t.getCellAsString(i, 4).length()>1) continue; /*如果DerivedAllele个数大于1，跳过*/
            posList[index].add(t.getCellAsInteger(i, 1));
            charList[index].add(t.getCellAsString(i, 4).charAt(0));
        }
        int[][] delePos = new int[chrNum][];
        char[][] deleChar = new char[chrNum][];
        for (int i = 0; i < chrNum; i++) {
            delePos[i] = posList[i].toArray();
            deleChar[i] = charList[i].toArray();
        }
        /**
         * Get taxa's name 并建立数组taxa[]
         * 建立double型数组addCount[]， 整型数组recCount[] siteWithMinDepthCount[]
         * 对每条染色体做多线程处理，对每个含有有害突变的位点delePos进行 addCount（加性效应0/1）计数，recCount (隐形效应计数)计数， 深度（AD（0，1））大于2的taxa数目计数。
         * 
         */
        String hmpChr10FileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/hmp321_agpv4_chr10_10000.vcf";
        //hmpChr10FileS = new File(hmpDirS, hmpChr10FileS).getAbsolutePath();
        File hmpChr10Flie = new File(hmpChr10FileS);
        String[] taxa = null;
        try {
            BufferedReader br = IOUtils.getTextReader(hmpChr10FileS);
            String temp = br.readLine();
            while ((temp = br.readLine()).startsWith("##")) {} //直到循环体结束，temp不以 ## 开头
            List<String> l = PStringUtils.fastSplit(temp, "\t");
            taxa = new String[l.size()-9];
            for (int i = 9; i < l.size(); i++) {
                taxa[i-9] = l.get(i);
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        int taxaNum = taxa.length;
        double[] addCount = new double[taxa.length];
        int[] recCount = new int[taxa.length];
        int[] siteWithMinDepthCount = new int[taxa.length];  //这里是什么意思？？
        chrList.parallelStream().forEach(chr -> {
            String hmpFileS = "hmp321_agpv4_chr"+String.valueOf(chr)+".vcf.gz";
            hmpFileS = new File(hmpDirS, hmpFileS).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(hmpFileS);
            int chrIndex = chr-1;
            try {
                String temp = br.readLine();
                while ((temp = br.readLine()).startsWith("##")) {} //当temp = #CHROM时，循环体挑出，此时执行下面语句块。
                int cnt = 0;
                while ((temp = br.readLine()) != null) { //表示temp又读了一行，已经进入到10 16296阶段了。。。不要误解为还在#CHROM阶段。
                    cnt++;
                    if (cnt%1000000 == 0) System.out.println(String.valueOf(cnt)+" lines on chr "+String.valueOf(chr));
                    List<String> l = PStringUtils.fastSplit(temp.substring(0, 50), "\t");
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(delePos[chrIndex], pos);
                    if (index < 0) continue; //前面建立的 delePos[][] 和deleChar[][] 都是为现在在vcf文件中找位置贡献的，不是有害突变的位点，都过滤。
                    l = PStringUtils.fastSplit(temp, "\t");
                    int[] idx = new int[2]; //idx是什么？？
                    if (l.get(3).charAt(0) == deleChar[chrIndex][index]) { //如果ref allele = deleChar allele
                        idx[0] = 0; idx[1] = 1;
                    }
                    else idx[0] = 1; idx[1] = 0;
                    for (int i = 0; i < taxaNum; i++) {
                        String genoS = l.get(i+9); //指的是 GT:AD:GL 信息
                        if (genoS.startsWith(".")) continue; //如果以.开头，说明没有基因型信息，此位点没有测到。
                        List<String> ll = PStringUtils.fastSplit(genoS, ":"); //分开为GT   AD   GL三类
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ","); //lll指将AD提取出来，并以"，"号分割。如 0/0:1,0:0,3,28中 ，1，0分别代表ref和alt的测序深度
                        int depth = Integer.valueOf(lll.get(0))+Integer.valueOf(lll.get(1)); //总得测序深度等于 ref + alt
                        if (depth < minDepth) continue; //最小测序深度是2，如果小于2，则弃用
                        lll = PStringUtils.fastSplit(ll.get(0), "/"); //这里lll指的是基因型GT，lll被重新赋值，之前代表的是AD
                        int v1 = Integer.valueOf(lll.get(0)); //v1等于 ref
                        int v2 = Integer.valueOf(lll.get(1)); // v2 等于 alt
                        int sum = 0;
                        if (v1 == idx[0]) sum++; 
                        if (v2 == idx[0]) sum++;
                        if (sum == 0) {}
                        else if (sum == 1) {
                            addCount[i] += 0.5; //????这一段实在是看不懂啊！！！！！！！！1
                        }
                        else {
                            addCount[i] += 1;
                            recCount[i] += 1;
                        }
                        siteWithMinDepthCount[i]++; //查看过滤了位点后，每个taxa有多少个位点被测到了！！！！
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        }); 
        try {
            BufferedWriter bw = IOUtils.getTextWriter(addCountFileS);
            bw.write("Taxa\tDeleteriousCountPerHaplotype\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < addCount.length; i++) {
                bw.write(taxa[i]+"\t"+String.valueOf(addCount[i])+"\t"+String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
            bw = IOUtils.getTextWriter(recCountFileS);
            bw.write("Taxa\tDeleteriousCountPerLine\tSiteCountWithMinDepth");
            bw.newLine();
            for (int i = 0; i < recCount.length; i++) {
                bw.write(taxa[i]+"\t"+String.valueOf(recCount[i])+"\t"+String.valueOf(siteWithMinDepthCount[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private void SiftGerp_Correlation(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/002_hmp321SiftTrans_filter/hmp321Info_filterbySift_chr010.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/006_SiftGerpCorrelation/siftGerpCorrelation.pdf";
        RowTable<String> t = new RowTable (infileS);
        int[] index = PArrayUtils.getRandomIntArray(t.getRowNumber(), 20000); //返回整型类的数组，抽取1万个位点进行测试！
        TDoubleArrayList siftList = new TDoubleArrayList();
        TDoubleArrayList gerpList = new TDoubleArrayList();
        for (int i = 0; i < index.length; i++) {
            if (t.getCell(index[i], 17).startsWith("N")) continue; //过滤sift值是NA的位点
            if (!t.getCell(i, 16).startsWith("NON")) continue; //过滤同义突变的位点
            if (t.getCell(index[i], 15).startsWith("N")) continue; //过滤Gerp是NA的位点
            siftList.add(Double.valueOf(t.getCell(index[i], 17))); //将SIFT值加入
            gerpList.add(Double.valueOf(t.getCell(index[i], 15))); // 将GERP值加入
        }
        ScatterPlot s = new ScatterPlot (siftList.toArray(), gerpList.toArray());
        s.setColor(255, 0, 0, 20); //
        s.setPlottingCharacter(16);
        s.setTitle("Conservation Vs effect of variants in coding sequence");
        s.setXLim(0, 1);
        s.setYLim(-10, 5);
        s.setXLab("SIFT score");
        s.setYLab("GERP score");
        s.showGraph();
        s.saveGraph(outfileS);
    }
    
    private void mkHmp321MafPlot_useR(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter/hmp321Info_filter_chr001_AGPv4_AnnoDB.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/005_hmp321MafPlot/maf.pdf";
        BufferedReader br = IOUtils.getTextReader(infileS);
        try{
            String header = br.readLine();
            String temp = null;
            List<String> l = null;
            TDoubleArrayList mafList = new TDoubleArrayList();
            int cnt = 0;
            int random = 0;
            while((temp = br.readLine()) != null){
                double r = Math.random();
                if (r > 0.002) continue;
                random++;
                l = PStringUtils.fastSplit(temp);
                if((temp = l.get(7)).equals("NA")){
                    cnt++;
                }
                else{
                    mafList.add(Double.parseDouble(l.get(7)));
                }
            }
            System.out.println(cnt + " is NA");
            System.out.println(random + " is random");
            br.close();
            double[] mafValue = mafList.toArray();
            Histogram dd = new Histogram(mafValue);
            dd.setTitle("");
            dd.setXLab("Minor allele frequency");
            dd.setYLab("Frequency");
            dd.setBreakNumber(20);
            dd.setXLim(0, 0.5);
            dd.saveGraph(outfileS);
            dd.showGraph();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
  
    private void mkHmp321MafPlot () {
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter/hmp321Info_filter_chr001_AGPv4_AnnoDB.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/005_hmp321MafPlot/hmp321MafFre.txt";
        TDoubleArrayList list = new TDoubleArrayList();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            int sum = 0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(7).contains(",")) continue;
                list.add(Double.valueOf(l.get(7)));
                sum++;
                if (sum%100000 == 0) System.out.println(sum);
            }
            int size = 20; //分成20个bin,每个bin 0.025
            double[] bound = new double[size];
            int[] count = new int[size];
            double[] ratio = new double[size];
            for (int i = 1; i < bound.length; i++) { //bound[1]=0.025; bound[2]=0.05;bound[3]=0.075;bound[19]=0.475; 注意，bound[0]=0(double类型默认是0);
                bound[i] = (double)0.5/size*i;
            }
            double[] value = list.toArray();// 将集合转化为数组，遍历数组，在bound中搜索定位区间，并统计各区间的个数，后续将各区间的个数除以总得数组，计算出频率ratio.
            for (int i = 0; i < value.length; i++) {
                int index = Arrays.binarySearch(bound, value[i]);
                if (index < 0) index = -index-2;
                count[index]++; 
            }
            for (int i = 0; i < ratio.length; i++) {
                ratio[i] = (double)count[i]/sum;
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("MAF\tFrequency");
            bw.newLine();
            for (int i = 0; i < ratio.length; i++) {
                bw.write(String.valueOf(bound[i])+"\t"+String.valueOf(ratio[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
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
            mafList[i] = new TDoubleArrayList(); 
            dafList[i] = new TDoubleArrayList(); 
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
                dafFrequency[i][j] = dafFrequency[i][j]/dafCount[i]; //因为daf值不是每个都有，有些pos是NA值，所以需要重新计算。
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
            for (int i = 0; i < mafFrequency[0].length; i++) { //i小于第一个文件的长度100，
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
    
    
    public void binarySearchtest(){
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
     * 
     */
    private void testarraylength(){
        int[][] a = new int[3][4];
        //int[] bound = new int[5];
        int bound[] = {1,2,3,4,5};
        //int[][] a = {1,2,3,4,5,6,7,8,9,10,11,12};
        for(int i = 0; i<a[0].length; i++){
            System.out.println(bound[i]);
        }
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
                tranNameList.add(t.getCellAsString(i, 0));
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
                    if (t1.getCell(j, 16).equals("SYNONYMOUS")) {   //16	17	18 Variant_type	SIFT_score	Transcript
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


    /**
     * 目的：
     */
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
        
        //*********************************** START1 ***********************************//
        //该段代码的作用是，通过读取每个基因，得到最长转录本的名字，计算该转录本的长度。
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int chrIndex = gf.getChromosomeOfGene(i)-1;
            /*这个地方是先过滤数据，将定位在11号12号染色体上的基因过滤掉，并且跳出循环*/
            if (chrIndex >9) {
                cntchr11and12++;
                continue;
            }
            cntchr1to10++; //能够得到1-10号染色体的基因数目
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String geneName = gf.getTranscriptName(i, longTransIndex); //得到最长的转录本的名字
            //genes[i] = geneName;
            genesList.add(geneName);
            List<Range> cdsList = gf.getCDSList(i, longTransIndex); /*得到基因的最长转录本的CDSList*/
            int cnt = 0;
         
            
        /*对于每一个基因的编码序列，还有很多个cds片段，即cdsList；我们对cdsList进行for循环，得到每个cds的起始和终止位置，从而计算出总长*/    
            for (int j = 0; j < cdsList.size(); j++) {
                int rStart = cdsList.get(j).start;
                int rEnd = cdsList.get(j).end;
                for (int k = rStart; k < rEnd; k++) {
                    /*posGeneMap是一个HashMap数组，一条染色体对应一个String类型的ArrayList；
                    故得到该位点的所属基因名字列表，如果该位点不含基因名，就将 genename赋值给该位点，完善posGeneMap
                    否则，如果该位点含有其他基因的基因名字，依旧把genename赋值给该位点*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(k); //建立map的关系，那个位点对应哪个list HashMap<Integer, ArrayList<String>>[] posGeneMap = new HashMap[chrNum]; 
                    if (geneNameList == null) {
                        geneNameList = new ArrayList();
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList);
                    }
                    else {
                        geneNameList.add(geneName);
                        posGeneMap[chrIndex].put(k, geneNameList); /*最终将posGeneMap绘图完成*/
                    }
                    cnt++; /*每一个CDS位点相加，最终得到这个cds的长度。*/
                }
                
                // 最终cnt是一个基因的所有cdslist中，每个cds的每个位点包含的基因数目的总和
            } //该循环是一个基因的所有cds循环
            geneCDSLengthMap.put(geneName, cnt); //
        }
        //*********************************** END1 ***********************************//
        
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
                    if (l.get(3).contains("<") || l.get(3).contains(",")) continue; //过滤含有2个alt和含有indel的位点
                    int pos = Integer.valueOf(l.get(1));
                    /*一个位置对应一个genelist， 该位置可能是很多基因的位点，存入geneNameList;开始循环，如果搜索到这个基因，index为第几个基因，就让该基因的计数加一*/
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos);
                    if (geneNameList == null) continue; //说明该变异位点不在基因区
                    for (int i = 0; i < geneNameList.size(); i++) {
                        int index = Arrays.binarySearch(genes, geneNameList.get(i));//在基因库的第i个位置，该基因数加一
                        snpCount[index]++; //第i个位置的基因含有的snp数目；
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
            snpPos[chrIndex] = snpPosList.toArray(); //第n条染色体的在基因区间的snp所对应的pos的集合；
            snps[chrIndex] = snpList.toArray(); //第n条染色体的在区间内的snp所对应的alt的集合
            snpAnc[chrIndex] = snpAncList.toArray(); //第n条染色体的在区间内的snp所对应的ancestral allele的集合
        });
        
        /*分为4类： 同义突变 非同义突变 有害突变 高GERP值的有害突变 NA值的计数*/
        int[] synCount = new int[genes.length];
        int[] nonCount = new int[genes.length]; //非同义突变的数目
        int[] delCount = new int[genes.length];
        int[] delHGCount = new int[genes.length];
        int[] naCount = new int[genes.length];  //非同义突变，但是sift值是NA的数目
        
        int[] b73SynCount = new int[genes.length];
        int[] b73NonCount = new int[genes.length];
        int[] b73DelCount = new int[genes.length];
        int[] b73DelHGCount = new int[genes.length];
        int[] noAncCount = new int[genes.length];
        
        TIntArrayList[] delPosList = new TIntArrayList[chrNum]; //有害变异的位点集合
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
                    if (l.get(16).startsWith("NA")) continue; //Variant_type	SIFT_score	Transcript  16 17 18
                    if (l.get(18).startsWith("NA")) continue; //没有变异类型和转录本的位点，都过滤掉。
                    String gene = l.get(18);
                    int geneIndex = Arrays.binarySearch(genes, gene); //在genes数组里搜索 sift中的基因
                    if (geneIndex < 0) continue;
                    int pos = Integer.valueOf(l.get(1));
                    int index = Arrays.binarySearch(snpPos[chrIndex], pos); // 
                    if (index < 0) continue;
                    if (snps[chrIndex][index] != l.get(3).getBytes()[0]) continue; //再次验证，如果该位点的alt和数据库中的alt不一致，则过滤掉。
                    byte ref = l.get(2).getBytes()[0];
                    byte derivedState = -1; //mean ancestral allele is not defined
                    if (snpAnc[chrIndex][index] == snps[chrIndex][index]) { //如果ancestral allele 等于alt的话，derived allele就等于1
                        derivedState = 1; //mean b73 carries derived allele
                    }
                    else if (snpAnc[chrIndex][index] == ref) { //如果ancestral allele 等于ref的话，derived allele就等于0
                        derivedState = 0;
                    }
                    
                    if (derivedState == -1) noAncCount[geneIndex]++; //如果ancestral allele 不存在的话，derived allele就等于-1
                    
                    String type = null;
                    if (l.get(16).equals("NA")) {
                        
                    }
                    else {
                        if (l.get(16).equals("SYNONYMOUS")) { //如果type等于syn，那么 该位点所属的基因的syn属性就加一
                            type = "Syn";
                            synCount[geneIndex]++;
                            if (derivedState == 1) {  //如果derived allele就等于1，说明ancestral allele 等于alt， derived allele 等于ref； 如何计算daf,判断da是major还是minor，如果da是major那么daf=1-maf，如何判断major allele和minor allele？ 如果ref allele frequency > alt allele frequence,那么major是ref; 反之亦然；
                                b73SynCount[geneIndex]++; //如果参考基因组是 derived allele,那么就加一，为什么？
                            }
                        }
                        else {
                            type = "Non";
                            nonCount[geneIndex]++;
                            if (derivedState == 1) {
                                b73NonCount[geneIndex]++;
                            }
                            if (l.get(17).startsWith("NA")) { //index 17列是sift的值，NON-SYNONYMOUS 存在的情况下，sift可能有也可能没有。
                                    naCount[geneIndex]++; //是nonsynonymous类型但是没有sift值的个数
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
        
        /**
         * 
         */
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
                    ArrayList<String> geneNameList = posGeneMap[chrIndex].get(pos); //根据pos信息，得到该pos对应的gene name的集合
                    if (geneNameList == null) continue;
                    for (int i = 0; i < geneNameList.size(); i++) {
                        String gene = geneNameList.get(i);
                        if (gene == null) continue;
                        int geneIndex = Arrays.binarySearch(genes, gene);
                        
                        double treeValue = Double.valueOf(l.get(14));
                        double scoreValue = Double.valueOf(l.get(15));
                        if (treeValue == 0) continue; //过滤枝长是0的数目
                        gerpAlignCount[geneIndex]++; //如果枝长不是0，说明该位点存在保守不保守
                        gerpTree[geneIndex]+=treeValue;
                        gerpScore[geneIndex]+=scoreValue; //第i个基因的gerpscore的总和是多少
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
                if (naCount[i] == nonCount[i]) ifSiftAligned = 0; //非同义突变，但是sift值是NA的数目 等于 非同义突变的数目，那么说明在这个基因内部没有sift突变
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
   
    private void summarizeTranscript_deprecated(){
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
            int chrIndex = gf.getChromosomeOfGene(i)-1;
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
                if(gf.getChromosomeOfGene(i) > 10) continue;
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
    
    /**
     * 主要是想看，由最开始的数据库，只进行sift过滤，过滤掉了多少位点，保留了多少的有sift值的位点。因为技术过滤（深度，杂合度过滤）时，把某些有sift值的位点也过滤掉了，我们拿到的最终的位点不是最开始计算的所有sift值位点。
     */
    private void filterHmp321Info_bysiftTrans_useDataBase(){
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/annoDB_new";
        String outDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/002_hmp321SiftTrans_filter_useDatabase";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "hmp");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            String outfileS = new File(outDirS, "hmp321Info_filterbySift_" + f.getName().split("Info_")[1].replaceFirst("_AGPv4_AnnoDB.txt.gz", ".txt")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
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
                    if(trans.equals("NA"))continue; /*如果siftascore的值为NA，则无法判断其为有害或是中性突变。我们要筛选即有sift变异类型又有sift值的sites*/
                    cnt++;
                    bw.write(temp);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt)+"\t"+ f.getName().split("Info_")[1].replaceFirst("_AGPv4_AnnoDB.txt", "") + " trans sites");
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
                
            }
        });
        
    }
    
    /*
    Conservation Vs effect of variants in coding sequence
    */
    
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
                    if(trans.equals("NA"))continue; /*如果siftascore的值为NA，则无法判断其为有害或是中性突变。我们要筛选即有sift变异类型又有sift值的sites*/
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
        String outDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter";
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
                        double indiDepth = Double.valueOf(l.get(11))/Double.valueOf(l.get(12)); // 第11列是该位点的测序深度。12列是该位点有几个taxa被测到。 //求得是平均测序深度
                        if (indiDepth>maxIndiDepth) continue;
                        bw.write(temp);
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(f.getName().split("_")[1] + ": " + cnt + "sites are NA about basic annotation.");
                System.out.println(f.getName().split("_")[1] + ": " + cntfilter_number + "sites are filter under this standard."); //统计出去NA后，将要面临技术过滤的snp数
                
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
            for(int i = 0; i< depth.length;i++){
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
