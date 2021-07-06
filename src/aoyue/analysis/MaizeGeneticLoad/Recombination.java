 /*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import pgl.infra.anno.gene.GeneFeature;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
//import deprecated.analysis.maizeGeneticLoad.CrossMapUtils;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

 /**
 *
 * @author Aoyue
 */
 class Recombination {

     Recombination() {
//        this.convertHotspotsStart_AGPV2toAGPV4_deprecated();
//        this.convertHotspotsStart_AGPV2toAGPV4_cn_deprecated();
//        this.convertHotspotsEnd_AGPV2toAGPV4_deprecated();
//        this.convertHotspotsEnd_AGPV2toAGPV4_cn_deprecated();
//        this.checkNumberofDifChr_deprecated();
//        this.checkNumberofDifChr_cn_deprecated();
//        this.makeRecombinationPointTable_deprecated();
//        this.mkBinTable_deprecated();
//        
//        
//        this.calculateMidValue();
//        this.converMidValue_V2V4();
//        this.FilterMissMidValue();
//        this.mkRecombinationPointTable();
        this.mkBinTable_useMidValue();
//        this.testgetSubsetsIndicesBySubsetSize_Method();
        
         
    }
     
     /**
      * bound的长度是7,每个长度是100
      */
     private void testgetSubsetsIndicesBySubsetSize_Method(){
         int binsize = 100;
         int a = 699;
         int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(a, binsize);
         for (int i = 0; i< bound.length;i++){
             System.out.println(bound[i][0]);
             
         }
     }
     private void mkBinTable_useMidValue () {
        String transcriptFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String crossOverFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/recombinationPoint.txt";
        String deleteriousSNPFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class/Non_Synonymous_Deleterious_High_GERP.txt";
        String synSNPFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class/Synonymous.txt";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String infoFileS = "/Users/Aoyue/Documents/Data/referenceGenome/position/ChrLenCentPosi_agpV4.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/binTable.txt";
        int chrNum = 10;
        int binSize = 10000000; // 10个million
        /**
         * 每条染色体都有一个bounds[]；都有一个crossover[] 和 cdsLength[];都有一个del[]; 都有一个syn[];
        将中心粒的文件读进表格，Chromosome	Length(V4)	CentromereS	CentromereE
               1	307041717	136770000	137120000
        进行Bin分区；
        */
        int[][] bounds = new int[chrNum][];
        int[][] crossover = new int[chrNum][];
        int[][] cdsLength = new int[chrNum][];
        double[][] del = new double[chrNum][];
        double[][] syn = new double[chrNum][];
        RowTable<String> t = new RowTable (infoFileS);
        for (int i = 0; i < chrNum; i++) {
            int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(Integer.valueOf(t.getCell(i, 1)), binSize);
            bounds[i] = new int[bound.length]; //此步骤的目的是确定二维数组的第二维长度
            cdsLength[i] = new int[bound.length];
            crossover[i] = new int[bound.length];
            del[i] = new double[bound.length];
            syn[i] = new double[bound.length];
            for (int j = 0; j < bound.length; j++) {
                bounds[i][j] = bound[j][0]; //等于bound[][0] ，第一个元素，即bound的左边边界
            }
        }
        t = new RowTable (transcriptFileS);
        ArrayList<String> tranList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i, 0).equals("0")) continue;
            tranList.add(t.getCell(i, 0));
        }
        String[] transName = tranList.toArray(new String[tranList.size()]);
        Arrays.sort(transName);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String name = gf.getTranscriptName(i, longTransIndex); /*得到最长的转录本的名字*/
            if (Arrays.binarySearch(transName, name) < 0) continue; /*过滤genefeature里那些假基因*/
            int chrIndex = gf.getChromosomeOfGene(i) - 1; /*获得该基因的染色体位置*/
            List<Range> cds = gf.getCDSList(i, longTransIndex); /*得到基因的最长转录本的CDSList*/
            for (int j = 0; j < cds.size(); j++) {
                int index = Arrays.binarySearch(bounds[chrIndex], cds.get(j).start); //在bounds里搜索 该基因的第j个cds的起始位点
                if (index < 0) index = -index-2;
                if (index == bounds[chrIndex].length-1) { //如果刚好位于bounds最后一个bin起始位点，则cds长度等于起始位点加上 cds起始终止的差值
                    cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                }
                else {
                    if (cds.get(j).end <= bounds[chrIndex][index+1]) { //如果
                        cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                    }
                    else {
                        cdsLength[chrIndex][index]+=(bounds[chrIndex][index+1]-cds.get(j).start);
                        cdsLength[chrIndex][index+1]+=(cds.get(j).end-bounds[chrIndex][index+1]);
                    }
                }
            }
        }
        t = new RowTable (crossOverFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            crossover[chrIndex][index]++;
        }
        t = new RowTable (deleteriousSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            del[chrIndex][index]++;
        }
        t = new RowTable (synSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            syn[chrIndex][index]++;
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tBinStart\tCrossoverCount\tDelCount\tSynCount\tDelCountPerSite\tSynCountPerSite\tDelSynRatio");
            bw.newLine();
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < bounds[i].length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(i+1).append("\t").append(bounds[i][j]).append("\t").append(crossover[i][j]).append("\t").append(del[i][j]).append("\t").append(syn[i][j]);
                    sb.append("\t").append((double)del[i][j]/cdsLength[i][j]).append("\t").append((double)syn[i][j]/cdsLength[i][j]);
                    if (syn[i][j] == 0) { 
                        sb.append("\t").append(Double.NaN);
                    }
                    else {
                        sb.append("\t").append((double)del[i][j]/syn[i][j]);
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
     
     private void mkRecombinationPointTable(){
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/003minValueofCrossOvers_V4_filtered";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/recombinationPoint.txt";
        int chrNum = 10;
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "Ro");
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        for (int i = 0; i < posList.length; i++) { //每一个新的list都要new一下！！！因为上文只说明建立一个数组posList,大小为10.没有针对一个posList进行初始化。
            posList[i] = new TIntArrayList();
        }
        try {
            for (int i = 0; i < fs.length; i++) {
                BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).startsWith("het")) continue;
                    int chrIndex = Integer.valueOf(l.get(0))-1;
                    int pos = (Integer.valueOf(l.get(4)));
                    posList[chrIndex].add(pos);
                }
                br.close();
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < posList.length; i++) {
                int[] pos = posList[i].toArray();
                Arrays.sort(pos);
                for (int j = 0; j < pos.length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(pos[j]));
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
     
     /**
      * Filter the value that has the -1 value on chr and the value that convert to other chr, eg: 1   2345 ->  3   7629
      */
    public void FilterMissMidValue(){
        String infileCNS1 = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/minValueofCrossOvers_V4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_midValue_agpv4.txt";
        String infileUSS2 = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/minValueofCrossOvers_V4/RodgersMelnick2015PNAS_usnamImputedXOsegments_midValue_agpv4.txt";
        String outfileCNS1 = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/003minValueofCrossOvers_V4_filtered/RodgersMelnick2015PNAS_cnnamImputedXOsegments_midValue_agpv4_filtered.txt";
        String outfileUSS2 = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/003minValueofCrossOvers_V4_filtered/RodgersMelnick2015PNAS_usnamImputedXOsegments_midValue_agpv4_filtered.txt";
        int sum = 0;
        try{
            BufferedReader br = IOUtils.getTextReader(infileCNS1);
            BufferedWriter bw = IOUtils.getTextWriter(outfileCNS1);
            String header = br.readLine();
            bw.write(header); bw.newLine();
            String temp = null;
            List<String> l = null;
            int cnt=1;
            int cntunmap = 0;
            while((temp = br.readLine()) != null){
                cnt++;
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(0));
                /**
                 * A=10 B=20
                 * >=	检查左操作数的值是否大于或等于右操作数的值，如果是那么条件为真。	（A> = B）为假。
                   <=	检查左操作数的值是否小于或等于右操作数的值，如果是那么条件为真。	（A <= B）为真。
                 */
                if((1 <  cnt) && (cnt <  9181)){
                    if(chr == 1){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((9180 <  cnt) && (cnt <  15844)){
                    if(chr == 2){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((15843 <  cnt) && (cnt <  22677)){
                    if(chr == 3){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((22676 <  cnt) && (cnt <  28753)){
                    if(chr == 4){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((28752 <  cnt) && (cnt <  35497)){
                    if(chr == 5){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((35496 <  cnt) && (cnt <  40765)){
                    if(chr == 6){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((40764 <  cnt) && (cnt <  46643)){
                    if(chr == 7){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((46642 <  cnt) && (cnt <  52056)){
                    if(chr == 8){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((52055 <  cnt) && (cnt <  57196)){
                    if(chr == 9){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((57195 <  cnt) && (cnt <  61596)){
                    if(chr == 10){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(String.valueOf(cnt) + "  lines in cn.txt");
            int map = cnt - cntunmap;
            System.out.println(String.valueOf(map) + "   in cn is mappped");
            System.out.println(String.valueOf(cntunmap) + "  in cn is unmap or map to other chr ");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
            
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileUSS2);
            BufferedWriter bw = IOUtils.getTextWriter(outfileUSS2);
            String header = br.readLine();
            bw.write(header); bw.newLine();
            String temp = null;
            List<String> l = null;
            int cnt=1;
            int cntunmap = 0;
            while((temp = br.readLine()) != null){
                cnt++;
                l = PStringUtils.fastSplit(temp);
                int chr = Integer.parseInt(l.get(0));
                /**
                 * A=10 B=20
                 * >=	检查左操作数的值是否大于或等于右操作数的值，如果是那么条件为真。	（A> = B）为假。
                   <=	检查左操作数的值是否小于或等于右操作数的值，如果是那么条件为真。	（A <= B）为真。
                 */
                if((1 <  cnt) && (cnt <  28337)){
                    if(chr == 1){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((28336 <  cnt) && (cnt <  50364)){
                    if(chr == 2){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((50363 <  cnt) && (cnt <  71599)){
                    if(chr == 3){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((71598 <  cnt) && (cnt <  90913)){
                    if(chr == 4){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((90912 <  cnt) && (cnt <  112121)){
                    if(chr == 5){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((112120 <  cnt) && (cnt <  128755)){
                    if(chr == 6){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((128754 <  cnt) && (cnt <  147304)){
                    if(chr == 7){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((147303 <  cnt) && (cnt <  164970)){
                    if(chr == 8){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((164969 <  cnt) && (cnt <  181806)){
                    if(chr == 9){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
                if((181805 <  cnt) && (cnt <  195737)){
                    if(chr == 10){
                        bw.write(temp);
                        bw.newLine();
                    }
                    else{
                        cntunmap++;
                    }
                }
            }
            bw.flush();
            bw.close();
            br.close();
            System.out.println(String.valueOf(cnt) + "  lines in us.txt");
            int map = cnt - cntunmap;
            System.out.println(String.valueOf(map) + "  in us is mappped");
            
            System.out.println(String.valueOf(cntunmap) + "  in us is unmap or map to other chr ");
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
            
        }
        
        /**
         * 61595  lines in cn.txt
            54709   in cn is mappped
            6886  in cn is unmap or map to other chr 
            195736  lines in us.txt
            174619  in us is mappped
            21117  in us is unmap or map to other chr 
            成功构建 (总时间: 0 秒)
         */
    }
     
     
    private void converMidValue_V2V4(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/midValueofCrossOvers";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/minValueofCrossOvers_V4";
        String tempInputbedFileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/temp";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "Ro");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
        String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", "_agpv4.txt")).getAbsolutePath();
        String tempInputbedFileS = new File(tempInputbedFileDirS, f.getName().replaceFirst(".txt", "_agpv4.bed")).getAbsolutePath();
        String outputbedFileS = new File(tempInputbedFileDirS, f.getName().replaceFirst(".txt", "_agpv4.bed.con")).getAbsolutePath();
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz";
        RowTable<String> t = new RowTable (f.getAbsolutePath());
        int[] chr = new int[t.getRowNumber()]; 
        int[] pos = new int[t.getRowNumber()]; 
        for (int i = 0; i < chr.length; i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 4));
        }
//        CrossMapUtils cm = new CrossMapUtils(chr, pos, tempInputbedFileS);
//        cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath, tempInputbedFileS, outputbedFileS );
//        cm.convert();
//        List<int[]> l = cm.getConvertedCoordinate();
//        chr = l.get(0);
//        pos = l.get(1);
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.setCell(i, 0, String.valueOf(chr[i]));
            t.setCell(i, 4, String.valueOf(pos[i]));
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
        });
        /**
         * /Users/Aoyue/miniconda3/bin/python /Users/Aoyue/miniconda3/bin/CrossMap.py bed /Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_midValue_agpv4.bed /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_midValue_agpv4.bed.con
/Users/Aoyue/miniconda3/bin/python /Users/Aoyue/miniconda3/bin/CrossMap.py bed /Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_midValue_agpv4.bed /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_midValue_agpv4.bed.con
Table is written to /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/minValueofCrossOvers_V4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_midValue_agpv4.txt
Table is written to /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/minValueofCrossOvers_V4/RodgersMelnick2015PNAS_usnamImputedXOsegments_midValue_agpv4.txt
成功构建 (总时间: 13 秒)

         */
    }
        
    private void calculateMidValue(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv2";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/whywhywhy/midValueofCrossOvers";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "txt");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS,f.getName().replaceFirst(".txt", "_midValue.txt")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            try{
                String header = br.readLine();
                List<String> l = null;
                l = PStringUtils.fastSplit(header);
                bw.write(l.get(0));bw.write("\t");bw.write(l.get(1));bw.write("\t");bw.write(l.get(2));bw.write("\t");bw.write(l.get(3));
                bw.write("\t");bw.write("MidValue");bw.newLine();
                String temp = null;
                int cnt = 1;
                while((temp = br.readLine()) != null){
                    l = PStringUtils.fastSplit(temp);
                    StringBuilder sb = new StringBuilder();
                    int midValue = (Integer.valueOf(l.get(4))+Integer.valueOf(l.get(5)))/2;
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(2)).append("\t").append(l.get(3))
                            .append("\t").append(String.valueOf(midValue));
                    bw.write(sb.toString());
                    bw.newLine();
                    cnt++;
                }
                bw.flush();
                bw.close();
                br.close();
                System.out.println(String.valueOf(cnt)+"\t" + "sites    " + f.getName());
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        });
        
       /**
        * 61595	sites    RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt
            195736	sites    RodgersMelnick2015PNAS_usnamImputedXOsegments.txt
            成功构建 (总时间: 0 秒)
        */
     }
    
     private void mkBinTable_deprecated () {
        String transcriptFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String crossOverFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/recombinationPoint.txt";
        String deleteriousSNPFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class/Non_Synonymous_Deleterious_High_GERP.txt";
        String synSNPFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/class/Synonymous.txt";
        String geneFeatureFileS = "/Users/Aoyue/Documents/Data/referenceGenome/GeneAnnotation/Zea_mays.AGPv4.38.pgf";
        String infoFileS = "/Users/Aoyue/Documents/Data/referenceGenome/position/ChrLenCentPosi_agpV4.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/binTable.txt";
        int chrNum = 10;
        int binSize = 10000000; // 10个million
        /**
         * 每条染色体都有一个bounds[]；都有一个crossover[] 和 cdsLength[];都有一个del[]; 都有一个syn[];
        将中心粒的文件读进表格，Chromosome	Length(V4)	CentromereS	CentromereE
               1	307041717	136770000	137120000
        进行Bin分区；
        */
        int[][] bounds = new int[chrNum][];
        int[][] crossover = new int[chrNum][];
        int[][] cdsLength = new int[chrNum][];
        double[][] del = new double[chrNum][];
        double[][] syn = new double[chrNum][];
        RowTable<String> t = new RowTable (infoFileS);
        for (int i = 0; i < chrNum; i++) {
            int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(Integer.valueOf(t.getCell(i, 1)), binSize);
            bounds[i] = new int[bound.length];
            cdsLength[i] = new int[bound.length];
            crossover[i] = new int[bound.length];
            del[i] = new double[bound.length];
            syn[i] = new double[bound.length];
            for (int j = 0; j < bound.length; j++) {
                bounds[i][j] = bound[j][0]; //等于bound[][0] ，第一个元素，即bound的左边边界
            }
        }
        t = new RowTable (transcriptFileS);
        ArrayList<String> tranList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCell(i, 0).equals("0")) continue;
            tranList.add(t.getCell(i, 0));
        }
        String[] transName = tranList.toArray(new String[tranList.size()]);
        Arrays.sort(transName);
        GeneFeature gf = new GeneFeature(geneFeatureFileS);
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            int longTransIndex = gf.getLongestTranscriptIndex(i);
            String name = gf.getTranscriptName(i, longTransIndex); /*得到最长的转录本的名字*/
            if (Arrays.binarySearch(transName, name) < 0) continue; /*过滤genefeature里那些假基因*/
            int chrIndex = gf.getChromosomeOfGene(i) - 1;
            List<Range> cds = gf.getCDSList(i, longTransIndex); /*得到基因的最长转录本的CDSList*/
            for (int j = 0; j < cds.size(); j++) {
                int index = Arrays.binarySearch(bounds[chrIndex], cds.get(j).start); //在bounds里搜索 该基因的第j个cds的起始位点
                if (index < 0) index = -index-2;
                if (index == bounds[chrIndex].length-1) { //如果刚好位于bounds最后一个bin起始位点，则cds长度等于起始位点加上 cds起始终止的差值
                    cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                }
                else {
                    if (cds.get(j).end <= bounds[chrIndex][index+1]) { //如果
                        cdsLength[chrIndex][index]+=(cds.get(j).end-cds.get(j).start);
                    }
                    else {
                        cdsLength[chrIndex][index]+=(bounds[chrIndex][index+1]-cds.get(j).start);
                        cdsLength[chrIndex][index+1]+=(cds.get(j).end-bounds[chrIndex][index+1]);
                    }
                }
            }
        }
        t = new RowTable (crossOverFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            crossover[chrIndex][index]++;
        }
        t = new RowTable (deleteriousSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            del[chrIndex][index]++;
        }
        t = new RowTable (synSNPFileS);
        for (int i = 0; i < t.getRowNumber(); i++) {
            int chrIndex = Integer.valueOf(t.getCell(i, 0))-1;
            int index = Arrays.binarySearch(bounds[chrIndex], Integer.valueOf(t.getCell(i, 1)));
            if (index < 0) index = -index-2;
            syn[chrIndex][index]++;
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tBinStart\tCrossoverCount\tDelCount\tSynCount\tDelCountPerSite\tSynCountPerSite\tDelSynRatio");
            bw.newLine();
            for (int i = 0; i < chrNum; i++) {
                for (int j = 0; j < bounds[i].length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(i+1).append("\t").append(bounds[i][j]).append("\t").append(crossover[i][j]).append("\t").append(del[i][j]).append("\t").append(syn[i][j]);
                    sb.append("\t").append((double)del[i][j]/cdsLength[i][j]).append("\t").append((double)syn[i][j]/cdsLength[i][j]);
                    if (syn[i][j] == 0) { 
                        sb.append("\t").append(Double.NaN);
                    }
                    else {
                        sb.append("\t").append((double)del[i][j]/syn[i][j]);
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
     
     private void makeRecombinationPointTable_deprecated () {
        String inDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_final";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/recombinationPoint.txt";
        int chrNum = 10;
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "txt");
        TIntArrayList[] posList = new TIntArrayList[chrNum];
        for (int i = 0; i < posList.length; i++) { //每一个新的list都要new一下！！！因为上文只说明建立一个数组posList,大小为10.没有针对一个posList进行初始化。
            posList[i] = new TIntArrayList();
        }
        try {
            for (int i = 0; i < fs.length; i++) {
                BufferedReader br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = br.readLine();
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).startsWith("het")) continue;
                    int chrIndex = Integer.valueOf(l.get(0))-1;
                    int pos = (Integer.valueOf(l.get(4))+Integer.valueOf(l.get(5)))/2;
                    posList[chrIndex].add(pos);
                }
                br.close();
            }
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Chr\tPos");
            bw.newLine();
            for (int i = 0; i < posList.length; i++) {
                int[] pos = posList[i].toArray();
                Arrays.sort(pos);
                for (int j = 0; j < pos.length; j++) {
                    bw.write(String.valueOf(i+1)+"\t"+String.valueOf(pos[j]));
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
     
     public void checkNumberofDifChr_cn_deprecated(){
         String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4.txt";
         String outfileS ="/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/cnnamImputedXOsegments_agpv4.txt";
         
         int cntNNNN =0;
         int cntYYNN =0;
         int cntNNYY =0;
         int cntYYYY =0;
         int cntsame =0;
         int cntdiff =0;
         RowTable<String> t = new RowTable<>(infileS);
         boolean[] ifOut = new boolean[t.getRowNumber()];
         for (int i = 0; i< t.getRowNumber(); i++){
             String chr1st = t.getCell(i, 0);
             String chr2nd = t.getCell(i, 6);
             String start = t.getCell(i, 4);
             String end = t.getCell(i, 5);
             if (chr1st.equals("-1") && chr2nd.equals("-1")) cntNNNN++;
             if ((!chr1st.equals("-1")) && chr2nd.equals("-1")) cntYYNN++;
             if (chr1st.equals("-1") && (!chr2nd.equals("-1"))) cntNNYY++;
             if ((!chr1st.equals("-1")) && (!chr2nd.equals("-1"))) {
                 cntYYYY++;
                 if (chr1st.equals(chr2nd)){
                     cntsame++;
                    //l.add(t.getRow(i).toArray(rowContent));
                    ifOut[i] = true;
                 }
                 else{
                     cntdiff++;
                 }
             }
         }
         t.removeColumn(6);
         t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
         System.out.println("The number of sites is: " + t.getRowNumber());
         System.out.println("The number of cntNNNN sites is: " + cntNNNN);
         System.out.println("The number of cntYYNN sites is: " + cntYYNN);
         System.out.println("The number of cntNNYY sites is: " + cntNNYY);
         System.out.println("The number of cntYYYY sites is: " + cntYYYY);
         System.out.println("The number of cntsame sites is: " + cntsame);
         System.out.println("The number of cntdiff sites is: " + cntdiff);
         
         /**
          * Table is written to /Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/cnnamImputedXOsegments_agpv4.txt
            The number of sites is: 61594
            The number of cntNNNN sites is: 540
            The number of cntYYNN sites is: 1015
            The number of cntNNYY sites is: 0
            The number of cntYYYY sites is: 60039
            The number of cntsame sites is: 57871
            The number of cntdiff sites is: 2168
          */
     }
     
     public void checkNumberofDifChr_deprecated(){
         String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4.txt";
         String outfileS ="/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/usnamImputedXOsegments_agpv4.txt";
         
         int cntNNNN =0;
         int cntYYNN =0;
         int cntNNYY =0;
         int cntYYYY =0;
         int cntsame =0;
         int cntdiff =0;
         RowTable<String> t = new RowTable<>(infileS);
         boolean[] ifOut = new boolean[t.getRowNumber()];
         for (int i = 0; i< t.getRowNumber(); i++){
             String chr1st = t.getCell(i, 0);
             String chr2nd = t.getCell(i, 6);
             String start = t.getCell(i, 4);
             String end = t.getCell(i, 5);
             if (chr1st.equals("-1") && chr2nd.equals("-1")) cntNNNN++;
             if ((!chr1st.equals("-1")) && chr2nd.equals("-1")) cntYYNN++;
             if (chr1st.equals("-1") && (!chr2nd.equals("-1"))) cntNNYY++;
             if ((!chr1st.equals("-1")) && (!chr2nd.equals("-1"))) {
                 cntYYYY++;
                 if (chr1st.equals(chr2nd)){
                     cntsame++;
                    //l.add(t.getRow(i).toArray(rowContent));
                    ifOut[i] = true;
                 }
                 else{
                     cntdiff++;
                 }
             }
         }
         t.removeColumn(6);
         t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
         System.out.println("The number of sites is: " + t.getRowNumber());
         System.out.println("The number of cntNNNN sites is: " + cntNNNN);
         System.out.println("The number of cntYYNN sites is: " + cntYYNN);
         System.out.println("The number of cntNNYY sites is: " + cntNNYY);
         System.out.println("The number of cntYYYY sites is: " + cntYYYY);
         System.out.println("The number of cntsame sites is: " + cntsame);
         System.out.println("The number of cntdiff sites is: " + cntdiff);
         /*
        The number of sites is: 195735
        The number of cntNNNN sites is: 1716
        The number of cntYYNN sites is: 4628
        The number of cntNNYY sites is: 0
        The number of cntYYYY sites is: 189391
        The number of cntsame sites is: 175505
        The number of cntdiff sites is: 13886
         */
        
         
     }
     
     
     
     public void convertHotspotsEnd_AGPV2toAGPV4_cn_deprecated () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4start.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4.txt";
        String tempInputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4end.bed";
        String outputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4end.bed.con";
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz";
        RowTable<String> t = new RowTable (infileDirS);
        int[] chr = new int[t.getRowNumber()]; 
        int[] pos = new int[t.getRowNumber()]; 
        for (int i = 0; i < chr.length; i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 5));
        }
//        CrossMapUtils cm = new CrossMapUtils(chr, pos, tempInputbedFileS);
//        /*chain file download position: ftp://ftp.ensemblgenomes.org/pub/plants/release-40/assembly_chain/zea_mays/ */
//        cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath, tempInputbedFileS, outputbedFileS );
//        cm.convert();
//        List<int[]> l = cm.getConvertedCoordinate();
        //cm.deleteBedFiles();
//        chr = l.get(0);
//        pos = l.get(1);
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.setCell(i, 5, String.valueOf(pos[i]));
            //t.setCell(i, 0, String.valueOf(chr[i]));
        }
        String[] chrS = new String[chr.length];
        for (int i=0; i< chr.length; i++){
            chrS[i] = String.valueOf(chr[i]);
        }
        t.insertColumn("chr", 6, Arrays.asList(chrS));
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
     
     public void convertHotspotsEnd_AGPV2toAGPV4_deprecated () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4start.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4.txt";
        String tempInputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4end.bed";
        String outputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4end.bed.con";
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz";
        RowTable<String> t = new RowTable (infileDirS);
        int[] chr = new int[t.getRowNumber()]; 
        int[] pos = new int[t.getRowNumber()]; 
        for (int i = 0; i < chr.length; i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 5));
        }
//        CrossMapUtils cm = new CrossMapUtils(chr, pos, tempInputbedFileS);
//        /*chain file download position: ftp://ftp.ensemblgenomes.org/pub/plants/release-40/assembly_chain/zea_mays/ */
//        cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath, tempInputbedFileS, outputbedFileS );
//        cm.convert();
//        List<int[]> l = cm.getConvertedCoordinate();
        //cm.deleteBedFiles();
//        chr = l.get(0);
//        pos = l.get(1);
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.setCell(i, 5, String.valueOf(pos[i]));
            //t.setCell(i, 0, String.valueOf(chr[i]));
        }
        String[] chrS = new String[chr.length];
        for (int i=0; i< chr.length; i++){
            chrS[i] = String.valueOf(chr[i]);
        }
        t.insertColumn("chr", 6, Arrays.asList(chrS));
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
     public void convertHotspotsStart_AGPV2toAGPV4_cn_deprecated () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv2/RodgersMelnick2015PNAS_cnnamImputedXOsegments.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4start.txt";
        String tempInputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4start.bed";
        String outputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_cnnamImputedXOsegments_agpv4start.bed.con";
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz";
        RowTable<String> t = new RowTable (infileDirS);
        int[] chr = new int[t.getRowNumber()]; 
        int[] pos = new int[t.getRowNumber()]; 
        for (int i = 0; i < chr.length; i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 4));
        }
//        CrossMapUtils cm = new CrossMapUtils(chr, pos, tempInputbedFileS);
//        /*chain file download position: ftp://ftp.ensemblgenomes.org/pub/plants/release-40/assembly_chain/zea_mays/ */
//        cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath, tempInputbedFileS, outputbedFileS );
//        cm.convert();
//        List<int[]> l = cm.getConvertedCoordinate();
        //cm.deleteBedFiles();
//        chr = l.get(0);
//        pos = l.get(1);
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.setCell(i, 0, String.valueOf(chr[i]));
            t.setCell(i, 4, String.valueOf(pos[i]));
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
     
     public void convertHotspotsStart_AGPV2toAGPV4_deprecated () {
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv2/RodgersMelnick2015PNAS_usnamImputedXOsegments.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/source_agpv4/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4start.txt";
        String tempInputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4start.bed";
        String outputbedFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/002recombination/temp/RodgersMelnick2015PNAS_usnamImputedXOsegments_agpv4start.bed.con";
        String myPythonPath = "/Users/Aoyue/miniconda3/bin/python";
        String myCrossMapPath = "/Users/Aoyue/miniconda3/bin/CrossMap.py";
        String myMaizeChainPath = "/Users/Aoyue/Documents/Data/referenceGenome/crossMap/AGPv2_to_AGPv4.chain.gz";
        RowTable<String> t = new RowTable (infileDirS);
        int[] chr = new int[t.getRowNumber()]; 
        int[] pos = new int[t.getRowNumber()]; 
        for (int i = 0; i < chr.length; i++) {
            chr[i] = Integer.parseInt(t.getCell(i, 0));
            pos[i] = Integer.parseInt(t.getCell(i, 4));
        }
//        CrossMapUtils cm = new CrossMapUtils(chr, pos, tempInputbedFileS);
//        /*chain file download position: ftp://ftp.ensemblgenomes.org/pub/plants/release-40/assembly_chain/zea_mays/ */
//        cm = new CrossMapUtils(myPythonPath, myCrossMapPath, myMaizeChainPath, tempInputbedFileS, outputbedFileS );
//        cm.convert();
//        List<int[]> l = cm.getConvertedCoordinate();
        //cm.deleteBedFiles();
//        chr = l.get(0);
//        pos = l.get(1);
        for (int i = 0; i < t.getRowNumber(); i++) {
            t.setCell(i, 0, String.valueOf(chr[i]));
            t.setCell(i, 4, String.valueOf(pos[i]));
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
}
