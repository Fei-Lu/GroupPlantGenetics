/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import com.google.common.collect.Table;
import static com.google.common.io.Files.map;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Expression {

    public Expression() {
        //this.checkHeader();
        //this.processRawPopBase();
        //this.method2();
        
        //this.conGeneModelV2toV4_deprecated();
        //this.conGeneModelV4toV2();
        //this.mergeConfidenceGene();
        //this.mergeConfidenceGene_phospho();
        //this.scatterData();
        //this.deleteriousAndGeneByTissue();
        
        /**
         * 把老师的文件重新处理一遍，看看最终结果是否一致。探究是否由于基因过滤引起的  非组织特异表达，含有有害突变多。
         */
        
        //this.FeimergeConfidenceGene();
        this.FeimergeConfidenceGene_phospho();
        this.FeiscatterData();
        
    }
    
    public void FeiscatterData () {
        String tranFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/transcriptome_withSift.txt";
        String proFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/proteome_withSift.txt";
        String phoFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/phosphoproteome_withSift.txt";
        String tranBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/003_scatterplot/transcriptome.bin.txt";
        String proBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/003_scatterplot/proteome.bin.txt";
        String phoBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/003_scatterplot/phosphoproteome.bin.txt";
        int binSize = 200;
        this.subBinData(tranFileS, tranBinFileS, binSize);
        this.subBinData(proFileS, proBinFileS, binSize);
        this.subBinData(phoFileS, phoBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(tranFileS, tranBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(proFileS, proBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(phoFileS, phoBinFileS, binSize);
    }
    
    public void FeimergeConfidenceGene_phospho(){
        String infileS = "/Users/Aoyue/Documents/maf_fei/003_expression/processedExpression/phosphoproteome.txt";
        String highConfidenceFileS = "/Users/Aoyue/Documents/maf_fei/highConfidence_transcript.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/phosphoproteome_withSift.txt";
        RowTable t = new RowTable (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 4).equals("0")) continue; //判断该基因是否有sift值，若没有，0，跳过
            String name = t.getCellAsString(i, 0);
            if(name.startsWith("G")){
                name = t.getCellAsString(i, 0).split("_")[0];
            }
            else{
                    name = this.getGeneName(name);
            }
            double del = Double.valueOf(t.getCellAsDouble(i, 13));
            double syn = Double.valueOf(t.getCellAsDouble(i, 6));
            double gerp = Double.valueOf(t.getCellAsDouble(i, 16));
            double v = Double.NaN; //not a number
            if (syn != 0) {
                v = del/syn;
            }
            geneSynMap.put(name, syn);
            geneDelMap.put(name, del);
            geneRioMap.put(name, v);
            geneGerpMap.put(name, gerp);
        }
        
        Set<String> geneSet = geneSynMap.keySet();
        try {
            BufferedReader br =IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            String temp = null;
            int cnt =0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                String name = l.get(0);
                if(name.startsWith("A") || name.startsWith("E")){
                    cnt++;
                    name = this.getGeneNameP(name);
                }
                if (!geneSet.contains(name)) continue; //如果geneSet不包含表达数据中的基因，就跳过。因此就错过了AC开头的
                if (Double.isNaN(geneRioMap.get(name))) continue; // 如果比率为0,也跳过。
                if (Double.isNaN(Double.valueOf(l.get(2)))) continue; // 标准偏差为0,也跳过。
                int m = l.size()-1;
                StringBuilder sb = new StringBuilder();
                sb.append(name);sb.append("\t");
                for(int i = 1; i< l.size(); i++){
                    sb.append(l.get(i));sb.append("\t");
                }
                
                sb.append(geneSynMap.get(name)).append("\t").append(geneDelMap.get(name)).append("\t").append(geneRioMap.get(name))
                .append("\t").append(geneGerpMap.get(name));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(cnt + "  A or E");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    private void FeimergeConfidenceGene () {
        String infileS = "/Users/Aoyue/Documents/maf_fei/003_expression/processedExpression/transcriptome.txt";
        String highConfidenceFileS = "/Users/Aoyue/Documents/maf_fei/highConfidence_transcript.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/transcriptome_withSift.txt";
        String infile2S = "/Users/Aoyue/Documents/maf_fei/003_expression/processedExpression/proteome.txt";
        String outfile2S = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/000_20181108feiresult/003_delAndExpression/001_merge/proteome_withSift.txt";
        RowTable t = new RowTable (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 4).equals("0")) continue; //判断该基因是否有sift值，若没有，0，跳过
            String name = t.getCellAsString(i, 0);
            if(name.startsWith("G")){
                name = t.getCellAsString(i, 0).split("_")[0];
            }
            else{
                    name = this.getGeneName(name);
            }
            double del = Double.valueOf(t.getCellAsDouble(i, 13));
            double syn = Double.valueOf(t.getCellAsDouble(i, 6));
            double gerp = Double.valueOf(t.getCellAsDouble(i, 16));
            double v = Double.NaN; //not a number
            if (syn != 0) {
                v = del/syn;
            }
            geneSynMap.put(name, syn);
            geneDelMap.put(name, del);
            geneRioMap.put(name, v);
            geneGerpMap.put(name, gerp);
        }
        
        Set<String> geneSet = geneSynMap.keySet();
        try {
            BufferedReader br =IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                if (!geneSet.contains(l.get(0))) continue; //如果geneSet不包含表达数据中的基因，就跳过。因此就错过了AC开头的
                if (Double.isNaN(geneRioMap.get(l.get(0)))) continue; // 如果比率为0,也跳过。
                if (Double.isNaN(Double.valueOf(l.get(2)))) continue; // 标准偏差为0,也跳过。
                StringBuilder sb = new StringBuilder(temp);
                sb.append("\t").append(geneSynMap.get(l.get(0))).append("\t").append(geneDelMap.get(l.get(0))).append("\t").append(geneRioMap.get(l.get(0)))
                .append("\t").append(geneGerpMap.get(l.get(0)));
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
    
    public void method2(){
//        this.phosphoproteomeCheck();
//        this.sortPhospho();
//        this.addAverage();
//        this.phosphoDupli();
        //this.spiltPhospho();
        
    }
    
    
    
    /**
     * 待解释？
     */
    
    void scatterData () {
        String tranFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/transcriptome_withSift.txt";
        String proFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/proteome_withSift.txt";
        String phoFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/phosphoproteome_withSift.txt";
        String tranBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/003_scatterplot/transcriptome.bin.txt";
        String proBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/003_scatterplot/proteome.bin.txt";
        String phoBinFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/003_scatterplot/phosphoproteome.bin.txt";
        int binSize = 200;
        this.subBinData(tranFileS, tranBinFileS, binSize);
        this.subBinData(proFileS, proBinFileS, binSize);
        this.subBinData(phoFileS, phoBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(tranFileS, tranBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(proFileS, proBinFileS, binSize);
//        this.subBinDataSortByExpressionValue(phoFileS, phoBinFileS, binSize);
    }
    
    void subBinData (String inputFileS, String outputFileS, int binSize) {
        RowTable<String> t = new RowTable<> (inputFileS);
        int[][] bound = PArrayUtils.getSubsetsIndicesBySubsetSize(t.getRowNumber(), binSize);
        double[][][] exAndError = new double[bound.length][3][2];
        double[][][] tsAndError = new double[bound.length][3][2];
        double[][][] vsAndError = new double[bound.length][3][2];
        t.sortAsNumber("Syn");
        System.out.println(t.getCellAsString(t.getRowNumber()-1, 0));
        double[] ex = t.getColumnAsDoubleArray(1);
        double[] ts = t.getColumnAsDoubleArray(2);
        double[] vs = t.getColumnAsDoubleArray(t.getColumnIndex("Syn"));
        
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][0][0] = mean;
            exAndError[i][0][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][0][0] = mean;
            tsAndError[i][0][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][0][0] = mean;
            vsAndError[i][0][1] = se;
        }
        
        t.sortAsNumber("Del");
        System.out.println(t.getCellAsString(t.getRowNumber()-1, 0));
        ex = t.getColumnAsDoubleArray(1);
        ts = t.getColumnAsDoubleArray(2);
        vs = t.getColumnAsDoubleArray(t.getColumnIndex("Del"));
        
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][1][0] = mean;
            exAndError[i][1][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][1][0] = mean;
            tsAndError[i][1][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][1][0] = mean;
            vsAndError[i][1][1] = se;
        }
        
        t.sortAsNumber("RatioDelVsSyn");
        System.out.println(t.getCellAsString(t.getRowNumber()-1, 0));
        ex = t.getColumnAsDoubleArray(1);
        ts = t.getColumnAsDoubleArray(2);
        vs = t.getColumnAsDoubleArray(t.getColumnIndex("RatioDelVsSyn"));
        for (int i = 0; i <  bound.length; i++) {
            TDoubleArrayList exList = new TDoubleArrayList();
            TDoubleArrayList tsList = new TDoubleArrayList();
            TDoubleArrayList vsList = new TDoubleArrayList();
            for (int j = bound[i][0]; j < bound[i][1]; j++) {
                exList.add(ex[j]);
                tsList.add(ts[j]);
                vsList.add(vs[j]);
            }
            double[] exSub = exList.toArray();
            double[] tsSub = tsList.toArray();
            double[] vsSub = vsList.toArray();
            DescriptiveStatistics d = new DescriptiveStatistics(exSub);
            double mean = d.getMean();
            double sd = d.getStandardDeviation();
            double se = sd/Math.sqrt(exSub.length);
            exAndError[i][2][0] = mean;
            exAndError[i][2][1] = se;
            d = new DescriptiveStatistics(tsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            tsAndError[i][2][0] = mean;
            tsAndError[i][2][1] = se;
            d = new DescriptiveStatistics(vsSub);
            mean = d.getMean();
            sd = d.getStandardDeviation();
            se = sd/Math.sqrt(exSub.length);
            vsAndError[i][2][0] = mean;
            vsAndError[i][2][1] = se;
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            bw.write("Group rank\tExpressSynMean\tExpressSynError\tExpressDelMean\tExpressDelError\tExpressRatioDVsSMean\tExpressRatioDVsSError\tTsSynMean\tTsSynError\tTsDelMean\tTsDelError\tTsRatioDVsSMean\tTsRatioDVsSError\tSynMean\tSynError\tDelMean\tDelError\tRatioDVsSMean\tRatioDVsSError");
            bw.newLine();
            for (int i = 0; i < bound.length-1; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(i+1);
                for (int j = 0; j < exAndError[i].length; j++) {
                    for (int k = 0; k < exAndError[i][j].length; k++) {
                        sb.append("\t").append(exAndError[i][j][k]);
                    }
                }
                for (int j = 0; j < tsAndError[i].length; j++) {
                    for (int k = 0; k < tsAndError[i][j].length; k++) {
                        sb.append("\t").append(tsAndError[i][j][k]);
                    }
                }
                for (int j = 0; j < vsAndError[i].length; j++) {
                    for (int k = 0; k < vsAndError[i][j].length; k++) {
                        sb.append("\t").append(vsAndError[i][j][k]);
                    }
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
    
    public void mergeConfidenceGene_phospho(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/001_V2/phosphoproteome.txt";
        String highConfidenceFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/highConfidence_geneV2.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/phosphoproteome_withSift.txt";
        RowTable t = new RowTable (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 4).equals("0")) continue; //判断该基因是否有sift值，若没有，0，跳过
            String name = t.getCellAsString(i, 0);
            double del = Double.valueOf(t.getCellAsDouble(i, 13));
            double syn = Double.valueOf(t.getCellAsDouble(i, 6));
            double gerp = Double.valueOf(t.getCellAsDouble(i, 16));
            double v = Double.NaN; //not a number
            if (syn != 0) {
                v = del/syn;
            }
            geneSynMap.put(name, syn);
            geneDelMap.put(name, del);
            geneRioMap.put(name, v);
            geneGerpMap.put(name, gerp);
        }
        
        Set<String> geneSet = geneSynMap.keySet();
        try {
            BufferedReader br =IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            String temp = null;
            int cnt =0;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                String name = l.get(0);
                if(name.startsWith("A") || name.startsWith("E")){
                    cnt++;
                    name = this.getGeneNameP(name);
                }
                if (!geneSet.contains(name)) continue; //如果geneSet不包含表达数据中的基因，就跳过。因此就错过了AC开头的
                if (Double.isNaN(geneRioMap.get(name))) continue; // 如果比率为0,也跳过。
                if (Double.isNaN(Double.valueOf(l.get(2)))) continue; // 标准偏差为0,也跳过。
                int m = l.size()-1;
                StringBuilder sb = new StringBuilder();
                sb.append(name);sb.append("\t");
                for(int i = 1; i< l.size(); i++){
                    sb.append(l.get(i));sb.append("\t");
                }
                
                sb.append(geneSynMap.get(name)).append("\t").append(geneDelMap.get(name)).append("\t").append(geneRioMap.get(name))
                .append("\t").append(geneGerpMap.get(name));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(cnt + "  A or E");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
    }
    
    String getGeneNameP(String a){
        //String a = "AC148167.6_FGP001";
        int c = a.indexOf("P");
        char[] d = a.toCharArray();
        List<Character> l = new ArrayList<>();
        for(int i = 0 ; i < d.length; i++){
            if(i == c) continue;
            l.add(d[i]);
        }
        String gene = StringUtils.join(l, "");
        System.out.println(gene);
        return(gene);
    }
    
    
    String getGeneName(String a){
        //String a = "AC148167.6_FGP001";
        int c = a.indexOf("T");
        char[] d = a.toCharArray();
        List<Character> l = new ArrayList<>();
        for(int i = 0 ; i < d.length; i++){
            if(i == c) continue;
            l.add(d[i]);
        }
        String gene = StringUtils.join(l, "");
        System.out.println(gene);
        return(gene);
    }
    
    private void mergeConfidenceGene () {
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/001_V2/transcriptome.txt";
        String highConfidenceFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/highConfidence_geneV2.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/transcriptome_withSift.txt";
        String infile2S = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/001_V2/proteome.txt";
        String outfile2S = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression/001_V2merge/proteome_withSift.txt";
        RowTable t = new RowTable (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 4).equals("0")) continue; //判断该基因是否有sift值，若没有，0，跳过
            String name = t.getCellAsString(i, 0);
            double del = Double.valueOf(t.getCellAsDouble(i, 13));
            double syn = Double.valueOf(t.getCellAsDouble(i, 6));
            double gerp = Double.valueOf(t.getCellAsDouble(i, 16));
            double v = Double.NaN; //not a number
            if (syn != 0) {
                v = del/syn;
            }
            geneSynMap.put(name, syn);
            geneDelMap.put(name, del);
            geneRioMap.put(name, v);
            geneGerpMap.put(name, gerp);
        }
        
        Set<String> geneSet = geneSynMap.keySet();
        try {
            BufferedReader br =IOUtils.getTextReader(infile2S);
            BufferedWriter bw = IOUtils.getTextWriter(outfile2S);
            bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
            bw.newLine();
            String temp = null;
            while ((temp = br.readLine()) != null) {
                List<String> l = PStringUtils.fastSplit(temp);
                if (!geneSet.contains(l.get(0))) continue; //如果geneSet不包含表达数据中的基因，就跳过。因此就错过了AC开头的
                if (Double.isNaN(geneRioMap.get(l.get(0)))) continue; // 如果比率为0,也跳过。
                if (Double.isNaN(Double.valueOf(l.get(2)))) continue; // 标准偏差为0,也跳过。
                StringBuilder sb = new StringBuilder(temp);
                sb.append("\t").append(geneSynMap.get(l.get(0))).append("\t").append(geneDelMap.get(l.get(0))).append("\t").append(geneRioMap.get(l.get(0)))
                .append("\t").append(geneGerpMap.get(l.get(0)));
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
     * 将highConfidence_transcript.txt 转化为 V2版本的gene model
     */
    
    public void conGeneModelV4toV2(){
        String dbfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/V4V2geneModel.txt";
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/highConfidence_geneV2.txt";
        RowTable<String> t = new RowTable<>(dbfileS);
        HashMap<String, String> hm = new HashMap<>();
        List<String> V4 = new ArrayList<>();
        int cnt = 0;
        for(int i =0; i< t.getRowNumber();i++){
            String V4geneModel = t.getCellAsString(i, 0);
            String V2geneModel = t.getCellAsString(i, 1);
            if(V2geneModel.equals("")){
                cnt++;
                //System.out.println(V4geneModel);
            }
            else{
                hm.put(V4geneModel, V2geneModel);
                V4.add(V4geneModel);
            }
        }
        System.out.println(cnt + "\tV4 has no V2"); //58	V4 has no V2
        String[] db = V4.toArray(new String[V4.size()]);
        Arrays.sort(db);
        System.out.println(db.length + "\tkey-value pairs"); //29370	key-value pairs
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS); //共有26016lines
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            int cntl =0;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                int m = l.size()-1;
                String query = l.get(0).split("_")[0];
                int index = Arrays.binarySearch(db, query);
                if (index < 0) continue;
                cntl++;
                bw.write(hm.get(query));bw.write("\t");
                for(int i = 1; i< m; i++){
                    bw.write(l.get(i));bw.write("\t");
                }
                bw.write(l.get(m));bw.newLine();
                
            }
            br.close();bw.flush();bw.close(); 
            System.out.println(cntl + "\tgenes are in high confidence"); //20595	genes are in high confidence
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void conGeneModelV2toV4_deprecated(){
        String v4v2geneFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/V4V2geneModel.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/V4geneModelwithNoV2.txt";
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String outfile2S = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/V4geneModelwithNoV2.txt";
        
        /**
         * 1.建立v2v4map,并且找到 #有V4没有V2# 的gene model，写入文件outfileS；
         * 2.将 #有V4没有V2# 转化为数组 query， 到high confi 文件中去搜索；58个
         */
        RowTable<String> t = new RowTable<>(v4v2geneFileS);
        HashMap<String, String> hm = new HashMap<>();
        List<String> lV4noV2 = new ArrayList<>();
        int cnt = 0;
        String temp = null;
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("V4genemodelwithNoV2");bw.newLine();
            for(int i =0; i< t.getRowNumber();i++){
                String V4geneModel = t.getCellAsString(i, 0);
                String V2geneModel = t.getCellAsString(i, 1);
                hm.put(V4geneModel, V2geneModel);
                if(V2geneModel.equals("")){
                    cnt++;
                    bw.write(V4geneModel);bw.newLine();
                    //System.out.println(V4geneModel);
                    lV4noV2.add(V4geneModel);
                }
            }
            
            System.out.println(cnt + "\tV4 has no V2");
            String[] query = lV4noV2.toArray(new String[lV4noV2.size()]);
                    
            /**
             * 1.在high confi 文件中找出不能转化的V4基因，一共26016lines。 41个
             * 存入数组V4noV2 中；
             */
            
            t = new RowTable<>(infileS);
            List<String> V4 = new ArrayList<>();
            int cntquery = 0;
            for(int i=0; i<t.getRowNumber();i++){
                String geneV4 = t.getCellAsString(i, 0).split("_")[0];
                V4.add(geneV4);
            }
            String[] V4db = V4.toArray(new String[V4.size()]);
            Arrays.sort(V4db);
            List<String> V4cannot = new ArrayList<>();
            for(int i = 0; i < query.length; i++){
                int index = Arrays.binarySearch(V4db, query[i]);
                if (index > 0){
                    V4cannot.add(query[i]);
                    cntquery++;
                    System.out.println(query[i]);
                }
            }
            System.out.println(cntquery + "\tcan not convert from V4 to V2");
            System.out.println(V4cannot.size());
            String[] V4noV2 = V4cannot.toArray(new String[V4cannot.size()]);
            Arrays.sort(V4noV2);
            
            /**
             * 
             */
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfile2S);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                int m = l.size()-1;
                String name = l.get(0).split("_")[0];
                String V2name = hm.get(name);
                //if(V2name.equals("")) continue;
                int index = Arrays.binarySearch(V4noV2, name);
                
                if (index <0 ) {
                    bw.write(V2name);bw.write("\t");
                    for(int i = 1; i< m; i++){
                        bw.write(l.get(i));bw.write("\t");
                    }
                    bw.write(l.get(m));bw.newLine();
                    }
            }
            bw.flush();bw.close();
            br.close();bw2.flush();bw2.close();
        }
        catch(Exception e){
            System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
        


        
     
        

        
//        File[] fs = new File(infileDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "txt");
//        for(int i =0; i<fs.length; i++){
//            String infileS = fs[i].getAbsolutePath();
//            String outfileS = new File(outfileDirS,fs[i].getName()).getAbsolutePath();
//            RowTable<String> q = new RowTable<>(infileS);
//            for(int j=0; j<q.getRowNumber(); j++){
//                String V3gene = q.getCellAsString(j, 0);
//                String V4gene = hm.get(V3gene);
//                q.setCell(j, 0, V4gene);
//            }
//            q.writeTextTable(outfileS, IOFileFormat.Text);
//            try{
//            BufferedReader br = IOUtils.getTextReader(infileS);
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            String header = br.readLine();
//            bw.write(header);bw.newLine();
//            String temp = null;
//            while((temp = br.readLine()) != null){
//                String key = PStringUtils.fastSplit(temp).get(0);
//                bw.newLine();
//            }
//            bw.flush();bw.close();br.close();
//            }
//            catch(Exception e){
//                //System.out.println(temp);
//                e.printStackTrace();
//                System.exit(1);
//
//            }
//
//
//        }
    }
    
    /**
     * 这一步的主要目的：
     * 1.统一第一列的基因名，将GRMZM2G123896_P04开头的，转化成GRMZM2G123896;并排序，按照AC AF AY EF （前几个，一个代表一个基因）GRMZM..._T.00.. 此开头，表示基因的第几个转录本,是重复的。
     *  如果还有4个重复，则选取一个？
     * 2.在每一行基因信息里，加3列信息：平均值、相对标准偏差、数组为0的个数。
     * 该方法思路复杂，请参照method2，弃用该方法。
     */
    
    private void processRawPopBase () {
        String inputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat";
        String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/001_V2";
        File[] fs = new File(inputDirS).listFiles();
        for (int i = 0; i <  fs.length; i++) {
            try {
                RowTable<String> t = new RowTable (fs[i].getAbsolutePath());
                String outfileS = new File (outputDirS, fs[i].getName().replaceFirst(".", "_V2.")).getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder("Gene\tMean\tRSD\tZeroCount");
                for (int k = 1; k < t.getHeader().size(); k++) {
                    sb.append("\t").append(t.getHeader().get(k));
                }
                bw.write(sb.toString());
                bw.newLine();
                HashSet<String> geneSet = new HashSet();
                
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.getCellAsString(j, 0);
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    geneSet.add(name);
                }
                String[] genes = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genes);
                //boolean[] isthere = new boolean[genes.length]; //isthere是什么意思？
                int[] indices = new int[genes.length]; //索引
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.getCellAsString(j, 0);
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    int index = Arrays.binarySearch(genes, name);
                    //if (isthere[index] == true) continue; //默认是都等于false，等于true时，跳过。
                    //isthere[index] = true; //把在genes库中第几个搜索到，就把第几个index赋值为true
                    indices[index] = j; //j等于0， index = 6的时候，说明第一行的基因在库里的第6个找到了，库6等于0。j等于8， index = 0的时候，说明第8行的基因在库里的第0个找到了，库0等于8。
                }
                for (int j = 0; j < genes.length; j++) { //对gene库的每个基因进行处理
                    double[] value = new double[t.getColumnNumber()-1]; // 有多少个变量值！！把组织表达的数值存储到一个数组value中
                    for (int k = 0; k < value.length; k++) { //对每一个gene的每一个value进行处理
                        int m = k +1; //要想知道这个基因的第k个值，就要知道这个基因在表格里的第几行，第几列！
                        value[k] = Double.valueOf(t.getCellAsString(indices[j], m)); //将String类型的值转化为double类型，库0对应的基因所在的行数，
                    }
                    bw.write(this.getValueString(genes[j], value)); // 需要输入 基因名字，基因的每一列值所在的数组
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }   
        }
    }  
    
    String getValueString (String gene, double[] value) {
        StringBuilder sb = new StringBuilder(gene);
        DescriptiveStatistics ds = new DescriptiveStatistics(value); //把数组名字输入进去
        double mean = ds.getMean(); //求平均数
        double sd = ds.getStandardDeviation(); //相对标准偏差，求这个有何意义？？是多个组织的相对标准偏差。
        int cnt = 0;
        for (int i = 0; i < value.length; i++) { //求该数组中为0的个数
            if (value[i] == 0) {
                cnt++;
            }
        }
        sb.append("\t").append(mean).append("\t").append(sd/mean).append("\t").append(cnt); //标准偏差除以平均数=变异系数,
        for (int i = 0; i < value.length; i++) {
            sb.append("\t").append(value[i]);
        }
        return sb.toString();
    }
    
    public void spiltPhospho(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_sortbyGene_addAverage_removedupli.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_sortbyGene_addAverage_removedupli_splitP.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            List<String> l = new ArrayList<>();
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp);
                int m = l.size()-1;
                String name = l.get(0);
                if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                bw.write(name);bw.write("\t");
                for(int i = 1; i< m; i++){
                    bw.write(l.get(i));bw.write("\t");
                }
                bw.write(l.get(m));bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    /**
     * 根据重复的基因列表，找出所有重复的转录本 37*2 +4=78个，进而打印出来进行一一筛选，选择标准的基因model，去掉40个重复的。
     */
    
    public void phosphoDupli(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage/phosphoproteome_sortbyGene_addAverage.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage/phosphoproteome_sortbyGene_addAverage_checkDupli.txt";
        String DupliS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage/phosphoproteomeDuplicata38.txt";
        RowTable<String> t = new RowTable<>(DupliS);
        
        List<String> l = t.getColumn(0);
        String[] genes = l.toArray(new String[l.size()]);
        Arrays.sort(genes);
        t = new RowTable<>(infileS);
        boolean[] ifout =new boolean[t.getRowNumber()];
        for(int i=0;i<t.getRowNumber();i++){
            String query = t.getCellAsString(i, 0).split("_")[0];
            int index = Arrays.binarySearch(genes, query);
            if (index >= 0){
                ifout[i] = true;
            }
        }
        t.writeTextTable(outfileS, IOFileFormat.Text, ifout);
    }
    
    public void addAverage(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/sort";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/addAverage";
        File[] fs = new File(infileDirS).listFiles();
        for(int i=0; i<fs.length;i++){
            try{
                String infileS = fs[i].getAbsolutePath();
                String outfileS = new File(outfileDirS, fs[i].getName().replaceFirst("Gene", "Gene_addAverage")).getAbsolutePath();
                RowTable<String> t = new RowTable<>(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder("Gene\tMean\tRSD\tZeroCount");
                for (int k = 1; k < t.getHeader().size(); k++) {
                    sb.append("\t").append(t.getHeader().get(k));
                }
                bw.write(sb.toString());
                bw.newLine();

                double[] value = new double[t.getColumnNumber()-1];
                for(int e =0; e <t.getRowNumber(); e++){
                    for(int k =0; k<value.length;k++){
                        int j = k+1;
                        value[k] = t.getCellAsDouble(e, j);
                    }
                    bw.write(this.getValueString(t.getCellAsString(e, 0), value));
                    bw.newLine();
                }
                bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        }
    }
    
    public void sortPhospho(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_removelist_Duplicata40.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome_removelist_Dupli40.txt";
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    /**
     * 把phospho变成list，再转成set，然后查看基因重复的转录本
     */
    public void phosphoproteomeCheck(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat/phosphoproteome.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(0);
        Collections.sort(l);
        List<String> lcuttrans = new ArrayList<>();
        for (int i =0; i<l.size(); i++){
            String trans = l.get(i);
            if(trans.startsWith("GRMZM")){
                String gene = trans.split("_")[0];
                lcuttrans.add(gene);
            }
            else{
                lcuttrans.add(trans);
            }
        }
        System.out.println(lcuttrans.size());
        Set<String> s = new HashSet<>(lcuttrans);
        System.out.println(s.size());
        System.out.println(s);
        int cnt = 0;
        for(String a : s){
            //System.out.println(a + "    " + Collections.frequency(lcuttrans, a));
            if(Collections.frequency(lcuttrans, a) == 2){
                System.out.println(a);
                cnt++;
            }
            else{
                if(Collections.frequency(lcuttrans, a) == 4){
                System.out.println(a + "\t4444444");
                // GRMZM2G179677	4444444
                cnt++;
                }
            }
        }
        System.out.println(cnt);
    }
    
    public void checkHeader(){
        String infileDirs = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat";
        File[] fs = new File (infileDirs).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "txt");
        for(int i = 0; i< fs.length; i++){
            String infileS = fs[i].getAbsolutePath();
            RowTable t = new RowTable (infileS);
            List l = t.getHeader();
            System.out.println(l.size() + "\tEs\t" + fs[i].getName().replaceFirst(".txt", ""));
            System.out.println(l);
        }
        System.out.println("Here");
    }
}
