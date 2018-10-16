/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import com.google.common.collect.Table;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Expression {

    public Expression() {
        //this.processRawPopBase();
        this.conGeneModelV3toV4();
        //this.mergeConfidenceGene();
        
    }
    
    private void mergeConfidenceGene () {
        String inputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression";
        String highConfidenceFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/004_snpclass/highConfidence_transcript.txt";
        String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/003_delAndExpression";
        new File (outputDirS).mkdir();
        RowTable t = new RowTable (highConfidenceFileS);
        HashMap<String, Double> geneSynMap = new HashMap();
        HashMap<String, Double> geneDelMap = new HashMap();
        HashMap<String, Double> geneRioMap = new HashMap();
        HashMap<String, Double> geneGerpMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 4).equals("0")) continue; //判断该基因是否有sift值，若没有，0，跳过
            String name = t.getCellAsString(i, 0);
            if (name.startsWith("GRM")) {
                name = name.split("_")[0];
            }
            double del = Double.valueOf(t.getCellAsDouble(i, 11));
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
//        for (int i = 0; i < t.getRowNumber(); i++) {
//            if (t.content[i][4].equals("0")) continue;
//            String name = t.content[i][0];
//            if (name.startsWith("GRM")) {
//                name = name.split("_")[0];
//            }
//            double del = Double.valueOf(t.content[i][27]);
//            double v = Double.NaN;
//            if (Double.valueOf(t.content[i][23]) != 0) {
//                v = del/Double.valueOf(t.content[i][23]);
//            }
//            geneSynMap.put(name, Double.valueOf(t.content[i][25]));
//            geneDelMap.put(name, del);
//            geneRioMap.put(name, v);
//            geneGerpMap.put(name, Double.valueOf(t.content[i][16]));
//        }
        Set<String> geneSet = geneSynMap.keySet();
        File[] fs = new File (inputDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            try {
                BufferedReader br =IOUtils.getTextReader(fs[i].getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File (outputDirS, fs[i].getName()).getAbsolutePath());
                bw.write(br.readLine()+"\tSyn\tDel\tRatioDelVsSyn\tMeanGerp");
                bw.newLine();
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (!geneSet.contains(l.get(0))) continue; //如果geneSet不包含表达数据中的基因，就跳过。因此就错过了AC开头的
                    if (Double.isNaN(geneRioMap.get(l.get(0)))) continue; // 如果比率为0.也跳过。
                    if (Double.isNaN(Double.valueOf(l.get(2)))) continue; // 标准偏差为0，也跳过。
                    StringBuilder sb = new StringBuilder(temp);
                    sb.append("\t").append(geneSynMap.get(l.get(0))).append("\t").append(geneDelMap.get(l.get(0))).append("\t").append(geneRioMap.get(l.get(0)));
                    sb.append("\t").append(geneGerpMap.get(l.get(0)));
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
    }
    
    private void conGeneModelV3toV4(){
        String infileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/001_V3";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/002_V3toV4GeneModel";
        String v3v4geneFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression/source/GeneModelv3v4.txt";
        RowTable<String> t = new RowTable<>(v3v4geneFileS);
        HashMap<String, String> hm = new HashMap<>();
        for(int i =0; i< t.getRowNumber();i++){
            String V4geneModel = t.getCellAsString(i, 0);
            String V3geneModel = t.getCellAsString(i, 1);
            hm.put(V3geneModel, V4geneModel);
        }
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "txt");
        for(int i =0; i<fs.length; i++){
            String infileS = fs[i].getAbsolutePath();
            String outfileS = new File(outfileDirS,fs[i].getName()).getAbsolutePath();
            RowTable<String> q = new RowTable<>(infileS);
            for(int j=0; j<q.getRowNumber(); j++){
                String V3gene = q.getCellAsString(j, 0);
                String V4gene = hm.get(V3gene);
                q.setCell(j, 0, V4gene);
            }
            q.writeTextTable(outfileS, IOFileFormat.Text);
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


        }
    }
    
    private void processRawPopBase () {
        String inputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/001_source/flat";
        String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/003_expression/002_processedExpression";
        File[] fs = new File(inputDirS).listFiles();
        for (int i = 0; i <  fs.length; i++) {
            try {
                RowTable<String> t = new RowTable (fs[i].getAbsolutePath());
                String outfileS = new File (outputDirS, fs[i].getName()).getAbsolutePath();
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
                boolean[] isthere = new boolean[genes.length];
                int[] indices = new int[genes.length];
                for (int j = 0; j < t.getRowNumber(); j++) {
                    String name = t.getCellAsString(j, 0);
                    if (name.startsWith("GRM")) {
                        name = name.split("_")[0];
                    }
                    int index = Arrays.binarySearch(genes, name);
                    if (isthere[index] == true) continue; //默认是都等于false
                    isthere[index] = true; //把在genes库中第几个搜索到，就把第几个index赋值为true
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
        double sd = ds.getStandardDeviation(); //标准偏差
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
    
}
