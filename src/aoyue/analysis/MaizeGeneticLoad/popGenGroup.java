/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import com.google.common.collect.Table;
import pgl.infra.table.RowTable;
import java.io.BufferedWriter;
import java.util.HashMap;
import java.util.Set;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class popGenGroup {

    public popGenGroup() {
        //this.getDepthOfHmp321();
        this.getHighDepthTaxaList();
        
        
    }
    
    /***** 结合hmp321.depth.txt，将测序深度大于5的taxa从文件geneticGroup中筛选出来，并建立高密度group ********/
    public void getHighDepthTaxaList () {
        String depthInfoFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/001_highDepthGroup/hmp321.depth.txt";
        String hmp321GeneticGroup = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/000_group/geneticGroup.txt";
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/001_highDepthGroup/hmp321.highDepth.group.txt";
        double minDepth = 5;
        RowTable t = new RowTable (depthInfoFileS);
        HashMap<String, Double> taxaDepthMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            if (t.getCellAsString(i, 2).startsWith("null")) continue;
            double depth = Double.valueOf(t.getCellAsString(i, 2));
            if (depth < minDepth) continue;
            taxaDepthMap.put(t.getCellAsString(i, 0), depth);
        }
        Set<String> keys = taxaDepthMap.keySet(); //把所有的key放入一个set中去
        t = new RowTable (hmp321GeneticGroup);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tPC1\tPC2\tGroupIndex\tGeneticGroup\tSequencingDepth");
            bw.newLine();
            for (int i = 0; i < t.getRowNumber(); i++) {
                String taxon = t.getCellAsString(i, 0);
                if (!keys.contains(taxon)) continue;
                StringBuilder sb = new StringBuilder();
                sb.append(taxon).append("\t").append(t.getCellAsString(i, 1)).append("\t").append(t.getCellAsString(i, 2)).append("\t").append(t.getCellAsString(i, 3)).append("\t").append(t.getCellAsString(i, 4)).append("\t").append(taxaDepthMap.get(taxon));
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
    
    
    /******************* 将原始Robert 的测序覆盖度 和 我们自己计算的平均SNP位点深度总结在一张表中 *********************************/
    public void getDepthOfHmp321 () {
        String robertDepthFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/001_highDepthGroup/source/depthByRobert.txt"; //每个taxa的原始测序coverage总结
        String siteDepthFileS = "/Users/Aoyue/Documents/maizeGeneticLoad/001_variantSummary/007_hmp321DeleCount/001_hmp321Depth/taxaDepth.summary.txt"; //每个taxa的平均深度
        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/001_highDepthGroup/hmp321.depth.txt";
        RowTable t = new RowTable (robertDepthFileS);
        HashMap<String, Double> robertMap = new HashMap();
        for (int i = 0; i < t.getRowNumber(); i++) {
            robertMap.put(t.getCellAsString(i, 0), Double.valueOf(t.getCellAsString(i, 1)));
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tSNPSiteDepth\tSequencingDepth");
            bw.newLine();
            t = new RowTable (siteDepthFileS);
            for (int i = 0; i < t.getRowNumber(); i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(t.getCellAsString(i, 0)).append("\t").append(t.getCellAsString(i, 2)).append("\t").append(robertMap.get(t.getCellAsString(i, 0)));
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
