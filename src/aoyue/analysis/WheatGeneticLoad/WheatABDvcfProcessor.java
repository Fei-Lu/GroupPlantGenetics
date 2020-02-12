/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class WheatABDvcfProcessor {

    public WheatABDvcfProcessor() {
        //this.subsetVCFdataRandom();
        new Treetest();
        
    }
    

    
    private void subsetVCFdataRandom () {
        String infileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001.vcf";
        String outfileS = "/data4/home/aoyue/vmap2/abd/rawVCF/chr001_subset.vcf";
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt =0;
            System.out.println(new SimpleDateFormat().format(new Date()) + "    program execution.\n");
            long startTime = System.nanoTime();
            while ((temp = br.readLine()) != null) {
                if(temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                else{
                    cnt++;
                    double r = Math.random();
                    if (r > 0.001) continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                    List<String> l = PStringUtils.fastSplit(temp);
                    if (l.get(3).contains(",")) continue; // 第3列是alt的信息，若有2个等位基因，则去除这一行
                    bw.write(temp);
                    bw.newLine();
                }
            }
            bw.flush();bw.close();br.close();
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");

            System.out.println("Chr 1 snp number is " + cnt + ".");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
