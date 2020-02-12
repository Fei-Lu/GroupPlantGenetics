/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author yxh
 */
public class STARalignment {
     
    public STARalignment(){
        this.alignment();
    }
    
    private void alignment(){
        String genomeindex = "/data1/home/xiaohan/rnaseq/reference/ABD/index";
//        String fastqpath = "/Users/yxh/Documents/RNA-seq/test/staralignementtest/fastq";
        String fastqpath = "/data1/home/xiaohan/rnaseq/root/splitedtest";
        String samoutputpath = "/data1/home/xiaohan/rnaseq/root/STARalignment";
        long startTimePoint = System.nanoTime();
        String subfastq = new File (fastqpath).getAbsolutePath();
        File[] fs = new File(subfastq).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq");
        List<File> fList = Arrays.asList(fs);
        List<String> fileList=new ArrayList<>();
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fileList.add(fList.get(i).getName().replace("_R1.fq", ""));
        }
        fileList.stream().forEach(p -> {
            StringBuilder sb = new StringBuilder();
            sb.append("STAR ").append("--runThreadN ").append("4 ").append("--genomeDir"+" ");
//            sb.append(genomeindex+"/SA ").append("--readFilesIn"+" ").append(fastqpath+"/"+p+"_R1.fq "+fastqpath+"/"+p+"_R2.fq");
            sb.append(genomeindex+" ").append("--readFilesIn"+" ").append(p+"_R1.fq "+p+"_R2.fq");
            sb.append(" --outFileNamePrefix ").append(new File(new File(samoutputpath,p+".sam").getAbsolutePath()));
            sb.append(" --outSAMtype SAM").append(" --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 80");
            String command = sb.toString();
            System.out.println(command);
            try {
                Runtime rt = Runtime.getRuntime();
                Process pro = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(pro.getInputStream()));
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    System.out.println(temp);
                }
                pro.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+p);
        });
        StringBuilder time = new StringBuilder();
        time.append("Alignment samples according to genome.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    
    public static void main(String[] args){
        new STARalignment();
    }
}


