/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 *
 * @author yxh
 */
public class Hisat2Alignment {
    
    public Hisat2Alignment(){
//        this.siglealignment();
        this.doublealignment();
    }
    
    private void siglealignment(){
        String genomeindex = "/data1/home/xiaohan/rnaseq/reference/ABDtest/ABD/genome";
        String fastqpath = "/data1/home/xiaohan/rnaseq/coleoptile/splitedtest";
        String samoutputpath = "/data1/home/xiaohan/rnaseq/coleoptile/BAMhisat";
        long startTimePoint = System.nanoTime();
        String subfastq = new File (fastqpath).getAbsolutePath();
        File[] fs = new File(subfastq).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "R2.fq");
        List<File> fList = Arrays.asList(fs);
        List<String> fileList=new ArrayList<>();
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fileList.add(fList.get(i).getName().replace("_R2.fq", ""));
        }
        fileList.stream().forEach(p -> {
            StringBuilder sb = new StringBuilder();
            sb.append("hisat2 ").append("-t ").append("-x ").append(genomeindex+" ");
            sb.append("-U ").append(new File(fastqpath,p+"_R2.fq "));
            sb.append("-S ").append(samoutputpath+"/").append(p).append(".sam");
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
    
    private void doublealignment(){
        String genomeindex = "/data1/home/xiaohan/rnaseq/reference/ABDtest/ABD/genome";
        String fastqpath = "/data1/home/xiaohan/rnaseq/root/splitedtest";
        String samoutputpath = "/data1/home/xiaohan/rnaseq/root/hisat2double";
        long startTimePoint = System.nanoTime();
        String subfastq = new File (fastqpath).getAbsolutePath();
        File[] fs = new File(subfastq).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "R2.fq");
        List<File> fList = Arrays.asList(fs);
        List<String> fileList=new ArrayList<>();
        for(int i=0;i<fs.length;i++){
            if (fs[i].length() < 6_500_000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            fileList.add(fList.get(i).getName().replace("_R2.fq", ""));
        }
        fileList.stream().forEach(p -> {
            StringBuilder sb = new StringBuilder();
            sb.append("hisat2 ").append("-t ").append("-x ").append(genomeindex+" ");
            sb.append("-1 ").append(new File(fastqpath,p+"_R1.fq "));
            sb.append("-2 ").append(new File(fastqpath,p+"_R2.fq "));
            sb.append("-S ").append(samoutputpath+"/").append(p).append(".sam");
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
        new Hisat2Alignment();
    }
}
