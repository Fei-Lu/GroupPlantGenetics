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
public class HTseqcount {

    public HTseqcount(){
        this.count();
    }

    private void count(){
        long startTimePoint = System.nanoTime();
        String input = "/data1/home/xiaohan/rnaseq/root/STARalignment";
        String inputDirS = new File(input).getAbsolutePath();
        String GTFDirS = "/data1/home/xiaohan/rnaseq/reference/ABD/annotation/wheat_v1.1_Lulab.gtf";
        File[] fs= new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".samAligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        List<String> fileList=new ArrayList<>();
        for(int i=0;i<fs.length;i++){
            fileList.add(fList.get(i).getName());
        }
        fileList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count ").append("-f sam -r name -s no ").append(" -m intersection-nonempty -s reverse ");
            sb.append(inputDirS+"/"+f+" ");
            sb.append(GTFDirS);
            sb.append(">");
            sb.append(f.replace(".samAligned.out.sam", ".count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File ("/data1/home/xiaohan/rnaseq/root/HTseqcount").getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Runtime rt = Runtime.getRuntime();
                Process pro = rt.exec(cmdarry,null,dir);
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
            System.out.println("Finished"+f);
        });
        StringBuilder time = new StringBuilder();
        time.append("HTseqcount according to sam.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    
    public static void main(String[] args){
        new HTseqcount();
    }
}
    

