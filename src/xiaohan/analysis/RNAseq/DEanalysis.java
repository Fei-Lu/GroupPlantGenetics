package xiaohan.analysis.RNAseq;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author yxh
 */
public class DEanalysis {
    public DEanalysis(){
        this.parseFqByIndexAndBarcode();
        this.STARalignment();
        this.Hisat2double();
        this.Hisat2single();
        this.HTseqcount();
        this.HTseqcountmerge();
        
        
    }
    
    private void parseFqByIndexAndBarcode () {
        HashMap barcodeStrain = new HashMap();
        HashMap indexStrain = new HashMap();
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rnaseq/RNA-seq-coleoptile-20180315-1.txt");
        RowTable<String> t = new RowTable<>("/Users/yxh/programs/RNA-seq-root-20180315-1.txt");
        
        List<String> barcodeList = new ArrayList<String>() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            barcodeStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0));
            barcodeList.add(t.getCell(i, 2).substring(1));
        }
        String inputDirS = "/Users/yxh/Documents/RNA-seq/test/splitedtest";
//        String inputDirS = "/data1/home/xiaohan/coleoptile/unsplitedtest";
        String outputDirS = "/Users/yxh/Documents/RNA-seq/test/splitedtest";
//        String outputDirS = "/data1/home/xiaohan/coleoptile/splitedtest";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".fastq.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.stream().forEach((String p) -> {      
            try {
                String infile1 = new File (inputDirS, p+"_R1_001.fastq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_R2_001.fastq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[barcodeList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[barcodeList.size()];
                for (int i = 0; i < barcodeList.size(); i++) {
                   bw[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R1.fq").getAbsolutePath());
                   bw1[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R2.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                int pos = -1 ;
                String temp = null;
                String seq = null;String index = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    index=temp.split(":")[9].substring(1,6);
                    if(indexStrain.get(index)!=null){
                        seq=br.readLine();
                        currentBarcode=seq.substring(1,8);
                        if(barcodeStrain.get(currentBarcode)!=null){
                            pos=barcodeList.indexOf(currentBarcode);
                            bw[pos].write(temp);bw[pos].newLine();
                            bw[pos].write(seq);bw[pos].newLine();
                            bw[pos].write(br.readLine());bw[pos].newLine();
                            bw[pos].write(br.readLine());bw[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());bw1[pos].newLine();
                        }else{
                            br.readLine();br.readLine();
                            br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                        }
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                    
                }          
                for(int i=0;i<barcodeList.size();i++){
                      bw[i].flush();bw[i].close();
                      bw1[i].flush();bw1[i].close();
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        StringBuilder time = new StringBuilder();
        time.append("Distinguish samples according to barcode and trim the barcode.").append("Took ").append(Benchmark.getTimeSpanSeconds(startTimePoint)).append(" seconds. Memory used: ").append(Benchmark.getUsedMemoryGb()).append(" Gb");
        System.out.println(time.toString());
    }
    public void HTseqcountmerge() {//顺序是HTSeq里面的顺序
        String DirS = "/Users/yxh/Documents/RNA-seq/006test/DEseq2/01-12";
//        String inputDirS = "/Users/yxh/Documents/RNA-seq/006test/HTseqcount/coleoptile";
//        String outputDirS = "/Users/yxh/Documents/RNA-seq/006test/DEseq2";
//        String inputDirS = "/data1/home/xiaohan/rnaseq/root/HTseqcount";
//        String outputDirS = "/data1/home/xiaohan/rnaseq/root/DESeq2";
//        String inputDirS = "/Users/yxh/Documents/RNA-seq/test/seq/coleoptile";
        List<String> nameList = new ArrayList<>();
        String subFqDirS = new File(DirS).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        Arrays.sort(fs);
        final List<File> fList = Arrays.asList(fs);
        int[][] count = new int[107891][fList.size()];
        
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < fList.size(); i++) {
            indexList.add(i);
        }
        indexList.stream().forEach(index -> {
            File f = fList.get(index);
            String temp = null;
            String[] tem = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                int cnt = 0;int counter = 0;
                while ((temp = br.readLine()) != null) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    if (tem[0].startsWith("TraesCS")) {
                        nameList.add(tem[0]);
                        int index1 = counter;
                        count[index1][index] = Integer.parseInt(tem[1]);
                        counter++;
                    }
                    cnt++;
                    if (cnt%10000 == 0) System.out.println(cnt);
                }
                br.close();
            } catch (Exception ex) {
                System.out.println(tem[0] + "\t1234");
                ex.printStackTrace();

            }
        });
        
                String outputFileS = new File(DirS, "count01-12.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            sb.append("Gene" + "\t");
            for (int i = 0; i < fList.size(); i++) {
                sb.append(fList.get(i).getName().replace(".txt", "") + "\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                sb = new StringBuilder();
                for (int j = 0; j < fList.size(); j++) {
                    if (j == 0) {
                        sb.append(nameList.get(i) + "\t");
                    }
                    sb.append(count[i][j] + "\t");
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

     private void HTseqcount(){
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
     
    private void Hisat2single(){
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
    
    private void Hisat2double(){
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
    private void STARalignment(){
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
        new DEanalysis();
    }
}


