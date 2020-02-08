/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Wheat200cleanDataProcessor {

    public Wheat200cleanDataProcessor() {
        this.GetIDlist();
        this.mergeIDlist();
        this.calFastqbp();
        //this.sampleFastQC();
        this.fastQC();
        this.mkParameterchr1_42();
        this.mkJavaCmdchr1_42();

    }
    
    /**
     * 本方法的目的是：建立42条染色体的java 运行脚本。
     */
    public void mkJavaCmdchr1_42() {
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/04_javaCMD/fastCall_chr1_42.sh";
        //nohup java -Xms200g -Xmx500g -jar FastCall.jar parameters_001_FastCall.txt > log_001_fastcall.txt &
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < 45; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                bw.write("java -Xms200g -Xmx500g -jar FastCall.jar parameters_");
                bw.write(chr);
                bw.write("_FastCall.txt > log_");
                bw.write(chr);
                bw.write("_fastcall.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 本方法的目的是：建立42条染色体的parameters文件。
     */
    public void mkParameterchr1_42() {
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/003_fastCall/000_uploadParameters/003_Parameters";
        try {
            for (int i = 0; i < 45; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_FastCall.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("Author: Fei Lu\n");
                sb.append("Email: flu@genetics.ac.cn; dr.lufei@gmail.com\n");
                sb.append("Homepage: https://plantgeneticslab.weebly.com/\n").append("\n");
                sb.append("#This WGS SNP discovery pipeline are designed for both heterozygous and inbred species, especially the depth is high (e.g. >10X each taxon).\n");
                sb.append("#To run and pipeline, the machine should have Java 8 and samtools installed. The lib directory should stay with FastCall.jar in the same folder. Command (e.g.): java -Xmx200g -jar FastCall.jar parameter.txt > log.txt &\n");
                sb.append("#To specify reference, bam direcotory, chromosome, and output direcotory, please edit the the parameters below. Also, please keep the order of parameters.\n").append("\n").append("\n");
                sb.append("#Parameter1:\tReference genome file with an index file (.fai). The reference should be in FastA format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).\n");
                sb.append("/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa").append("\n").append("\n");
                sb.append("#Parameter2:\tTaxa bam information file, including the info about what bams are included for each taxon\n");
                sb.append("/data4/home/aoyue/vmap2/abd/taxaBamMap.txt").append("\n").append("\n");
                sb.append("#Parameter3:\tChromosome or region on which genotyping will be performed (e.g. chromosome 1 is designated as 1. Region 1bp to 100000bp on chromosome 1 is 1:1,100000)\n");
                sb.append(i).append("\n").append("\n");
                sb.append("#Parameter4:\tVCF output directory\n");
                sb.append("/data4/home/aoyue/vmap2/abd/rawVCF/").append("\n").append("\n");
                sb.append("#Parameter5:\tNumber of threads for pileup\n");
                sb.append("32");
                bw.write(sb.toString());bw.newLine();
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    /**
     * 对抽样的数据进行质控，并试图将质控结果进行合并。
     */
    public void fastQC(){
        String inputDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/003_sampleFastq";
        String outputDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/004_fastQC";
        try {
            StringBuilder sb = new StringBuilder("/Users/Aoyue/software/FastQC/fastqc");
            File[] fs = new File (inputDirS).listFiles();
            fs =  IOUtils.listFilesEndsWith(fs, ".gz");
            for (int i = 0; i < fs.length; i++) {
                sb.append(" ").append(fs[i].getAbsoluteFile());
            }
            sb.append(" -o ").append(outputDirS);
            String cmd = sb.toString();
            System.out.println(cmd);
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            System.out.println("Fastqc evalutation is finished at " + outputDirS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
        
    }
    
    /**
     * 随机抽查1个样品，包含16个1端和16个2端，每个样品10万条reads，进行质控。
     */
    
    public void sampleFastQC(){
        String infileDirS = "/Volumes/LuLab4T_03/CleanData_15Samples/BT01374/";
        String outputDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/003_sampleFastq";
        int readNum = 100000;
        int startPoint = 1000000;
        
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        fs = IOUtils.listFilesEndsWith(fs, "_1.clean.fq.gz");
        
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_1.clean.")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            
            String infileDir1 = new File (infileDirS, name).getAbsolutePath();
            String infile1 = new File (infileDir1, name + "_1.clean.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDir1, name+"_2.clean.fq.gz").getAbsolutePath();
            
            String outfile1 = new File (outputDirS, name+"_1.clean.sample.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.clean.sample.fq.gz").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                String temp = null;
                int cnt = 0;
                while ((temp = br1.readLine()) != null) {
                    cnt++;
                    if (cnt < startPoint) {
                        br1.readLine();br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                    else {
                        bw1.write(temp+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                        bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                        for (int i = 0; i < readNum-1; i++) {
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        br1.close();
                        br2.close();
                        break;
                    }
                }
                System.out.println(name+ " completed");
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            
        });
    }
    
    /**
     * 随机挑选一个样品，进行fastq文件的reads数统计，并计算符合 10X 样品的size大小，从而确定数据量不够的样品。
     */
    
    public void calFastqbp(){
        //String infileDirS = "/Volumes/LuLab4T_03/CleanData_15Samples/BT01411";
        String infileDirS = "/Volumes/LuLab4T_03/CleanData_15Samples/BT01407";
      
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/002_countReads/";
        File f = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(f);
        fs = IOUtils.listFilesEndsWith(fs, "_1.clean.fq.gz");
        int cntReads = 0;
        /************************* Method1:单线程运行 *****************************************/
//        try{
//            for(int i=0; i<fs.length; i++){
//                String infileS = fs[i].getAbsolutePath();
//                BufferedReader br = IOUtils.getTextGzipReader(infileS);
//                String temp = null;
//                while((temp=br.readLine()) != null){
//                    cntReads++;
//                    br.readLine();br.readLine();br.readLine();
//                }
//                br.close();
//                cntReads += cntReads;
//            }
//            System.out.println( f.getName() + " totalreads count is  " + cntReads);
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }

        /************************* Methods2:多线程运行 *****************************************/
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(p -> {
            try{
                String infileS = p.getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(infileS);
                String temp = null;
                int cnt = 0;
                while((temp=br.readLine()) != null){
                    cnt++;
                    br.readLine();br.readLine();br.readLine();
                }
                br.close();
                System.out.println(p.getName() + "  reads count is  " + cnt);
            }
            catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            }
        });
        
        int totalReads = 1666383652;
        
        int wheatGenome = 145000000; //这里将bp替换成了reads,除以100
        
        
        /**
         * 规律，双端reads的数目是一致的，故只需要计算一端的reads即可。
         * 成功构建 (总时间: 24 分钟 54 秒)
         */
    }
    
    /**
     * 将6个盘的数据进行合并，并作为一个库DB，让ori文件中的每个样品作为query，进行搜索判断。
     * 若全部搜到，说明盘里数据齐全；若没有搜到，说明数据缺失。
     * 对库文件进行重复性验证，结果BT01505有重复。
     */
    public void mergeIDlist(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/001_mergeID/000_mergeIDlist.txt";
        String oriFileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/origin_wheat200IDlist.txt";
        File[] fs = new File(infileDirS).listFiles();
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab8T_1.IDlist.txt");
            bw.write(br.readLine() + "\n");
            br.close();
            //bw.write("ID_BGI\tSequenceCount\tSize(G)\n");
            for(int i =0; i<fs.length; i++){
                br = IOUtils.getTextReader(fs[i].getAbsolutePath());
                String temp = null;
                while((temp = br.readLine()) != null){
                    if(temp.startsWith("ID_BGI")) continue;
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(oriFileS);
            RowTable<String> t = new RowTable<>(outfileS);
            List<String> l = t.getColumn(0);
            Collections.sort(l);
            int cnt =0;
            String temp = null;
            while((temp = br.readLine()) != null){
                String query = PStringUtils.fastSplit(temp).get(0);
                if(query.startsWith("ID"))continue;
                int index = Collections.binarySearch(l, query);
                if(index <0) {
                    System.out.println("NONE    " + query);
                }
                else{
                    cnt ++;
                }
            }
            br.close();
            System.out.println(cnt);
            /******************* 查找重复 ********************************/
            HashMap<String, Integer> hashMap = new HashMap<String, Integer>();
            for (String string : l) {
                if (hashMap.get(string) != null) {
                    Integer value = hashMap.get(string);
                    hashMap.put(string, value+1);
                    System.out.println("the element: "+string+" is repeat");
                } else {
                    hashMap.put(string, 1);
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        //NONE    BT01504 199
    }
    
    public static ArrayList<String> getDirectory (String path) {
        ArrayList<String> files = new ArrayList<String>();
        File file = new File(path);
        File[] tempList = file.listFiles();

        for (int i = 0; i < tempList.length; i++) {
            if (tempList[i].isFile()) {
                //files.add(tempList[i].toString());
            }
            if (tempList[i].isDirectory()) {
                files.add(tempList[i].toString());
            }
        }
        return files;
    }

    /**
     * 将硬盘中的数据进行一一统计，检查200样品是否全部寄到。
     */
    public void GetIDlist() {
        //String infileDirs = "/Volumes/LuLab8T_1";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab8T_1.IDlist.txt";
        
       //String infileDirs = "/Volumes/LuLab8T_2";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab8T_2.IDlist.txt";
        
        //String infileDirs = "/Volumes/LuLab8T_3";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab8T_3.IDlist.txt";

        //String infileDirs = "/Volumes/LuLab4T_01/";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab4T_01.IDlist.txt";
        
        //String infileDirs = "/Volumes/LuLab4T_02/";
        //String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab4T_02.IDlist.txt";
        
        String infileDirs = "/Volumes/LuLab4T_03/CleanData_15Samples";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/000_cleandata/000_dataCheck/000_LuLab4T_03.IDlist.txt";
        
        File f = new File(infileDirs);
        String[] x = f.list();
        int cnt = 0;
        /****** 获取此路径的文件夹名字 ********/
        String storagename = f.getName();
        System.out.println(storagename);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("ID_BGI\tStorage\tSequenceCount\tSize_G" + "\n");
            for (int i = 0; i < x.length; i++) {
                if (!x[i].startsWith("BT")) {
                    continue;
                }
                
                cnt++; //对taxa进行计数
                //System.out.println(x[i]);
                //******************************开始对这一个taxa进行操作****************************************
                File subf = new File(infileDirs,x[i]);
                String subfpath = subf.getAbsolutePath();
                int subfDirectory = getDirectory(subfpath).size(); //得到每个样本的测序次数
                File[] fs = IOUtils.listRecursiveFiles(subf);
                double filesize = 0;
                for(int j = 0; j < fs.length; j++){
                    if(fs[j].getName().endsWith(".fq.gz")){
                        filesize = filesize + fs[j].length()/1000/1000/1000;
                    }
                }
                bw.write(x[i] + "\t" + storagename + "\t" + String.valueOf(subfDirectory) + "\t" + String.format("%.0f", filesize)); bw.newLine();
            }
            bw.flush();bw.close();
            System.out.println(infileDirs + " contains " + cnt + " taxa");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
