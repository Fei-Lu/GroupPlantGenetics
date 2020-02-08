/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

/**
 *
 * @author kanglipeng
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */


import lipengKang.analysis.DataProcessor;
import lipengKang.analysis.PlateAndID;
import com.google.common.primitives.Bytes;
import pgl.infra.table.RowTable;
import pgl.infra.table.TableInterface;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
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
public class WheatBWA {
    public WheatBWA(){
//        this.mkMd5();
//        this.checkMd5();
//        this.sampleData(); //单线程单个抽样
//        this.sampleFastq(); //多线程多个抽样
//        this.fastQC();
//        this.alignBWA();
//        this.listSpecificalFiles();
//        this.samtoolsSort();
//        this.samtoolsMerge();
        this.pipeline();
//        this.statReadLength();
//        this.statScore();
        //this.sampleBamfile();
        
        
    }
    
    /*
    本方法的目的是，对一个文件夹内的bam文件进行抽样,多线程运行。
    1.给出一个table,按照t里的taxa名字进行抽样；
    2.给出一个文件夹，把文件夹内所有的bam文件进行抽样。
    3.给出一个String taxa字符，对一个文件进行抽样。
    */
    public void sampleBamfile(){
        
        String outfileDirs = "/data2/aoyue/fastCall_project/bam/";
        String infileDirS = "/data2/sharedData/wheat_Jiao/splitBamfile/001/";
        String scriptDirS = "/data2/aoyue/fastCall_project/";
        
        /********************** 根据表格进行抽样 **********************************/
        String tfileS = "/data2/aoyue/fastCall_project/taxaFastCalltest.t.txt";
        RowTable<String> t = new RowTable<>(tfileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        
        /********************** 根据文件夹进行抽样 **********************************/
//        File[] fs = new File(infileDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".bam");
//        List<String> namelist = new ArrayList<>();
//        for(int i=0; i<fs.length; i++){
//            String name = fs[i].getName().split(".chr")[0];
//            namelist.add(name);
//        }
//        Collections.sort(namelist);
        
        /********************** 根据字符串进行抽样 **********************************/
//        String name = "AFG-L1";
//        List<String> namelist = new ArrayList<>();
//        namelist.add(name);
        
        /********************** 抽样1：run java包直接进行抽样 **********************************/
        namelist.parallelStream().forEach(taxa ->{
            //写命令脚本
            String infileS = new File(infileDirS, taxa + ".chr001.bam").getAbsolutePath();
            String outfileS = new File (outfileDirs,taxa + ".chr001.test.bam").getAbsolutePath();
            StringBuilder sb = new StringBuilder("samtools view -h ");
            sb.append(infileS).append(" 1:10000-150000 -o ").append(outfileS).append(" \nsamtools index ").append(outfileS);
            String cmd = sb.toString();
            try{
                //写脚本到文件中--将命令写入指定的脚本文件中
                String scriptS = new File(scriptDirS,taxa + ".sh").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                bw.write(cmd);
                bw.newLine();
                bw.flush();bw.close();
                //调用系统命令 - 新建一个缓冲，用于运行脚本
                sb = new StringBuilder("sh ");
                sb.append(scriptS);
                Runtime run = Runtime.getRuntime();
                Process p = run.exec(sb.toString());
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                //BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                StringBuilder ssb = new StringBuilder(); //有时候需要将屏幕信息输出到一个文件中，此时建立输出文本路径。
                String line = null;
                while ((line = br.readLine()) != null) {
                    System.out.println(line);
                    ssb.append(line + "\n");
                }
                String result = ssb.toString(); //有必要时输出出来
                p.waitFor();
                System.out.println(sb.toString());
                System.out.println("Sample Bamfile " + taxa + ".chr001.bam is finished.");
                new File (scriptS).delete();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        });
        
        
        /********************** 抽样2：所有信息写到一个脚本，上传脚本文件，执行 **********************************/
//        for(int i =0; i < namelist.size(); i++){ //一共有60个循环
//            String bamName = namelist.get(i);
//            //String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/splitbam/" + bamName + "_split.sh";
//            String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/splitbam/" + bamName + "_split.sh";
//            //samtools view -h mergeWheat24SM.bam 44 -o mergeWheat24SM.chr44.bam
//            try{
//                String inputDirS = "/data2/aoyue/output/bamfile/";
//                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
//                for(int j=0; j<45;j++){
//                    String chr = PStringUtils.getNDigitNumber(3, j);
//                    String outputDirS = "/data2/aoyue/output/splitBamfile/" + chr + "/";
//                    bw.write("samtools view -h " + inputDirS + bamName + ".rmdup.bam " + j +" -o " + outputDirS + bamName + ".chr" + chr +".bam && samtools index " + 
//                            outputDirS + bamName + ".chr" + chr + ".bam");
//                    bw.newLine();
//                }
//                bw.flush();bw.close();
//            }
//            catch(Exception e){
//                e.printStackTrace();
//                System.exit(1);
//            }
//        }
        
    }
    

    public void statScore () {
        String infileDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/001_sampleFq/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_1")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            try{
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                String temp = null;
                int allreads =0;
                int unqualifiedReads =0;
                while ((temp = br1.readLine()) != null) {
                    allreads++;
                    String line1 = br2.readLine();
                    String reads1 = br1.readLine(); String reads2 = br2.readLine();
                    String jia1 = br1.readLine(); String jia2 = br2.readLine(); 
                    String quality1 = br1.readLine(); String quality2 = br2.readLine();
                     
                    /****对1端的reads进行统计，将质量值转化为字符数组，进行ASCII转换，再排序，取前10个进行比较*****/
                    char[] chars = quality1.toCharArray(); 
                    byte[] b = new byte[chars.length];
                    for(int i=0;i <chars.length; i++){
                        b[i] = (byte)chars[i] ;
                        //System.out.println(b[i]);
                    }
                    Arrays.sort(b);
                    byte[] byte10 = Arrays.copyOfRange(b, 0, 10);
                    
                    int cal =0; //每次进行内部定义，判断这10个值的大小并统计；
                    for(byte i : byte10){
//                        if(i < 36) { //Phred值 3 + 33 
//                            cal++;
//                        }
                        if(i < 43) { //Phred值 10 + 33 
                            cal++;
                        }
//                        if(i < 48) { //Phred值 15 + 33 
//                            cal++;
//                        }
//                        if(i < 53) { //Phred值 20 + 33 
//                            cal++;
//                        }
//                        if(i < 55) { //Phred值 22 + 33 
//                            cal++;
//                        }
//                        if(i < 63) { //Phred值 30 + 33 
//                            cal++;
//                        }
                        //System.out.println(i);
                    }
                    if(cal > 9){
                        unqualifiedReads++;
                        
                        //System.out.println(cal);
                        //System.out.println(name + "_1.score.......... this read 不合格");
                    }
                }
                double ratio = unqualifiedReads/allreads;
                System.out.println(unqualifiedReads + "\t" + name.split("-")[1]);
                //System.out.println(ratio + "\t小于20" + name + "_1.score");
               br1.close();
               br2.close();   
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    
    public void statReadLength () {
        String infileDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/001_sampleFq/";
        String outputDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/003_readsLength/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_1")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.readLength.txt").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.readLength.txt").getAbsolutePath();
            try{
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextWriter(outfile2);
                String temp = null;
                int cnt =0;
                bw1.write("ID\tReadsLength");bw1.newLine();
                bw2.write("ID\tReadsLength");bw2.newLine();
                while ((temp = br1.readLine()) != null) {
                    cnt ++;
                    String line1 = br2.readLine();
                    String reads1 = br1.readLine(); String reads2 = br2.readLine();
                    String jia1 = br1.readLine(); String jia2 = br2.readLine(); 
                    String quality1 = br1.readLine(); String quality2 = br2.readLine();
                    int length1 = reads1.length(); 
                    int length2 = reads2.length();
                    bw1.write(cnt + "\t" + length1); bw1.newLine();
                    bw2.write(cnt + "\t" + length2); bw2.newLine();
                }
                bw1.flush();bw1.close();br1.close();
                bw2.flush();bw2.close();br2.close();   
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    
    
    
    
    /**
     * 本程序适合于一个样品的深测序质控；每个脚本单独运行完毕后，开始运行下一个脚本。
     * 输出脚本目录：ScriptParentS
     * 数据输出父目录：HPCS
     * Name列表：infileS
     * 原始数据文件夹：rawdataDirS
     */
    
    public void pipeline(){
        
        String parentFileS = "/data1/home/lipeng/output20190227wheat-pcr-free/output/"; //在本地主机上 建立工作文件夹
        String HPCS = "/data1/home/lipeng/result/output20190227wheat-pcr-free/output/"; //集群工作父目录
        String ScriptParentS = "/data1/home/lipeng/output20190227wheat-pcr-free/script/"; //脚本生成文件 父目录
        String bamS = "bamfile/";  //子文件夹
        String mergeS = "mergefile/"; //子文件夹
        String statS = "statsfile/"; //子文件夹
        String testS = "testbamfile/"; //子文件夹
        
        /*********************** 建立工作区(在本地路径中) **************************/
        new File(parentFileS).mkdirs();
        new File(parentFileS,bamS).mkdirs();
        new File(parentFileS,mergeS).mkdirs();
        new File(parentFileS,statS).mkdirs();
        new File(parentFileS,testS).mkdirs();
        
        /*********************** 建立工作区(在HPC中) **************************/
        String mkdirScriptS = ScriptParentS +"000_mkdir.sh";
        try{
            BufferedWriter bw = IOUtils.getTextWriter(mkdirScriptS);
            bw.write("mkdir " + HPCS);bw.newLine();
            bw.write("mkdir " + HPCS + bamS);bw.newLine();
            bw.write("mkdir " + HPCS + mergeS);bw.newLine();
            bw.write("mkdir " + HPCS + statS);bw.newLine();
            bw.write("mkdir " + HPCS + testS);bw.newLine();
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        String infileS = "/data1/home/lipeng/20190227wheat-pcr-free.txt"; // 建立Name列表
        //String infileS = "/Users/Aoyue/Documents/test2M.txt"; //测试数据
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        String[] names = namelist.toArray(new String[namelist.size()]);
        Arrays.sort(names);
        
        
        String bwaThread = "8";  //根据实际情况修改，使用的线程数
        String indexFileS = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz"; //比对参考基因组index library
        //String rawdataDirS = "/data2/aoyue/a-wheatRawdata1219";
        String rawdataDirS = "/data1/home/lipeng/database/clean/"; //原始数据文件夹
        String bamfileDirS = HPCS + bamS; 
        String bwaScriptS = ScriptParentS +"001_bwa.sh";
        
        
        String sortScriptS = ScriptParentS + "002_sort.sh";
        String sortMemory = "10G"; //根据实际情况修改，每个线程使用的最大内存
        String sortThread = "10"; //根据实际情况修改，使用的线程数
        
        String mergefileDirS = HPCS + mergeS;
        String mergeName = "mergeWheat" + String.valueOf(names.length) + "SM.bam";
        String mergeScriptS = ScriptParentS + "003_merge.sh";
        
        String RefSeqS = "/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa"; //统计stat需要解压的ref fa文件
        String bamQCScriptS = ScriptParentS + "004_bamQC.sh";
        String statfileDirS = HPCS + statS;
        
        
        
        
        
        /**************************** bwa **************************************/
        try {
            BufferedWriter bw = IOUtils.getTextWriter(bwaScriptS);
            for (int i = 0; i < names.length; i++) {
            //bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
                StringBuilder sb = new StringBuilder("bwa mem -t ");
                sb.append(bwaThread).append(" -R ").append("'@RG\\tID:").append(names[i]).append("\\t").append("PL:PACBIO").append("\\t").append("SM:").append("w1").append("\\t").append("LB:").append(names[i]).append("' ");
                sb.append(indexFileS).append(" ");
                sb.append(new File(rawdataDirS, names[i]+".fastq.gz").getAbsolutePath()).append(" ");
               // sb.append(new File(rawdataDirS, names[i]+"_2.fq.gz").getAbsolutePath());
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(bamfileDirS, names[i]+".pe.bam").getAbsolutePath()).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            
            BufferedReader br = IOUtils.getTextReader(bwaScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        /**************************** sort **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(sortScriptS);
            for (int i = 0; i < names.length; i++) {
       //samtools sort -m 1G -o /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.pos.bam -O bam -@ 10 /data2/aoyue/output/testbamfile/HRV-L1.pe.fixmate.bam &
                StringBuilder sb = new StringBuilder("samtools sort -m ");
                sb.append(sortMemory).append(" -o ").append(bamfileDirS).append(names[i]).append(".sorted.bam").append(" -O bam -@ ").append(sortThread).append(" ").append(bamfileDirS).append(names[i]).append(".pe.bam && ");
                sb.append("samtools index ").append(bamfileDirS).append(names[i]).append(".sorted.bam &");
                bw.write(sb.toString());bw.newLine();
            }
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(sortScriptS);
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        /**************************** merge **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(mergeScriptS);
            
            //samtools merge out.bam in1.bam in2.bam in3.bam
            StringBuilder sb = new StringBuilder("samtools merge ");
            sb.append(bamfileDirS).append(mergeName).append(" ");
            for (int i = 0; i < names.length; i++) {
                sb.append(bamfileDirS).append(names[i]).append(".sorted.bam").append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            sb.append(" && mv ").append(bamfileDirS).append(mergeName).append(" ").append(mergefileDirS);
            sb.append(" && ");
            sb.append("samtools index ").append(mergefileDirS).append(mergeName);
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(mergeScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        /**************************** bam QC **************************************/
        try{
            BufferedWriter bw = IOUtils.getTextWriter(bamQCScriptS);
            //samtools stats -r /data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa 180803.pe.sorted.bam > 180803.pe.sorted.bam.bc
            StringBuilder sb = new StringBuilder("samtools stats -r ");
            sb.append(RefSeqS).append(" ").append(mergefileDirS).append(mergeName).append(" > ");
            sb.append(mergefileDirS).append(mergeName).append(".bc");
            sb.append(" && ").append("mv ").append(mergefileDirS).append(mergeName).append(".bc ").append(statfileDirS);
            sb.append(" && ");
            sb.append("perl /data1/programs/samtools-1.8/misc/plot-bamstats -p ").append(statfileDirS).append(mergeName).append(".bc ")
                    .append(statfileDirS).append(mergeName).append(".bc");
            bw.write(sb.toString());
            bw.newLine();
            
            sb = new StringBuilder();
            sb.append("samtools flagstat ").append(mergefileDirS).append(mergeName).append(" > ").append(mergefileDirS).append(mergeName).append(".flagstat.txt &");
            bw.write(sb.toString());
            bw.newLine();
            
            
            bw.flush();bw.close();
            
            BufferedReader br = IOUtils.getTextReader(bamQCScriptS);
            System.out.println(br.readLine());
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    public void samtoolsMerge(){
        //String inputDirS = "/Users/Aoyue/Downloads/huadagene/bamfile";
        //String outputPerlS = "/Users/Aoyue/Downloads/huadagene/runMerge.sh";
        String inputDirS = "/data1/home/aoyue/huada/bamfile";
        String outputPerlS = "/data1/home/aoyue/huada/runMerge.sh";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".bam");
        Arrays.sort(fs);
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            bw.write("samtools merge mergeWheat24SM.bam ");
            StringBuilder sb = new StringBuilder();
            for(int i = 0; i < fs.length; i++){
                String inBam = fs[i].getAbsolutePath();
                System.out.println(inBam);
                sb.append(inBam).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.write(" &");
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void samtoolsSort(){
        //samtools sort -m 6G -o 180803.pe.sorted.bam -O bam -@ 8 180803.pe.bam && echo "** BAM sort done **"
        String inputDirS = "/data1/home/aoyue/huada/bamfile/";
        String outputPerlS = "/data1/home/aoyue/huada/runSort.pl";
        File[] fs = new File (inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".bam");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            String fullName = fs[i].getName();
            String[] all = fullName.split(".pe");
            nameSet.add(all[0]);
        }
        String[] names = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(names);
        int memory = 8;
        int numthreads = 8;
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outputPerlS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder("samtools sort");
                sb.append(" -m ").append(memory).append("G").append(" -o ");
                sb.append(new File(inputDirS,names[i] + ".pe.sorted.bam" ).getAbsolutePath());
                sb.append(" -O ").append("bam").append(" -@ ").append(numthreads).append(" ");
                sb.append(new File (inputDirS,names[i] + ".pe.bam" ).getAbsolutePath());
                sb.append(" && ");
                sb.append("echo \"** BAM sort done **\"");
                sb.append(" && ");
                sb.append("samtools index ");
                sb.append(new File(inputDirS,names[i] + ".pe.sorted.bam" ).getAbsolutePath());
                sb.append(" && ");
                sb.append("echo \"** BAM index done **\"").append(" &");
                String cmd = sb.toString();
                //String command = "system(\""+cmd+"\");";
                String command = cmd;
                bw.write(command);
                bw.newLine();
            }
            bw.flush();
            bw.close();
//            StringBuilder sb = new StringBuilder("perl ");
//            sb.append(outputPerlS);
//            Runtime run = Runtime.getRuntime();
//            Process p = run.exec(sb.toString());
//            System.out.println(sb.toString());
//            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//            p.waitFor();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void testsplit(){
        String inDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "1.clean.fq.gz");
        String outfileS = "/Users/Aoyue/Documents/list.sh";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fs.length; i++) {
                sb = new StringBuilder();
                String fullName = fs[i].getName();
                String[] all = fullName.split("_",4);
                String sample = all[3];
                sb.append(fs[i].getName());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(fs.length);  
    }
    
    private void listSpecificalFiles() {
        String inDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        File[] fs = new File(inDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        String outfileS = "/Users/Aoyue/Documents/Datalist.txt";
        //String[] header = {"FileName", "FileSize"};
        String[] header = {"FileName"};
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < header.length; i++) {
                sb.append(header[i]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < fs.length; i++) {
                sb = new StringBuilder();
                
                double qq = (double)fs[i].length()/1024/1024/1024;
                DecimalFormat df = new DecimalFormat ("0.00");

                sb.append(fs[i].getName());
                        //.append("\t").append(df.format(qq));
//                sb.append(subFs[i].getName()).append("\t").append(subFs[i].getAbsolutePath()).append("\t").append((double)subFs[i].length()/1024/1024/1024).append("\t").append("350bp");
                
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(fs.length);    
    }
    
    private void alignBWA () {
        String infileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/001_script/SM8.t.txt";
        String thread = "16";  //根据实际情况修改，使用的线程数
        String indexFileS = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz";
        String inputDirS = "/data2/aoyue/a-wheatRawdata1219";
        String outputDirS = "/data2/aoyue/a-output/bamfile/";
        String ScriptS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/001_script/001_bwa.sh";
        
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        String[] names = namelist.toArray(new String[namelist.size()]);
        
        //bwa mem -t 100 -R '@RG\tID:HRV-L1\tPL:illumina\tSM:HRV-L1\tLB:HRV-L1' /data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz /data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz /data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz | samtools view -S -b - > /data2/aoyue/output/bamfile/HM.pe.bam && echo "** bwa mapping done **" &
        try {
            BufferedWriter bw = IOUtils.getTextWriter(ScriptS);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder("bwa mem -t ");
                sb.append(thread).append(" -R ").append("'@RG\\tID:").append(names[i]).append("\\t").append("PL:illumina").append("\\t").append("SM:").append(names[i]).append("\\t").append("LB:").append(names[i]).append("' ");
                sb.append(indexFileS).append(" ");
                sb.append(new File(inputDirS, names[i]+"_1.fq.gz").getAbsolutePath()).append(" ");
                sb.append(new File(inputDirS, names[i]+"_2.fq.gz").getAbsolutePath());
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(outputDirS, names[i]+".pe.bam").getAbsolutePath()).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(ScriptS);
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    private void fastQC () {
//        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/";
//        String outputDirS = "/Users/Aoyue/Downloads/huadagene/fastqc/";
        String inputDirS = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A/";
        String outputDirS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/000_fastqc";
        
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
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getErrorStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            System.out.println("Fastqc evalutation is finished at" + outputDirS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }
    }
    
    public void sampleFastq () {
        String infileDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A/";
        String outputDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/001_sampleFq/";
        String outputFastaDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/SampleOne20181224/data/002_sampleFasta/";
        int readNum = 100000; //抽取10万条reads
        int startPoint = 100000; //从第10万条后开始抽样
        int fastaNum = 1000; //抽取1000条fasta序列
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_1")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfile1 = new File (outputDirS, name+"_1.fq.gz").getAbsolutePath();
            String outfile2 = new File (outputDirS, name+"_2.fq.gz").getAbsolutePath();
            String outfileFasta = new File (outputFastaDirS, name+"_1.fa").getAbsolutePath();
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1);
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2);
                BufferedWriter bwf = IOUtils.getTextGzipWriter(outfileFasta);
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
                            bw1.write(br1.readLine()+"\n");
                            temp = br1.readLine();bw1.write(temp+"\n");
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                            if (i > fastaNum) continue;
                            bwf.write(">"+String.valueOf(i));
                            bwf.newLine();
                            bwf.write(temp);
                            bwf.newLine();
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        bwf.flush();bwf.close();
                        br1.close();
                        br2.close();
                        break;
                    }
                }
                System.out.println(String.valueOf(readNum) + " reads are sampled from"+ name);
            }
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
    
    /*
    本方法适用于单个样品的抽样，不适合多线程抽样。
    */
        public void sampleData(){
        String infile1S = "/data2/sharedData/Jiao/ABD/HRV-L1_1.fq.gz";
        String infile2S = "/data2/sharedData/Jiao/ABD/HRV-L1_2.fq.gz";
        String outfile1S = "/data2/aoyue/test/HRV-L1_test_1.fq.gz";
        String outfile2S = "/data2/aoyue/test/HRV-L1_test_2.fq.gz";
        String fastaS = "/data2/aoyue/HRV-L1_fasta_1.fa";
        
        int readNum = 100000;
        int startPoint = 100000;
        int fastaNum = 1000;
        try{
            BufferedReader br1 = IOUtils.getTextGzipReader(infile1S);
            BufferedReader br2 = IOUtils.getTextGzipReader(infile2S);
            BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfile1S);
            BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfile2S);
            BufferedWriter bwf = IOUtils.getTextGzipWriter(fastaS);
            String temp = null;
            int cnt = 0;
            while((temp = br1.readLine()) != null){
                cnt++;
                if(cnt < startPoint){
                    br1.readLine();br1.readLine();br1.readLine();
                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                }
                else{
                    bw1.write(temp+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                    bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                    for (int i = 0; i < readNum-1; i++) {
                            bw1.write(br1.readLine()+"\n");
                            temp = br1.readLine();bw1.write(temp+"\n");
                            bw1.write(br1.readLine()+"\n");bw1.write(br1.readLine()+"\n");
                            bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");bw2.write(br2.readLine()+"\n");
                            if (i > fastaNum) continue;
                            bwf.write(">"+String.valueOf(i));
                            bwf.newLine();
                            bwf.write(temp);
                            bwf.newLine();
                        }
                        bw1.flush();bw1.close();
                        bw2.flush();bw2.close();
                        bwf.flush();bwf.close();
                        br1.close();
                        br2.close();
                        break;
                }
            }
            System.out.println(String.valueOf(readNum) + " reads are sampled from"+ "HRV-L1");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("。。。。。。。");
    }
    
    private void checkMd5() {
        //ori 文件 54de998c7883a4592a1683bec2590d64  K16HL0119_1_clean.fq.gz
        //ori 文件 0dc05f49a968937197f20841ef8427b2  Clean/WHYD18111796_A/V300010959_L1_B5RDWHEsikRAAAAA-533_1.fq.gz
        //des文件 MD5 (/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/180803_I13_V100004234_L1_WHEkapRAADT-585_1.clean.fq.gz) = 472e1633b40b0ae3866a820f5b8637a5
        //des文件 MD5 (/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A/V300010959_L1_B5RDWHEsikRAAAAA-533_1.fq.gz) = 0dc05f49a968937197f20841ef8427b2
        String ori = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/md5.txt"; //原始md5
        String des = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/Wheat20181219.md5.txt"; //mac生成的md5文件
        HashMap<String, String> fmd5Map = new HashMap<>(); //建立一个键值对应的hashmap,此时hashmap为空。下文会把原始的ori文件放入hashmap中去
        TableInterface oT = new RowTable(ori, "  "); //分隔符是2个空格，将ori文件读进表格
        TableInterface dT = new RowTable(des, " "); //分隔符是1个空格，将des文件读进表格
        String choice = "YES";
        try{
            BufferedReader br = IOUtils.getTextReader(ori); //先读进原始生成的md5文件
            String temp;
            List<String> l = null;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, "  ");
                //fmd5Map.put(l.get(1), l.get(0));
                fmd5Map.put(l.get(1).split("_A/")[1], l.get(0));
            }
            br.close();
            int cnt = 0;
            br = IOUtils.getTextReader(des); //读进mac生成md5文件
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, " ");
                //String key = l.get(1).split("00010496/")[1].replaceFirst("\\)", "");
                String key = l.get(1).split("_A/")[1].replaceFirst("\\)", "");
                String value = fmd5Map.get(key);
                /********* 文件全部下载完毕，可以进行一一对应 ***********/
                if (value == null) { //如果value为空，则为真。
                System.out.println(key+"\tdoesn't exist");
                continue;
                }
                if (value.equals(l.get(3))) {
                    cnt++;
                    continue;
                }
                System.out.println(key + "\t is incorrect");
            }
            System.out.println(cnt + " key is correct");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mkMd5(){
        /*该方法只用将文件输入路径和输出文件名进行修改，即可使用*/
//        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
//        String outfileS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/checkmd5.txt";
        //String inputDirS = "/Users/Aoyue/Documents/Data/project/WheatGBSVI/02IlluminaSeq/data/2.cleandata/20180601-51-NEB12_TKD180600155";
        //String inputDirS = "/Users/Aoyue/Documents/P18Z12200N0030_WHEsikR/Clean/WHYD18111796_A";
        String inputDirS = "/Volumes/LuLab3T_42/a-wheatRawdata1219";
        String outfileS = "/Users/Aoyue/Documents/output20181219wheat-pcr-free/Wheat20181219.md5.txt";
        String suffix = ".gz";
        
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            /*1#新建存放数据的文件对象fall，将该目录下所有以fq.gz结尾的文件名列出来，并存放在fs文件数组中*/
            File fall = new File (inputDirS);
            File[] fs = IOUtils.listRecursiveFiles(fall);
            fs =  IOUtils.listFilesEndsWith(fs, suffix);
            /*2# 根据fq.gz结尾的文件，生成md5计算的脚本，输入到StringBuilder sb中（这里用循环的方法）*/
            for (int i = 0; i < fs.length; i++) {
                sb.append("md5").append(" ").append(fs[i].getAbsoluteFile()).append("\n");
            }
            /*3# 将生成的sb脚本转换成String类型，对象名为cmd，并打印出来检查是否无误*/
            String cmd = sb.toString();
            System.out.println(cmd);
            /*4# 调用终端命令，运行md5命令，并将输出的结果写入输出文件outfileS中*/
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);           
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream())); 
            //p.getInputStream()返回InputStream类型；InputStreamReader构造器函数需要加 InputStream 类型的参数； 
            //InputStreamReader类继承了Reader类. BufferedReader构造器需要加Reader类的参数。
            StringBuilder ssb = new StringBuilder();
            String line;
            while ((line = br.readLine()) != null) {
            ssb.append(line).append("\n");
            }
            String result = ssb.toString();
            System.out.println(result);  
            p.waitFor();
            bw.write(result);
            bw.flush();
            bw.close();            
            System.out.println("md5Check is finished at " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);       
        }       
    }
    
    
    public static void main (String[] args){
        new WheatBWA();
        //new ProcessVCF();
        System.out.println("\n**********************************" );
        System.out.println("Here is the main class of WheatBWA" );
        
    }  
}

