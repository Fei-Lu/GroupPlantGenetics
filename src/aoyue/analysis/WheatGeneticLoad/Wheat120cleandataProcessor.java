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
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Wheat120cleandataProcessor {
    public Wheat120cleandataProcessor(){
        //this.getName();
        //this.sort();
        //this.checkDupli();
        //this.newCheckSecondData();
        //this.sampleData();
        this.splitBam();
//        this.splitBam_chr0_mit_chl();
//        this.calSingleChrSize();
        this.mergeCalFile();
        //this.mkIndexFile();
        //this.cp_First();
//        this.cpFirstChr();
//        this.cp_Second();
        this.mergeSplitbam();
        
        
    }
    
    public Wheat120cleandataProcessor(String inDirS, String outS, String suffix, String unit){
        new Statistic(inDirS, outS, suffix, unit);
    }
    
    public void mergeSplitbam(){
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt";
        String outfileDirS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/mv_secondSplitbam/";
        String oriPathDirS = "/data2/aoyue/splitBamfile_Second/";
        String desPathDirS = "/data2/sharedData/Jiao/splitBamfile/";
        
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        
        for(int i =0; i < namelist.size(); i++){
            String bamName = namelist.get(i);
            String scriptS = outfileDirS + bamName + "-mv_secondSplitbam.sh";
            try{
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                for(int j=0; j<45;j++){
                    String chr = PStringUtils.getNDigitNumber(3, j);
                    String oriPathS = oriPathDirS + chr + "/";
                    String desPathchrDirS = desPathDirS + chr + "/";
                    bw.write("mv " + oriPathS + bamName + ".chr" + chr + ".bam "+ desPathchrDirS + " && " 
                            + "mv " + oriPathS + bamName + ".chr" + chr + ".bam.bai "+ desPathchrDirS);
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
    
    public void cp_Second(){
        
    }
    
    public void cpFirstChr(){
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String outfileS = doc + "First-cp_split.bam_chr000-016.sh"; //写出的脚本文件路径 (需要改)
//        String outParentS = "/mnt/usb/";
//        String output = "ABD001_split.bam_chr000-016/";
//        String outputS = outParentS + output;//脚本中的命令，目标文件夹 （需要改）
        
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String outfileS = doc + "First-cp_split.bam_chr017-038.sh"; //写出的脚本文件路径 (需要改)
//        String outParentS = "/mnt/usb/";
//        String output = "ABD001_split.bam_chr017-038/";
//        String outputS = outParentS + output;//脚本中的命令，目标文件夹 （需要改）
        
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String outfileS = doc + "First-cp_split.bam_chr039-040.sh"; //写出的脚本文件路径 (需要改)
//        String outParentS = "/mnt/usb/";
//        String output = "ABD001_split.bam_chr039-040/";
//        String outputS = outParentS + output;//脚本中的命令，目标文件夹 （需要改）
//        
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String outfileS = doc + "First-cp_split.bam_chr041-044.sh"; //写出的脚本文件路径 (需要改)
        String outParentS = "/mnt/usb/";
        String output = "ABD001_split.bam_chr041-044/";
        String outputS = outParentS + output;//脚本中的命令，目标文件夹 （需要改）
        
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for(int j=0; j<45;j++){
                //if(j>16) continue;
                
//                if(j<17) continue;
//                if(j>38) continue;
                
//                if(j<39) continue;
//                if(j>40) continue;
                
                if(j<41) continue;
                if(j>44) continue;
                
                
                String chr = PStringUtils.getNDigitNumber(3, j);
                String infileDirS = "/data2/aoyue/output/splitBamfile/" + chr + "/";
                StringBuilder sb = new StringBuilder();
                sb.append("cp Rf ").append(infileDirS).append(" ").append(outputS);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(outfileS);
            System.out.println(br.readLine());
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
        
    
    
    public void cp_First(){
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "First-cp_fixmate.pos.bam_1-28.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
//        String output = "ABD001_fixmate.pos.bam_1-28/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名
        
        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
        String outfileS = doc + "First-cp_fixmate.pos.bam_29-55.sh"; //写出的脚本文件路径 (需要改)
        String infileDirS = "/data2/aoyue/output/bamsorted/"; //脚本中的命令，查找数据的地方 (需要改)
        String output = "ABD001_fixmate.pos.bam_29-55/"; //脚本中的命令，目标文件夹 （需要改）
        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
        String suffix = ".fixmate.pos.bam"; //脚本中的命令，list中的后缀名
        
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "First-cp_rmdup.bam_1-25.sh"; //写出的脚本文件路径
//        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
//        String output = "ABD001_rmdup.bam_1-25/"; //脚本中的命令，目标文件夹
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名
        
//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "First-cp_rmdup.bam_26-55.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
//        String output = "ABD001_rmdup.bam_26-55/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名

//        String doc = "/Users/Aoyue/Documents/"; //进入doc目录下操作
//        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //读取的list文件路径
//        String outfileS = doc + "First-cp_rmdup.bam_56-60.sh"; //写出的脚本文件路径 (需要改)
//        String infileDirS = "/data2/aoyue/output/bamfile/"; //脚本中的命令，查找数据的地方
//        String output = "ABD001_rmdup.bam_56-60/"; //脚本中的命令，目标文件夹 （需要改）
//        String outfileDirS = "/mnt/usb/" + output; //脚本中的命令，存放目标文件的文件夹路径
//        String suffix = ".rmdup.bam"; //脚本中的命令，list中的后缀名


        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = br.readLine(); //读入表头
            int cnt =0;
            while((temp = br.readLine()) != null){
                cnt++;
                
 //               if(cnt > 25) continue;
//                if(cnt < 26) continue;
//                if(cnt > 55) continue;
                //if(cnt > 28) continue;
                if(cnt <29) continue;
                if(cnt >55) continue;
                
                StringBuilder sb = new StringBuilder();
                sb.append("cp -Rf ").append(infileDirS).append(temp).append(suffix).append(" ").append(outfileDirS);
                sb.append(" && ");
                sb.append("cp -Rf ").append(infileDirS).append(temp).append(suffix).append(".bai").append(" ").append(outfileDirS);
                
                bw.write(sb.toString());bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            BufferedReader br = IOUtils.getTextReader(outfileS);
            System.out.println(br.readLine());
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mkIndexFile(){
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        for(int i =0; i < namelist.size(); i++){ //一共有60个循环
            String bamName = namelist.get(i);
            String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/index/" + bamName + "_index.sh";
            //samtools index mergeWheat24SM.bam
            try{
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                for(int j=0; j<45;j++){
                    String chr = PStringUtils.getNDigitNumber(3, j);
                    String outputDirS = "/data2/aoyue/output/splitBamfile/" + chr + "/";
                    bw.write("samtools index " + outputDirS + bamName + ".chr" + chr + ".bam");
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
    
    public void mergeCalFile(){
//        String inDirS = "/Users/Aoyue/Documents/calFirstSize/";  //没条染色体的文件存储
//        String OriginS = "/Users/Aoyue/Documents/calFirst.chrmerge.txt";  //文件第一列的取样名字
//        String outfileS = "/Users/Aoyue/Documents/calFirst.chr0-44.txt";  //合并后的文件夹
        
        String inDirS = "/Users/Aoyue/Documents/calSecondSize/";
        String OriginS = "/Users/Aoyue/Documents/calSecond.chrmerge.txt"; //到时候copy其中一个文件中的一个
        String outfileS = "/Users/Aoyue/Documents/calSecond.chr0-44.txt";
        
        RowTable<String> t = new RowTable<>(OriginS);
        List[] l = new List[45];
        for(int i = 0; i < 45;i++){
            l[i]=new ArrayList();
        }
        for(int i=0; i < 45; i++){
            String chr = PStringUtils.getNDigitNumber(3, i);
            //String infileS = inDirS + "calFirst.chr" + chr + ".bam.txt";
            String infileS = inDirS + "calSecond.chr" + chr + ".bam.txt";
            RowTable<String> tt = new RowTable<>(infileS);
            l[i] = tt.getColumn(1);
            if((i == 0) || (i == 43) || (i == 44)){
                t.addColumn(chr + "_Size(MB)", l[i]);
            }
            else{
                t.addColumn(chr + "_Size(GB)", l[i]);
            }
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    public void calSingleChrSize(){
//        String inDirS = "/data2/aoyue/output/spiltBamfile/";
//        String outDirS = "/data2/aoyue/calFirstSize/";
        
        String inDirS = "/data2/aoyue/output/splitBamfile/";
        String outDirS = "/data2/aoyue/calSecondSize/";
        try{
            BufferedWriter bw = IOUtils.getTextWriter("/Users/Aoyue/Documents/Java-jar-script.txt");
            for(int i = 0; i < 45; i++){
                //java -jar /data2/aoyue/calFileSize.jar /data2/aoyue/output/spiltBamfile/000 /data2/aoyue/calFirstSize/calFirst.chr000.bam.txt .chr000.bam MB
                String chr = PStringUtils.getNDigitNumber(3, i);
                String chrDirS = inDirS + chr;
                String suffix = ".chr" + chr + ".bam";
                //String outS = outDirS + "calFirst" + suffix + ".txt";
                String outS = outDirS + "calSecond" + suffix + ".txt";
                if((i == 0) || (i == 43) || (i == 44)){
                    bw.write("java -jar /data2/aoyue/calFileSize.jar ");
                    bw.write(chrDirS);bw.write(" ");bw.write(outS);bw.write(" ");bw.write(".chr");bw.write(chr);bw.write(".bam ");bw.write("MB");
                    bw.newLine();
                }
                else{
                    bw.write("java -jar /data2/aoyue/calFileSize.jar ");
                    bw.write(chrDirS);bw.write(" ");bw.write(outS);bw.write(" ");bw.write(".chr");bw.write(chr);bw.write(".bam ");bw.write("GB");
                    bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        try{
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/Documents/Java-jar-script.txt");
            System.out.println(br.readLine());
            System.out.println(br.readLine());
            br.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    public void splitBam_chr0_mit_chl(){
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //nameSet的路径
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        for(int i =0; i < namelist.size(); i++){ //一共有60个循环
            String bamName = namelist.get(i);
            String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/splitbam/" + bamName + "_spilt.sh";
            //samtools view -h mergeWheat24SM.bam 44 -o mergeWheat24SM.chr44.bam
            try{
                String inputDirS = "/data2/aoyue/output/bamfile/";
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                String chr0 = PStringUtils.getNDigitNumber(3, 0);
                String outputDirS = "/data2/aoyue/output/spiltBamfile/" + chr0 + "/";
                bw.write("samtools view -h " + inputDirS + bamName + ".rmdup.bam " + 0 +" -o " + outputDirS + bamName + ".chr" + chr0 +".bam");
                bw.newLine();
                String chr43 = PStringUtils.getNDigitNumber(3, 43);
                outputDirS = "/data2/aoyue/output/spiltBamfile/" + chr43 + "/";
                bw.write("samtools view -h " + inputDirS + bamName + ".rmdup.bam " + 43 +" -o " + outputDirS + bamName + ".chr" + chr43 +".bam");
                bw.newLine();
                String chr44 = PStringUtils.getNDigitNumber(3, 44);
                outputDirS = "/data2/aoyue/output/spiltBamfile/" + chr44 + "/";
                bw.write("samtools view -h " + inputDirS + bamName + ".rmdup.bam " + 44 +" -o " + outputDirS + bamName + ".chr" + chr44 +".bam");
                bw.newLine();
                bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
        }     
    }
    
    public void splitBam(){
        //String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/First-1-60SM.t.txt"; //nameSet的路径
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/Second-1-60SM.t.txt";
        //String infileS = "/Users/Aoyue/Documents/DKM.t.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        Collections.sort(namelist);
        for(int i =0; i < namelist.size(); i++){ //一共有60个循环
            String bamName = namelist.get(i);
            //String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/First1-60/splitbam/" + bamName + "_split.sh";
            String scriptS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/Jiao/002_script/Second1-60/splitbam/" + bamName + "_split.sh";
            //samtools view -h mergeWheat24SM.bam 44 -o mergeWheat24SM.chr44.bam
            try{
                String inputDirS = "/data2/aoyue/output/bamfile/";
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                for(int j=0; j<45;j++){
                    String chr = PStringUtils.getNDigitNumber(3, j);
                    String outputDirS = "/data2/aoyue/output/splitBamfile/" + chr + "/";
                    bw.write("samtools view -h " + inputDirS + bamName + ".rmdup.bam " + j +" -o " + outputDirS + bamName + ".chr" + chr +".bam && samtools index " + 
                            outputDirS + bamName + ".chr" + chr + ".bam");
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
 
        System.out.println("这是BWAJiao出口");
    }
    
    /**
     * 查找第二批数据是否和上次所得的 Jiao_nameList_NOdata.txt 信息一致，如果一致，说明数据正确无误。
     * 把 dbfileS 当成库，转化为数组，然后将infileS读进去，一一搜索。
     */
    public void newCheckSecondData(){
        String infileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60taxaname_later.txt";
        String dbfileS = "/Users/Aoyue/Documents/Jiao/001_db/Jiao_nameList_NOdata.txt";
        
        RowTable<String> t = new RowTable<>(dbfileS);
        List<String> l = t.getColumn(0);
        String[] db = l.toArray(new String[l.size()]);
        Arrays.sort(db);  /*千万不能忘了排序啊！！！*/
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            int cnt =0;
            while((temp = br.readLine()) != null){
                String query = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(db, query);
                if(index <0){
                    System.out.println(query);
                    cnt++;  
                }  
            }
            System.out.println(cnt);
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * 最后结果，一一对应，说明数据没有问题。
         */
    }
    
    public void checkDupli(){
        String infileS = "/Users/Aoyue/Documents/Jiao120taxa_country.txt";
        String Jiao60fileS = "/Users/Aoyue/Documents/Jiao60taxaname.txt";
        String outfileS = "/Users/Aoyue/Documents/Jiao_nameList_NOdata.txt";
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsText(0);
        List<String> l = t.getColumn(0);
        Set s = new HashSet(l);
        System.out.println(l.size());
        System.out.println(s.size());
        
        t = new RowTable<String>(Jiao60fileS);
        l = t.getColumn(0);
        String[] taxa60 = l.toArray(new String[l.size()]);
        Arrays.sort(taxa60);
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String query = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa60, query);
                if (index < 0) {
                    bw.write(temp);bw.newLine();
                }
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
        
        t = new RowTable<String>(outfileS);
        List<String> country = t.getColumn(1);
        Set<String> sco = new HashSet<>(country);
        for(String a : sco){
            System.out.println(a + "    " + Collections.frequency(country, a));
        }

    }
    
    public void sort(){
//        String infileS = "/Users/Aoyue/Documents/Jiao60name.txt";
//        String outfileS = "/Users/Aoyue/Documents/Jiao60taxaname.txt";
        
        String infileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60_later.txt";
        String outfileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60taxaname_later.txt";
        
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    public void getName(){
//        String infileDirS = "/Volumes/WheatAABBDD/ABD001";
//        String outfileS = "/Users/Aoyue/Documents/Jiao1.txt";
        
//        String infileDirS = "/Volumes/ABD002/ABD002";
//        String outfileS = "/Users/Aoyue/Documents/Jiao2.txt";
        
//        String infileDirS = "/Volumes/AABBDD003/ABD003";
//        String outfileS = "/Users/Aoyue/Documents/Jiao3.txt";
        
        String infileDirS = "/Volumes/Seagate Backup Plus Drive/NHT151096_60s/release_1";
        String outfileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60_later.txt";
        Set<String> s = new HashSet<>();
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        int cnt = 0;
        for(int i=0; i< fs.length;i++){
            String name = fs[i].getName().split("_")[0];
            System.out.println(name);
            s.add(name);
            cnt++;
        }
        System.out.println(cnt + "\tfq.gz 文件数");
        System.out.println(s.size() + "\tname文件数");
        
        String[] names = s.toArray(new String[s.size()]);
        Arrays.sort(names);
        
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Name");bw.newLine();
            for(int i=0; i< names.length;i++){
                bw.write(names[i]);bw.newLine();
            }
            bw.flush();bw.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 
     * @param p 
     */
    //public static void main (String[] p){
//        System.out.println("To use this calFileSize jar, please add 4 params\n@param1: infileDirS @param2: outfileS @param3: suffix @param4: unit, You can choose TB GB MB KB");
//        new Wheat120cleandataProcessor(p[0], p[1], p[2], p[3]);
        
    //new Wheat120cleandataProcessor();
    //}  
}
