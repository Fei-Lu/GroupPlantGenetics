/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import com.koloboke.collect.map.hash.HashDoubleIntMap;
import com.koloboke.collect.map.hash.HashDoubleIntMaps;
import com.koloboke.collect.map.hash.HashFloatIntMap;
import com.koloboke.collect.map.hash.HashFloatIntMaps;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
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
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import pgl.infra.range.Range;
import pgl.infra.table.ColumnTable;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import rcaller.RCaller;
import rcaller.RCode;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author xujun
 */
public class SiPAS {
    int[] barcodeLengths = null;
    List<String> indexList = null;
    List<String>[] barcodeLists = null;

    HashMap<String, String>[] barcodeMethodMaps = null;
    //    HashMap<Integer,Integer>
    HashMap barcodeStrain = new HashMap();
    HashMap indexStrain = new HashMap();
    List<String>[] methodLists = null;

    List<String> allMethodList = new ArrayList<String>();
    String inputDirS=null;
    String outputDirS=null;
    int geneNumber=92;

    public SiPAS(){//String parameterFileS
//        this.parseFq();
//        this.parseFqByBarcodePlateseq();//两端文件都读入 但是只分出了不带ployT的那一端！
//            this.checkFq();
//        this.ployTLengthAndQUMI();//注意有没有去除barcode
//        this.ployTLengthAndQ();
//        this.averageQ();
//        this.HTSeqCountPair();
//        this.parseFqByIndexAndBarcode();
//        this.HTSeqCountSingle();
//        this.HTSeqCountSingle();
//        this.HTSeqCountMerge();
//        this.HTSeqCountMergeRPKM();
//        this.parseFqByIndex();
//        this.findIndex();
//        this.removeUnmappedStartEnd();
//        this.avergeTranscriptLengthInWheat();
//        this.splitSample();//双端文件读入，R1和R2都输出，用于双端比对
//        this.dataVolumn();//使用seqtk对fq文件进行双端reads的抽取
//        this.findTriad();
//        this.shScoreParse();
//        this.getList();
//        this.highDiverganceGene();
//        this.isSame();
//        this.geneNumber();
//        this.UMIGo();
//        this.SiPASReverse();
//        this.changeQ();
//        this.annotatedBAM();
//        this.removeDuplicate();
//        this.countToBed();
//        this.removeBandT();
    }

    public void removeBandTUMI(){
//         String inputFile="/data1/home/junxu/BRBUMI/UMItest/SiPAS/difValumn/subFastqs/500000020191226am_R1.fq";
//         String outputFile="/data1/home/junxu/BRBUMI/UMItest/SiPAS/difValumn/5000000_20191226am_trim_R1.fq.gz";
        String inputDirS = "/data2/junxu/SiPASResult/200415/SiPASU/subFastqs/";
        String outputDirS="/data2/junxu/SiPASResult/200415/SiPASU/withoutBandT/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f->{
            try{
                BufferedReader br1 = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outputDirS+f.getName());
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outputDirS+f.getName().replace("_R1.fq.gz", "_R2.fq.gz"));
                String temp=null;String seq=null; int a=0;
                while((temp=br1.readLine())!=null){
                    seq=br1.readLine();
                    for(int i=23;i<seq.length();i++){//UMI的需要改动
                        if(seq.charAt(i)!='T'){
                            a++;
                        }else{
                            a=0;//这里如果等于0的话对ployT的定义更加严格 不等于0的话就宽松一点儿
                        }
                        if(a>2){
                            bw1.write(temp+":"+seq.substring(8,18));bw1.newLine();
                            bw1.write(seq.substring(i-3));bw1.newLine();
                            bw1.write(br1.readLine());bw1.newLine();
                            bw1.write(br1.readLine().substring(i-3));bw1.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            break;
                        }
                    }
                    if(a<=2){
                        br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                }
                br1.close();br2.close();
                bw1.flush();bw1.close();
                bw2.flush();bw2.close();
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        });

    }
    public void removeBandT(){
//         String inputFile="/data1/home/junxu/BRBUMI/UMItest/SiPAS/difValumn/subFastqs/500000020191226am_R1.fq";
//         String outputFile="/data1/home/junxu/BRBUMI/UMItest/SiPAS/difValumn/5000000_20191226am_trim_R1.fq.gz";
        String inputDirS = "/data2/junxu/SiPASResult/200415/SiPASU/subFastqs/";
        String outputDirS="/data2/junxu/SiPASResult/200415/SiPASU/withoutBandT/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f->{
            try{
                BufferedReader br1 = IOUtils.getTextGzipReader(f.getAbsolutePath());
                BufferedReader br2 = IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                BufferedWriter bw1 = IOUtils.getTextGzipWriter(outputDirS+f.getName());
                BufferedWriter bw2 = IOUtils.getTextGzipWriter(outputDirS+f.getName().replace("_R1.fq.gz", "_R2.fq.gz"));
                String temp=null;String seq=null; int a=0;
                while((temp=br1.readLine())!=null){
                    seq=br1.readLine();
                    for(int i=23;i<seq.length();i++){//UMI的需要改动
                        if(seq.charAt(i)!='T'){
                            a++;
                        }else{
                            a=0;//这里如果等于0的话对ployT的定义更加严格 不等于0的话就宽松一点儿
                        }
                        if(a>2){
                            bw1.write(temp);bw1.newLine();
                            bw1.write(seq.substring(i-3));bw1.newLine();
                            bw1.write(br1.readLine());bw1.newLine();
                            bw1.write(br1.readLine().substring(i-3));bw1.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            bw2.write(br2.readLine());bw2.newLine();
                            break;
                        }
                    }
                    if(a<=2){
                        br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }
                }
                br1.close();br2.close();
                bw1.flush();bw1.close();
                bw2.flush();bw2.close();
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        });

    }

    public void getRawSingle(){
        String inputDirS="/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/data_P101SC18112845-01-F004-B4-21/1.rawdata/S20190118S1-7-20190118S1-7_FKDL190665547-1a/";
        String outputDirS ="/data1/home/junxu/eQTL/phenotype/";
//         RowTable rt = new RowTable("/data2/junxu/SiPASData/index.txt");
//         HashMap indexName = new HashMap();
//         List<String> nameList = new ArrayList();
//         for (int i=0;i<rt.getRowNumber();i++){
//             indexName.put(rt.getCellAsString(i, 1), rt.getCellAsString(i, 0));
//             nameList.add(rt.getCellAsString(i, 0));
//         }
//        File[] subDirS = new File(inputDirS).listFiles();
//        HashSet<String> nameSet = new HashSet();
//        for(int i =0;i<subDirS.length;i++){
//            if(subDirS[i].isDirectory()){
//                File[] fs = subDirS[i].listFiles();
//                fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
//                for(int j=0;j<fs.length;j++){
//                    if (fs[j].isHidden()) continue;
//                    nameSet.add(fs[j].getName().replace("_"+fs[j].getName().split("_")[fs[j].getName().split("_").length-1],""));
//                }
//            }
//        }
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().replace("_"+fs[i].getName().split("_")[fs[i].getName().split("_").length-1],""));
        }
//        List<File> fList = Arrays.asList(fs);
//        BufferedWriter [] bw = new BufferedWriter[nameList.size()];
//        BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
//        for(int i=0;i<nameList.size();i++){
//            bw[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R1.fq.gz");
//            bw1[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R2.fq.gz");
//        }
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS,p+"_1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_2.fq.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw=IOUtils.getTextGzipWriter(outputDirS+"/"+p+"_R1.fq.gz");
                BufferedWriter bw1=IOUtils.getTextGzipWriter(outputDirS+"/"+p+"_R2.fq.gz");
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String index=temp.split(":")[9];
                    if(index.equals("TCCCGA")){
                        bw.write(temp);bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw.write(br.readLine());bw.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                        bw1.write(br1.readLine());bw1.newLine();
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                }
                br.close();br1.close();
                bw.flush();bw.close();
                bw1.flush();bw1.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void getRaw(){
        String inputDirS="/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/data_P101SC18112845-01-F004-B4-21/1.rawdata/";
        String outputDirS ="/data2/junxu/SiPASData/20190319/Rawdata/";
        RowTable rt = new RowTable("/data2/junxu/SiPASData/index.txt");
        HashMap indexName = new HashMap();
        List<String> nameList = new ArrayList();
        for (int i=0;i<rt.getRowNumber();i++){
            indexName.put(rt.getCellAsString(i, 1), rt.getCellAsString(i, 0));
            nameList.add(rt.getCellAsString(i, 0));
        }
        File[] subDirS = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for(int i =0;i<subDirS.length;i++){
            if(subDirS[i].isDirectory()){
                File[] fs = subDirS[i].listFiles();
                fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
                for(int j=0;j<fs.length;j++){
                    if (fs[j].isHidden()) continue;
                    nameSet.add(fs[j].getName().replace("_"+fs[j].getName().split("_")[fs[j].getName().split("_").length-1],""));
                }
            }
        }
//        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_1.fq.gz");
//        HashSet<String> nameSet = new HashSet();
//        for (int i = 0; i < fs.length; i++) {
//            if (fs[i].isHidden()) continue;
//            nameSet.add(fs[i].getName().split("_")[0]);
//        }
//        List<File> fList = Arrays.asList(fs);
        BufferedWriter [] bw = new BufferedWriter[nameList.size()];
        BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
        for(int i=0;i<nameList.size();i++){
            bw[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R1.fq.gz");
            bw1[i]=IOUtils.getTextGzipWriter(outputDirS+"/"+nameList.get(i)+"_R2.fq.gz");
        }
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"/"+p+"_1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"/"+p+"_2.fq.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    String index=temp.split(":")[9];
                    if(indexName.get(index)!=null){
                        int num =nameList.indexOf(indexName.get(index));
                        bw[num].write(temp);bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw[num].write(br.readLine());bw[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                        bw1[num].write(br1.readLine());bw1[num].newLine();
                    }else{
                        br.readLine();br.readLine();br.readLine();
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                    }
                }
                br.close();br1.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try{
            for (int i=0;i<nameList.size();i++){
                bw[i].flush();bw[i].close();
                bw1[i].flush();bw1[i].close();
            }
        }
        catch(Exception ex){
            ex.printStackTrace();
        }
    }

    public void dataCommand() {
        String subFqDirS = new File ("/data1/home/junxu/wheat/SiPAS1.1/test/subFastqs/").getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, "R1.fq.gz");
        for(int i=0;i<fs.length;i++){
//            if (fs[i].length() < 450000000) continue; // Assume 19,000 genes expressed with 1X coverage, to ignore low-quality sample
            if (fs[i].length() < 159000000) continue;
            fList.add(fs[i]);
        }
        try{
            BufferedWriter bw =IOUtils.getTextWriter("/data1/home/junxu/wheat/SiPAS1.1/test/2M/file.txt");
            for(int i=0;i<fList.size();i++){
                bw.write(fList.get(i).getName().replace("_R1.fq.gz",""));
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception ex){
            ex.printStackTrace();
        }
    }
    public void subFastqRatio(){
        String inputDirS="/data2/junxu/SiPASResult/SiPASReverse/subFastqs/";
        String outputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/subFastqs/";
//        String inputDirS="/Users/xujun/Desktop/eQTL/N344/";
//        String outputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/subFastqs/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
//        fs = IOUtils.listFilesStartsWith(fs, "homoGene.nominals");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        RowTable rtAm = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/subFastqs/countam.txt");
        HashMap hmAm = new HashMap ();
        for(int i=0;i<rtAm.getRowNumber();i++){
            hmAm.put("am"+(i+1), rtAm.getCell(i, 0));
        }
        RowTable rtPm = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/subFastqs/countpm.txt");
        HashMap hmPm = new HashMap ();
        for(int i=0;i<rtPm.getRowNumber();i++){
            hmPm.put("pm"+(i+1), rtAm.getCell(i, 0));
        }
        try{
            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
            for(int i=0;i<fList.size();i++){
                if((i+1)/2==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                    fL1.add(fList.get(i));
                }else{
                    fL1.add(fList.get(i));
                    fL1.parallelStream().forEach(f -> {
                        try{
                            int num=0;
                            BufferedWriter bwR1 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName()).getAbsolutePath());
                            if(hmAm.get(f.getName().split("_")[0].split("-")[1])!=null){
                                num=Integer.valueOf(hmAm.get(f.getName().split("_")[0].split("-")[1]).toString());
                            }else{
                                num=Integer.valueOf(hmPm.get(f.getName().split("_")[0].split("-")[1]).toString());
                            }
                            int fqNumber=num/4;
                            List<Integer> total =new ArrayList();
                            for(int a=1;a<fqNumber+1;a++){
                                total.add(a);
                            }
                            Collections.shuffle(total);
                            int sub = (int)Math.round(fqNumber*0.4);
                            System.out.println(f.getName()+"\t"+fqNumber+"\t"+sub);
                            List subList = total.subList(0, sub+1);
                            String tempR1=null; int rowR1=1;
                            BufferedReader brR1 =IOUtils.getTextGzipReader(f.getAbsolutePath());
                            while((tempR1=brR1.readLine())!=null){
                                if(subList.contains(rowR1)){
                                    bwR1.write(tempR1);bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                    bwR1.write(brR1.readLine());bwR1.newLine();
                                }
                                rowR1++;
                            }
                            brR1.close();
                            bwR1.flush();bwR1.close();

                            BufferedWriter bwR2 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().replace("_R1", "_R2")).getAbsolutePath());
                            String tempR2=null; int rowR2=1;
                            BufferedReader brR2 =IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                            while((tempR2=brR2.readLine())!=null){
                                if(subList.contains(rowR2)){
                                    bwR2.write(tempR2);bwR1.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                    bwR2.write(brR2.readLine());bwR2.newLine();
                                }
                                rowR2++;
                            }
                            brR2.close();
                            bwR2.flush();bwR2.close();
                        }
                        catch(Exception ex){
                            ex.getStackTrace();
                        }
                    });
                    count++;
                    fL1.clear();
                }
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }
    public void subFastqNum(){
//        String inputDirS="/data2/junxu/SiPASResult/SiPASAm2Pm1/difVolumn/";
//        String inputDirS="/data2/junxu/SiPASResult/SiPASUR/difVolumn/";
//        String inputDirS="/data2/junxu/SiPASResult/SiPASReverse/2G/difVolumn/";
        String inputDirS="/data1/home/junxu/wheat/SiPAS1.1/test/";
        String outputDirS="/data1/home/junxu/wheat/SiPAS1.1/testR/repTest/subFastqs";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASReverse/2G/difVolumn/count1.txt");
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASAm2Pm1/difVolumn/count.txt");
        RowTable rt = new RowTable("/data1/home/junxu/wheat/SiPAS1.1/test/count.txt");
//        RowTable rt = new RowTable("/data2/junxu/SiPASResult/SiPASUR/difVolumn/count.txt");
        HashMap hm = new HashMap ();
        for(int i=0;i<rt.getRowNumber();i++){
            hm.put(rt.getCell(i, 0), rt.getCell(i, 1));
        }
        try{
            int count=0;List<File> fL1 = new ArrayList(Arrays.asList());
            List<Integer>total =new ArrayList();
            for(int i=0;i<fList.size();i++){
                if((i+1)/1==count && i!=fList.size()-1){//实现了同时跑10个HTSeq
                    fL1.add(fList.get(i));
                }else{
                    fL1.add(fList.get(i));
                    fL1.parallelStream().forEach(f -> {
                        try{
                            int num=0;
                            if(hm.get(f.getName().split("_")[0])!=null){
                                num=Integer.valueOf(hm.get(f.getName().split("_")[0]).toString());
                            }
                            int fqNumber=num/4;
                            for(int a=1;a<fqNumber+1;a++){
                                total.add(a);
                            }
                            for(int rep=1;rep<4;rep++){
                                Collections.shuffle(total);
                                int sub = 2000000;
                                List<Integer> subList = total.subList(0, sub);
                                Set<Integer> set = new HashSet<Integer>(subList);
                                BufferedWriter bwR1 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().split("_")[0]+rep+"_R1.fq.gz").getAbsolutePath());
                                String tempR1=null; int rowR1=1;
                                BufferedReader brR1 =IOUtils.getTextGzipReader(f.getAbsolutePath());
                                while((tempR1=brR1.readLine())!=null){
                                    if(set.contains(rowR1)){
                                        bwR1.write(tempR1);bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                        bwR1.write(brR1.readLine());bwR1.newLine();
                                    }else{
                                        brR1.readLine();brR1.readLine();brR1.readLine();
                                    }
                                    rowR1++;
                                }
                                brR1.close();
                                bwR1.flush();bwR1.close();

                                BufferedWriter bwR2 =IOUtils.getTextGzipWriter(new File(outputDirS,f.getName().split("_")[0]+rep+"_R2.fq.gz").getAbsolutePath());
                                String tempR2=null; int rowR2=1;
                                BufferedReader brR2 =IOUtils.getTextGzipReader(f.getAbsolutePath().replace("_R1.fq.gz", "_R2.fq.gz"));
                                while((tempR2=brR2.readLine())!=null){
                                    if(set.contains(rowR2)){
                                        bwR2.write(tempR2);bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                        bwR2.write(brR2.readLine());bwR2.newLine();
                                    }else{
                                        brR2.readLine();brR2.readLine();brR2.readLine();
                                    }
                                    rowR2++;
                                }
                                brR2.close();
                                bwR2.flush();bwR2.close();
                            }
                        }
                        catch(Exception ex){
                            ex.getStackTrace();
                        }
                    });
                    count++;
                    fL1.clear();
                }
            }

        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }

    public String tail( File file, int lines) {
        java.io.RandomAccessFile fileHandler = null;
        try {
            fileHandler =
                    new java.io.RandomAccessFile( file, "r" );
            long fileLength = fileHandler.length() - 1;
            StringBuilder sb = new StringBuilder();
            int line = 0;
            for(long filePointer = fileLength; filePointer != -1; filePointer--){
                fileHandler.seek( filePointer );
                int readByte = fileHandler.readByte();
                if( readByte == 0xA ) {
                    if (filePointer < fileLength) {
                        line = line + 1;
                    }
                } else if( readByte == 0xD ) {
                    if (filePointer < fileLength-1) {
                        line = line + 1;
                    }
                }
                if (line >= lines) {
                    break;
                }
                sb.append( ( char ) readByte );
            }
            String lastLine = sb.reverse().toString();
            return lastLine;
        } catch( java.io.FileNotFoundException e ) {
            e.printStackTrace();
            return null;
        } catch( java.io.IOException e ) {
            e.printStackTrace();
            return null;
        }
        finally {
            if (fileHandler != null )
                try {
                    fileHandler.close();
                } catch (IOException e) {
                }
        }
    }
    public void parseBed (){
        String inputFile="/data1/home/junxu/eQTL/7FastQTL/7_nor_box_sorted.bed.gz";
        String outputDirS="/data1/home/junxu/eQTL/7FastQTL/7ParseBed/";
        try{
            BufferedReader br = IOUtils.getTextGzipReader(inputFile);
            String temp=null;
            String [] tem = null;List<String> tList=new ArrayList();
            try{
                temp=br.readLine();
                BufferedWriter[] bw = new BufferedWriter[45];
                for(int i=0;i<45;i++){
                    bw[i] = IOUtils.getTextWriter(new File(outputDirS,i+".sorted.bed").getAbsolutePath());
                    bw[i].write(temp);bw[i].newLine();
                }
                while((temp=br.readLine())!=null){
                    tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    bw[Integer.valueOf(tem[0])].write(temp);bw[Integer.valueOf(tem[0])].newLine();
                }
                br.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                }
            }
            catch(Exception ex){
                ex.getStackTrace();
            }
        }
        catch(Exception ex){
            ex.getStackTrace();
        }
    }


    public void countToBed(){
//         String sampleInfor="/data1/home/junxu/eQTL/FastQTL2/P7_sampleName.txt";
        String inputFiles="/Users/xujun/Desktop/eQTL/N344/7_nor_box.txt";
        String outputFiles="/Users/xujun/Desktop/eQTL/N344/P7_nor_box.bed";
        GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/eQTL/N344/wheat_v1.1_Lulab.gff3");

        gf.isWithinThisGene(geneNumber, geneNumber, geneNumber);
        try{
            BufferedReader br = IOUtils.getTextReader(inputFiles);
            BufferedWriter bw = IOUtils.getTextWriter(outputFiles);
            String temp=null; String [] tem = null;
            StringBuilder header = new StringBuilder();
            List<String> tList=new ArrayList();HashSet nameList=new HashSet();
            temp=br.readLine();tList= PStringUtils.fastSplit(temp);
            tem = tList.toArray(new String[tList.size()]);
            header.append("#chr"+"\t"+"start"+"\t"+"end"+"\t"+"gene"+"\t");
            for(int i=1;i<tem.length;i++){
                header.append(tem[i]+"\t");
            }
            bw.write(header.toString().replaceAll("\\s+$", ""));
            bw.newLine();
            while((temp=br.readLine())!=null){
                tList= PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                StringBuilder out = new StringBuilder();
                int index=gf.getGeneIndex(tem[0]);
                out.append(gf.getGeneChromosome(index)+"\t"+gf.getGeneStart(index)+"\t"+gf.getGeneEnd(index)+"\t");
                for(int i =0;i<tem.length;i++){
                    out.append(tem[i]+"\t");
                }
                bw.write(out.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            br.close();
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

    public static double sd(double a[],double n) {
        double sum = 0; double mean=0;
        for (int i = 0; i < n; i++)  {
            sum+=a[i];
        }
        mean=sum/n;
        sum=0;
        for (int i = 0; i < n; i++){
            sum+=Math.pow(a[i]-mean, 2);
        }
        return Math.sqrt(sum / (n -1));
    }
    public void annotatedBAM(){
        String inputDirS="/data2/junxu/SiPASResult/200415/SE/";
        String fqDirS="/data2/junxu/SiPASResult/200415/SiPASU/downSample/without/subFastqs/";
//        String inputDirS="/Users/xujun/Desktop/";
//        String fqDirS="/Users/xujun/Desktop/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.bam");
        fs=IOUtils.listFilesContains(fs,"U");
        List<File> fList = Arrays.asList(fs);
        int count=0;
        List<File> fSet1 = new ArrayList(Arrays.asList());
        for(int i=0;i<fList.size();i++){
            if((i+1)/10==count && i!=fList.size()-1){//实现了n个同时跑
                fSet1.add(fList.get(i));
            }else{
                fSet1.add(fList.get(i));
                fSet1.parallelStream().forEach(f -> {
                    try{
                        BufferedReader br =IOUtils.getTextGzipReader(fqDirS+f.getName().replace("Aligned.out.bam", "R1.fq.gz"));
                        HashMap<String,String>[] hm = new HashMap[11];
                        HashSet<String>[] hs=new HashSet[11];
                        for(int j=0;j<hm.length;j++){
                            hm[j]=new HashMap();
                            hs[j]=new HashSet();
                        }
                        String temp=null;String QNAME=null;String currentUMI=null;
                        while((temp=br.readLine())!=null){
                            QNAME=temp.split(" ")[0].replace("@", "");
                            currentUMI=temp.split("_")[1];
                            int num=10-currentUMI.replaceAll("A", "").length();
                            hm[num].put(QNAME, currentUMI);
                            hs[num].add(QNAME);
                            br.readLine();br.readLine();br.readLine();
                        }
                        SamReader sr = SamReaderFactory.makeDefault().validationStringency(ValidationStringency.SILENT).open(f);
                        SAMFileHeader header=sr.getFileHeader();
                        SAMFileWriterFactory samWriterFactory = new SAMFileWriterFactory();
                        SAMFileWriter samWriter = samWriterFactory.makeBAMWriter(header, false, new File(inputDirS + f.getName().replace("Aligned.out.bam", "annotated.bam")));
                        SAMRecordIterator r = sr.iterator();
                        while(r.hasNext()) {
                            SAMRecord tem=r.next();
//                            System.out.print(tem.getReadName());
                            for(int j =0;j<hs.length;j++){
                                if(hs[j].contains(tem.getReadName())){
                                    tem.setAttribute("UI", hm[j].get(tem.getReadName()));
                                }
                            }
                            samWriter.addAlignment(tem);
                        }
                        samWriter.close();
                        System.out.println("Finished"+f);
                    }
                    catch (Exception e) {
                        e.printStackTrace();
                    }

                });
                count++;
                fSet1.clear();
            }
        }
    }

    public void SiPASReverse(){
//        String inputDirS="/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/Rawdata/";
//        String outputDirS="/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/result/";
        String inputDirS="/data2/junxu/SiPASData/15/Rawdata/R/";
        String outputDirS="/data1/home/junxu/BRBUMI/SiPASReverse/subFastqs/";
        File[] subDirS = new File(inputDirS).listFiles();
        HashSet<String> nameList = new HashSet();
        for(int i =0;i<subDirS.length;i++){
            if(subDirS[i].isDirectory()){
                File[] fs = subDirS[i].listFiles();
                fs = IOUtils.listFilesEndsWith(fs, ".fq.gz");
                for(int j=0;j<fs.length;j++){
                    if (fs[j].isHidden()) continue;
                    nameList.add(fs[j].getName().split("_")[0]);
                }
            }
        }
//        RowTable rt =new RowTable("/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/IndexBarcode.txt");
        RowTable rt =new RowTable("/data1/home/junxu/BRBUMI/IndexBarcode_15.txt");
        HashSet indexSet = new HashSet ();
        HashMap<String, BufferedWriter> barcodeWriterMap1 = new HashMap<>();
        HashMap<String, BufferedWriter> barcodeWriterMap2 = new HashMap<>();
        BufferedWriter[] bws1 = new BufferedWriter[rt.getRowNumber()];
        BufferedWriter[] bws2 = new BufferedWriter[rt.getRowNumber()];
        for(int i =0;i<rt.getRowNumber();i++){
            bws1[i] = IOUtils.getTextGzipWriter(new File(outputDirS+rt.getCell(i, 0)+"-"+rt.getCell(i, 3)+"_R1.fq.gz").getAbsolutePath());
            bws2[i] = IOUtils.getTextGzipWriter(new File(outputDirS+rt.getCell(i, 0)+"-"+rt.getCell(i, 3)+"_R2.fq.gz").getAbsolutePath());
            barcodeWriterMap1.put(rt.getCell(i, 1)+"_"+rt.getCell(i, 2), bws1[i]);
            barcodeWriterMap2.put(rt.getCell(i, 1)+"_"+rt.getCell(i, 2), bws2[i]);
            indexSet.add(rt.getCell(i, 1));
        }
        int count=0;HashSet <String> nameSet1 = new HashSet();
        String [] names=null;
        names = nameList.toArray(new String[nameList.size()]);
        BufferedWriter bwindex =IOUtils.getTextGzipWriter(new File(outputDirS+"error_index.txt.gz").getAbsolutePath());
        BufferedWriter bwRam =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Ram.txt.gz").getAbsolutePath());
        BufferedWriter bwRpm =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Rpm.txt.gz").getAbsolutePath());
        BufferedWriter bwR =IOUtils.getTextGzipWriter(new File(outputDirS+"error_R.txt.gz").getAbsolutePath());//index和barcode都没有匹配
        for(int i=0;i<names.length;i++){
            if((i+1)/2==count && i!=nameList.size()-1){//实现了同时跑10个HTSeq
                nameSet1.add(names[i]);
            }else{
                nameSet1.add(names[i]);
                nameSet1.parallelStream().forEach(f -> {
                    BufferedReader br1 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R2.fq.gz");
                    BufferedReader br2 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R1.fq.gz");
                    String temp =null;String seq=null; String q=null;
                    String currentIndex=null;String currentBarcode=null; String currentUMI=null;
                    BufferedWriter tw = null;
                    try{
                        while((temp=br1.readLine())!=null){
                            currentIndex=temp.split(":")[9];
                            if(!indexSet.contains(currentIndex)){
                                bwindex.write(temp);bwindex.newLine();
                                bwindex.write(br1.readLine());bwindex.newLine();
                                bwindex.write(br1.readLine());bwindex.newLine();
                                bwindex.write(br1.readLine());bwindex.newLine();
                                br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                                continue;
                            }
                            seq = br1.readLine();
                            currentBarcode = seq.substring(0, 8);
                            if (barcodeWriterMap1.get(currentIndex+"_"+currentBarcode)!=null) {
                                tw=barcodeWriterMap1.get(currentIndex+"_"+currentBarcode);
                                tw.write(temp);tw.newLine();
                                tw.write(seq);tw.newLine();
                                tw.write(br1.readLine());tw.newLine();
                                tw.write(br1.readLine());tw.newLine();
                                tw=barcodeWriterMap2.get(currentIndex+"_"+currentBarcode);
                                temp=br2.readLine();seq=br2.readLine();
                                tw.write(temp);tw.newLine();
                                tw.write(seq);tw.newLine();
                                tw.write(br2.readLine());tw.newLine();
                                tw.write(br2.readLine());tw.newLine();
                            }else {
                                if(currentIndex=="GGCTAC" || currentIndex=="CGGAAT" || currentIndex=="CTAGCT" || currentIndex=="CTATAC"){
                                    bwRam.write(temp);bwRam.newLine();
                                    bwRam.write(seq);bwRam.newLine();
                                    bwRam.write(br1.readLine());bwRam.newLine();
                                    bwRam.write(br1.readLine());bwRam.newLine();
                                }else{
                                    if(currentIndex=="CTCAGA" || currentIndex=="GACGAC"|| currentIndex=="TAATCG" || currentIndex=="TACAGC"){
                                        bwRpm.write(temp);bwRpm.newLine();
                                        bwRpm.write(seq);bwRpm.newLine();
                                        bwRpm.write(br1.readLine());bwRpm.newLine();
                                        bwRpm.write(br1.readLine());bwRpm.newLine();
                                    }else{
                                        bwR.write(temp);bwR.newLine();
                                        bwR.write(seq);bwR.newLine();
                                        bwR.write(br1.readLine());bwR.newLine();
                                        bwR.write(br1.readLine());bwR.newLine();
                                    }
                                }
                                br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                            }
                        }
                        br1.close();br2.close();
                    }
                    catch (Exception e) {
                        //                System.out.println(tw.toString());
                        e.printStackTrace();
                    }
                });
                count++;
                nameSet1.clear();
            }
        }
        try{
            for(int i=0;i<bws1.length;i++){
                bws1[i].flush();bws1[i].close();
                bws2[i].flush();bws2[i].close();
            }
            bwindex.flush();bwindex.close();
            bwRam.flush();bwRam.close();bwRpm.flush();bwRpm.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void UMIGo(){
//        String inputDirS="/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/Rawdata/";
//        String outputDirS="/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/result/";
        String inputDirS="/data2/junxu/SiPASData/ANNO_ANCBJ170529_PM-ANCBJ170529-14_2020-03-10/temp/";
        String outputDirS="/data2/junxu/SiPASData/ANNO_ANCBJ170529_PM-ANCBJ170529-14_2020-03-10/parse/";
        File[] subDirS = new File(inputDirS).listFiles();
        HashSet<String> nameList = new HashSet();
        for(int i =0;i<subDirS.length;i++){
            if(subDirS[i].isDirectory()){
                File[] fs = subDirS[i].listFiles();
                fs = IOUtils.listFilesEndsWith(fs, ".fq.gz");
                for(int j=0;j<fs.length;j++){
                    if (fs[j].isHidden()) continue;
                    nameList.add(fs[j].getName().split("_")[0]);
                }
            }
        }
//        RowTable rt =new RowTable("/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/IndexBarcode.txt");
        RowTable rt =new RowTable("/data1/home/junxu/BRBUMI/UMItest/IndexBarcode.txt");
        HashSet indexSet = new HashSet ();
        HashMap<String, BufferedWriter> barcodeWriterMap1 = new HashMap<>();
        HashMap<String, BufferedWriter> barcodeWriterMap2 = new HashMap<>();
        BufferedWriter[] bws1 = new BufferedWriter[rt.getRowNumber()];
        BufferedWriter[] bws2 = new BufferedWriter[rt.getRowNumber()];
        for(int i =0;i<rt.getRowNumber();i++){
            bws1[i] = IOUtils.getTextGzipWriter(new File(outputDirS+rt.getCell(i, 0)+"-"+rt.getCell(i, 3)+"_R1.fq.gz").getAbsolutePath());
            bws2[i] = IOUtils.getTextGzipWriter(new File(outputDirS+rt.getCell(i, 0)+"-"+rt.getCell(i, 3)+"_R2.fq.gz").getAbsolutePath());
            barcodeWriterMap1.put(rt.getCell(i, 1)+"_"+rt.getCell(i, 2), bws1[i]);
            barcodeWriterMap2.put(rt.getCell(i, 1)+"_"+rt.getCell(i, 2), bws2[i]);
            indexSet.add(rt.getCell(i, 1));
        }
        System.out.println(barcodeWriterMap1.size());
        int count=0;HashSet <String> nameSet1 = new HashSet();
        String [] names=null;
        names = nameList.toArray(new String[nameList.size()]);
        BufferedWriter bwindex =IOUtils.getTextGzipWriter(new File(outputDirS+"error_index.txt.gz").getAbsolutePath());
        BufferedWriter bwRam =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Ram.txt.gz").getAbsolutePath());
        BufferedWriter bwRpm =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Rpm.txt.gz").getAbsolutePath());
        BufferedWriter bwUam =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Uam.txt.gz").getAbsolutePath());
        BufferedWriter bwUpm =IOUtils.getTextGzipWriter(new File(outputDirS+"error_Upm.txt.gz").getAbsolutePath());
        BufferedWriter bwR =IOUtils.getTextGzipWriter(new File(outputDirS+"error_R.txt.gz").getAbsolutePath());//index和barcode都没有匹配
        BufferedWriter bwU =IOUtils.getTextGzipWriter(new File(outputDirS+"error_U.txt.gz").getAbsolutePath());
        for(int i=0;i<names.length;i++){
            if((i+1)/2==count && i!=nameList.size()-1){//实现了同时跑10个HTSeq
                nameSet1.add(names[i]);
            }else{
                nameSet1.add(names[i]);
                nameSet1.parallelStream().forEach(f -> {
                    if(f.toString().contains("R")){
                        BufferedReader br1 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R2.fq.gz");
                        BufferedReader br2 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R1.fq.gz");
                        String temp =null;String seq=null; String q=null;
                        String currentIndex=null;String currentBarcode=null; String currentUMI=null;
                        BufferedWriter tw = null;
                        try{
                            while((temp=br1.readLine())!=null){
                                currentIndex=temp.split(":")[9];
                                if(!indexSet.contains(currentIndex)){
                                    bwindex.write(temp);bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                                    continue;
                                }
                                seq = br1.readLine();
                                currentBarcode = seq.substring(0, 8);
                                currentUMI=seq.substring(8,18);
                                if (barcodeWriterMap1.get(currentIndex+"_"+currentBarcode)!=null) {
                                    tw=barcodeWriterMap1.get(currentIndex+"_"+currentBarcode);
                                    //                            tw.write(temp+":"+currentUMI);tw.newLine();
                                    tw.write(temp+":"+currentBarcode+":"+currentUMI);tw.newLine();
                                    tw.write(currentBarcode+seq.substring(23));tw.newLine();//把ployT前面的UMI、5个V都去掉了
                                    tw.write(br1.readLine());tw.newLine();
                                    q=br1.readLine();
                                    tw.write(q.substring(0, 8)+q.substring(23));tw.newLine();
                                    tw=barcodeWriterMap2.get(currentIndex+"_"+currentBarcode);
                                    temp=br2.readLine();seq=br2.readLine();currentUMI=seq.substring(8,18);
                                    tw.write(temp+":"+currentBarcode+":"+currentUMI);tw.newLine();
                                    tw.write(seq);tw.newLine();
                                    tw.write(br2.readLine());tw.newLine();
                                    tw.write(br2.readLine());tw.newLine();
                                }else {
                                    if(currentIndex=="GACGAC" || currentIndex=="TAATCG"){
                                        bwRam.write(temp);bwRam.newLine();
                                        bwRam.write(seq);bwRam.newLine();
                                        bwRam.write(br1.readLine());bwRam.newLine();
                                        bwRam.write(br1.readLine());bwRam.newLine();
                                    }else{
                                        if(currentIndex=="TACAGC" || currentIndex=="TATAAT"){
                                            bwRpm.write(temp);bwRpm.newLine();
                                            bwRpm.write(seq);bwRpm.newLine();
                                            bwRpm.write(br1.readLine());bwRpm.newLine();
                                            bwRpm.write(br1.readLine());bwRpm.newLine();
                                        }else{
                                            bwR.write(temp);bwR.newLine();
                                            bwR.write(seq);bwR.newLine();
                                            bwR.write(br1.readLine());bwR.newLine();
                                            bwR.write(br1.readLine());bwR.newLine();
                                        }
                                    }
                                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                                }
                            }
                            br1.close();br2.close();
                        }
                        catch (Exception e) {
                            //                System.out.println(tw.toString());
                            e.printStackTrace();
                        }
                    }else{
                        BufferedReader br1 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R1.fq.gz");
                        BufferedReader br2 =IOUtils.getTextGzipReader(inputDirS+f+"/"+f+"_R2.fq.gz");
                        String temp =null;String seq=null; String q=null;
                        String currentIndex=null;String currentBarcode=null; String currentUMI=null;
                        BufferedWriter tw = null;
                        try{
                            while((temp=br1.readLine())!=null){
                                currentIndex=temp.split(":")[9];
                                if(!indexSet.contains(currentIndex)){
                                    bwindex.write(temp);bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    bwindex.write(br1.readLine());bwindex.newLine();
                                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                                    continue;
                                }
                                seq = br1.readLine();
                                currentBarcode = seq.substring(0, 8);
                                currentUMI=seq.substring(8,18);
                                if (barcodeWriterMap1.get(currentIndex+"_"+currentBarcode)!=null) {
                                    tw=barcodeWriterMap1.get(currentIndex+"_"+currentBarcode);
                                    //                            tw.write(temp+":"+currentUMI);tw.newLine();
                                    tw.write(temp+":"+currentBarcode+":"+currentUMI);tw.newLine();
                                    tw.write(currentBarcode+seq.substring(23));tw.newLine();
                                    tw.write(br1.readLine());tw.newLine();
                                    q=br1.readLine();
                                    tw.write(q.substring(0,8)+q.substring(23));tw.newLine();
                                    tw=barcodeWriterMap2.get(currentIndex+"_"+currentBarcode);
                                    temp=br2.readLine();seq=br2.readLine();currentUMI=seq.substring(8,18);
                                    tw.write(temp+":"+currentBarcode+":"+currentUMI);tw.newLine();
                                    tw.write(seq);tw.newLine();
                                    tw.write(br2.readLine());tw.newLine();
                                    tw.write(br2.readLine());tw.newLine();
                                }else {
                                    if(currentIndex=="GGCTAC" || currentIndex=="CGGAAT"){
                                        bwUam.write(temp);bwUam.newLine();
                                        bwUam.write(seq);bwUam.newLine();
                                        bwUam.write(br1.readLine());bwUam.newLine();
                                        bwUam.write(br1.readLine());bwUam.newLine();
                                    }else{
                                        if(currentIndex=="CTAGCT" || currentIndex=="CTATAC"){
                                            bwUpm.write(temp);bwUpm.newLine();
                                            bwUpm.write(seq);bwUpm.newLine();
                                            bwUpm.write(br1.readLine());bwUpm.newLine();
                                            bwUpm.write(br1.readLine());bwUpm.newLine();
                                        }else{
                                            bwU.write(temp);bwU.newLine();
                                            bwU.write(seq);bwU.newLine();
                                            bwU.write(br1.readLine());bwU.newLine();
                                            bwU.write(br1.readLine());bwU.newLine();
                                        }

                                    }
                                    br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                                }
                            }
                            br1.close();br2.close();
                        }
                        catch (Exception e) {
                            //                System.out.println(tw.toString());
                            e.printStackTrace();
                        }
                    }
                });
                count++;
                nameSet1.clear();
            }
        }
        try{
            for(int i=0;i<bws1.length;i++){
                bws1[i].flush();bws1[i].close();
                bws2[i].flush();bws2[i].close();
            }
            bwindex.flush();bwindex.close();
            bwRam.flush();bwRam.close();bwRpm.flush();bwRpm.close();
            bwUam.flush();bwUam.flush();bwUpm.flush();bwUpm.close();
            bwR.flush();bwR.close();bwU.flush();bwU.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void geneNumber(){
        RowTable rt = new RowTable("/Users/xujun/Desktop/tempory/readme.txt");
        HashMap subChrLength=new HashMap();
        for( int i=0;i< rt.getRowNumber();i++){
            subChrLength.put(rt.getCellAsInteger(i, 0), rt.getCellAsInteger(i, 2));
        }
        BufferedReader br =IOUtils.getTextReader(new File("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneDis.txt").getAbsolutePath());
        int [][] geneNumber = new int [45][10];
//         GeneFeature gf = new GeneFeature("/Users/xujun/Desktop/IGVmaterial/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3");
        int erro=0;
        try{
            String [] tem =null;String temp=null;
            while((temp=br.readLine())!=null){
                tem=temp.split("\t");
                if(tem[2].startsWith("gene")){
                    if(Integer.valueOf(tem[0])%2!=0 && Integer.valueOf(tem[0])!=0){
                        int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])+1).toString());
                        double pos = (double)Integer.valueOf(tem[3])/chrLength;
                        bw.write(pos+"");bw.newLine();
                    }else{
                        if(Integer.valueOf(tem[0])!=0){
                            int chrLength=Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])).toString())+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString());
                            double pos =(double) (Integer.valueOf(tem[3])+Integer.valueOf(subChrLength.get(Integer.valueOf(tem[0])-1).toString()))/chrLength;
//                            geneNumber[Integer.valueOf(tem[0])-1][pos]++;
                            bw.write(pos+"");bw.newLine();
                        }
                    }
                }

            }
            bw.flush();bw.close();br.close();
//            for(int i=1;i<geneNumber.length;i=i+2){
//                for(int j=0;j<geneNumber[i].length;j++){
//                    System.out.println(i+"\t"+geneNumber[i][j]);
//                }
//            }
        }
        catch (Exception e) {
            System.out.println(erro);
            e.printStackTrace();
        }
    }

    public void isSame(){
        BufferedReader br1 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAB.txt").getAbsolutePath());
        BufferedReader br2 =IOUtils.getTextReader(new File("/Users/xujun/Desktop/Tempory/homoAD.txt").getAbsolutePath());
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/samA.txt").getAbsolutePath());
        List<String> nameList1 = new ArrayList();
        List<String> nameList2 = new ArrayList();
        try{
            String [] tem =null;String temp=null;
            while((temp=br1.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList1.add(tem[1].split("\t")[1]);
                }
            }
            while((temp=br2.readLine())!=null){
                if(temp.split(" ").length>1){
                    tem=temp.split(" ");
                    nameList2.add(tem[1].split("\t")[1]);
                }
            }
            for(int i=0;i<nameList2.size();i++){
                if(nameList1.contains(nameList2.get(i))){
                    bw.write(nameList2.get(i));bw.newLine();
                }
            }
            br1.close();br2.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getList(){//to get the data from the distance matrix
        ColumnTable ct = new ColumnTable("/Users/xujun/Desktop/Tempory/g.txt");
        BufferedWriter bw =IOUtils.getTextWriter(new File("/Users/xujun/Desktop/Tempory/geneticList.txt").getAbsolutePath());
        try{
            for(int i =0 ; i < ct.getColumnNumber() ; i++){
                for(int j =i ; j < ct.getRowNumber();j++){
                    bw.write(ct.getCellAsString(i, j));bw.newLine();
                }
            }
            bw.flush();bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void shScoreParse(){
        String inputFileS = "/Users/xujun/summarys.txt";
        String outputFileS = "/Users/xujun/shScore.txt";
        String temp=null;String shtemp=null;
        double sh=0;double shM=0;int k=0;int tempK=0;String num=null;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            while((temp=br.readLine())!=null){
                if(temp.startsWith("logical")){
                    bw.write(String.valueOf(k));bw.newLine();
                    sh=0;shM=0;k=0;tempK=0;
                }
                if(temp.startsWith("Silhouette")){
                    tempK=Integer.valueOf(temp.split(" ")[5]);
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();br.readLine();
                    shtemp=br.readLine();
                    if(!shtemp.contains("Individual")){
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }else{
                        br.readLine();shtemp=br.readLine();
                        if(shtemp.contains("  ")){
                            sh=Double.valueOf(shtemp.split("  ")[3]);
                        }else{
                            sh=Double.valueOf(shtemp.split(" ")[3]);
                        }
                    }

                    if(sh>shM){
                        shM=sh;k=tempK;
                        if(k==6){

                        }
                    }
                }
            }
            br.close();
            bw.flush();bw.close();
        }
        catch (Exception e) {
            System.out.println(sh);
            e.printStackTrace();
        }
    }
    public void shScore (){
        String inputfsK="/Users/xujun/Desktop/eQTL/total/Homology/shScore";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".txt");
        for(int i=0;i<fsK.length;i++){
            if(!fsK[i].getName().contains("temp")){
                fListK.add(fsK[i]);
            }
        }
        RCaller callerL = new RCaller();
        callerL.setRscriptExecutable("/usr/local/bin/R");
        RCode rCodeL = new RCode();
        rCodeL.addRCode("library(cluster);library(factoextra);library(vegan);install.packages(\"Runiversal\");");
        callerL.setRCode(rCodeL);
        callerL.runOnly();
        List KPool =new ArrayList();
        fListK.stream().forEach(f -> {
            double sh=0;double shM=0;int k=0;
            for(int i=1;i<11;i++){
                RCaller caller = new RCaller();
                caller.setRscriptExecutable("//usr//local//bin//R");
                RCode rCode = new RCode();
                rCode.addRCode("a<-read.csv(\""+f+"\",sep=\"\\t\",header = F)");
                rCode.addRCode("results=kmeans(a,");
                rCode.addRCode(i+");dis = dist(a^2);");
                rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
                caller.setRCode(rCode);
                caller.runAndReturnResult("sh");
                sh=caller.getParser().getAsDoubleArray("sh")[0];
                if(sh>shM){
                    shM=sh;
                }else{
                    k=i;
                    break;
                }
            }
            KPool.add(k);
        });
        for(int i=0;i<KPool.size();i++){
            System.out.println(KPool.get(i));
        }
    }
    public void countTPM(){
        RCaller caller = new RCaller();
        caller.setRscriptExecutable("/usr/local/bin/R");
        RCode rCode = new RCode();
        rCode.addRCode("library(cluster);library(factoextra);library(vegan);results=kmeans(");
        rCode.addRCode(""+",");
        for(int i=1;i<11;i++){
            rCode.addRCode(i+");dis = dist(a");
            rCode.addRCode(""+")^2;");
        }
        rCode.addRCode("sil = silhouette (results$cluster, dis);sum=summary(sil);sh<-sum$avg.width");
        caller.setRCode(rCode);
        caller.runAndReturnResult("sh");
    }
    public void findTriad(){
        String inputFileS="/Users/xujun/Desktop/eQTL/total/DESeq/TaS-D/750MSiPAS-TPM.txt";
        RowTable rt = new RowTable(inputFileS);
        HashMap AExpre = new HashMap();
        HashMap BExpre = new HashMap();
        HashMap DExpre = new HashMap();
        for (int i=0;i<rt.getRowNumber();i++){
            if(rt.getCellAsDouble(i, 1)>0.3 ){
                if(rt.getCellAsString(i, 0).charAt(8)=='A'){
                    AExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                }else{
                    if(rt.getCellAsString(i, 0).charAt(8)=='B'){
                        BExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                    }else{
                        DExpre.put(rt.getCell(i, 0), rt.getCellAsDouble(i, 1));
                    }
                }
            }
        }
        Iterator iterA = AExpre.entrySet().iterator();
        while (iterA.hasNext()) {
            HashMap.Entry entryA = (HashMap.Entry)iterA.next();
            Object keyA = entryA.getKey();
            Object valA = entryA.getValue();
            Iterator iterB = BExpre.entrySet().iterator();
            while(iterB.hasNext()){
                HashMap.Entry entryB = (HashMap.Entry)iterB.next();
                Object keyB = entryB.getKey();
                if(keyB.toString().charAt(7)==keyA.toString().charAt(7)){
                    Object valB = entryB.getValue();
                    if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valB.toString()))){
                        Iterator iterD = DExpre.entrySet().iterator();
                        HashMap.Entry entryD = (HashMap.Entry)iterD.next();
                        Object keyD = entryD.getKey();
                        if(keyA.toString().charAt(7)==keyD.toString().charAt(7)){
                            Object valD = entryD.getValue();
                            while(iterD.hasNext()){
                                if(isEqual(Double.parseDouble(valA.toString()),Double.parseDouble(valD.toString()))){
                                    System.out.println(keyA.toString()+"\n"+keyB.toString()+"\n"+keyD.toString());
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    public boolean isEqual(double a,double b){
        double eps = 1E-1;
        if (Math.abs(a - b) / a <= eps) {
            return true;
        }else{
            return false;
        }
    }
    public void splitSample(){
//         String inputDirS="/data2/junxu/SiPASData/20200324/Rawdata/20190118aptest/";
//         String outputDirS ="/data1/home/junxu/wheat/SiPAS1.1/testR/subFastqs";
//         RowTable rt = new RowTable("/data1/home/junxu/wheat/SiPAS1.1/testR/SampleInformation_test.txt");
        String inputDirS="/data2/junxu/SiPASData/20200324/Rawdata/20190118S1-7/";
        String outputDirS ="/data1/home/junxu/eQTL/7phenotype/subFastqs";
        RowTable rt = new RowTable("/data1/home/junxu/eQTL/7phenotype/SampleInformation.txt");
        HashMap barcodeName = new HashMap();
        List<String> nameList = new ArrayList();
        for (int i=0;i<rt.getRowNumber();i++){
//             barcodeName.put(rt.getCellAsString(i, 2), rt.getCellAsString(i, 0)+"-"+rt.getCell(i,3));
//             nameList.add(rt.getCellAsString(i, 0)+"-"+rt.getCell(i,3));
            barcodeName.put(rt.getCellAsString(i, 2), rt.getCellAsString(i, 0)+"_"+rt.getCell(i,1));
            nameList.add(rt.getCellAsString(i, 0)+"_"+rt.getCell(i,1));
        }
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"_R1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_R2.fq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[nameList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[nameList.size()];
                BufferedWriter bwE =IOUtils.getTextGzipWriter(new File(outputDirS,"error.txt.gz").getAbsolutePath());
                for (int i = 0; i < nameList.size(); i++) {
                    bw[i]=IOUtils.getTextGzipWriter(new File(outputDirS, nameList.get(i)+"_R1.fq.gz").getAbsolutePath());
                    bw1[i]=IOUtils.getTextGzipWriter(new File(outputDirS, nameList.get(i)+"_R2.fq.gz").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                String temp = null;String seq=null;
                String currentBarcode=null;
                while ((temp = br.readLine()) != null) {
                    seq=br.readLine();
                    if(barcodeName.get(seq.substring(0, 8))!=null){
                        int index = nameList.indexOf(barcodeName.get(seq.substring(0, 8)));
                        bw[index].write(temp);bw[index].newLine();
                        bw[index].write(seq);bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw[index].write(br.readLine());bw[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                        bw1[index].write(br1.readLine());bw1[index].newLine();
                    }else{
                        br.readLine();br.readLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                        bwE.write(br1.readLine());bwE.newLine();
                    }
                }
                br.close();br1.close();
                for(int i=0;i<bw.length;i++){
                    bw[i].flush();bw[i].close();
                    bw1[i].flush();bw1[i].close();
                }

            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        System.out.println("finished");
    }
    public void avergeTranscriptLengthInWheat(){
        String inputFile="/Users/xujun/Desktop/wheat/rightchangeiwgsc_refseqv1.0_HighConf_2017Mar13.gff3";
        GeneFeature gf = new GeneFeature(inputFile);
        long totalTranscriptLength = 0;
        int transcriptNumber = 0;
        for(int i=0;i<gf.genes.length;i++){
            for(int j=0;j<gf.genes[i].ts.size();j++){
                totalTranscriptLength+=gf.genes[i].ts.get(j).getTranscriptEnd()-gf.genes[i].ts.get(j).getTranscriptStart()+1;
                transcriptNumber++;
            }
        }
        System.out.println(totalTranscriptLength/transcriptNumber);
    }
    private void removeUnmappedStartEnd () {
        List<String> diffValumnMethodList = new ArrayList<String>();
        String inputDirS = new File("/data1/home/junxu/analysis0215/test-result/test-withoutERCC").getAbsolutePath();
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sam");
        List<File> fList = Arrays.asList(fs);
        for(int i =0;i<fList.size();i++){
            diffValumnMethodList.add(fList.get(i).getName().replace(".sam", ""));
        }
        Collections.sort(diffValumnMethodList);
        fList.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(inputDirS+"/change"+f.getName());
                String temp = null;
                String [] tem = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("[")) continue;
                    if(temp.startsWith("@")){
                        bw.write(temp);bw.newLine();
                    }else{
                        List<String> tList= FStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        if(!tem[2].equals("*")){
                            bw.write(temp);bw.newLine();
                        }
                    }
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void findIndex(){
        String inputFile="/data1/home/junxu/wheat/splitByIndex/indexSiPAS-all.fq";
        List <String> indexList=new ArrayList<>();
        try{
            BufferedReader br = IOUtils.getTextReader(inputFile);
            BufferedWriter bw =IOUtils.getTextGzipWriter(inputFile.replace(".fq", ".txt"));
            String temp=null;
            while((temp = br.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                if(!indexList.contains(index)){
                    indexList.add(index);
                }
            }
            BufferedReader br1 = IOUtils.getTextReader(inputFile);
            int [] indexNumber = new int[indexList.size()];
            while((temp = br1.readLine()) != null){
                String index=temp.split(" ")[1].split(":")[3];
                indexNumber[indexList.indexOf(index)]++;
            }
            for(int i=0;i<indexList.size();i++){
                bw.write(indexList.get(i)+"\t"+indexNumber[i]);
                bw.newLine();
            }
            br.close();br1.close();
            bw.flush();bw.close();
        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
    }
    public void HTSeqCountMergeRPKM(){
        List <String> nameList=new ArrayList<>();
        List<String> fileList=new ArrayList<>();
//        String subFqDirS = new File (this.outputDirS).getAbsolutePath();
        String subFqDirS = new File ("/Users/xujun/Desktop/TEP/TEPOut/HTSeqLibrary").getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
        List<File> fList = Arrays.asList(fs);
        for(int i=0;i < fList.size();i++){
            fileList.add(fList.get(i).getName().replace("Count.txt", ""));
        }
        Collections.sort(fileList);
        int[][] geneCount = new int[this.geneNumber][]; // chr, geneByChr, taxon
        double [][] RPKM=new double[this.geneNumber][];
        int [] allReadsNumber = new int [fList.size()];
        for (int i = 0; i < geneNumber; i++) {
            geneCount[i] = new int[fList.size()];
            RPKM[i] = new double[fList.size()];
        }
        fList.stream().forEach(f -> {
            String temp=null;String[] tem = null;
            try{
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                while((temp = br.readLine()) != null){
                    List<String> tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
                    if(tem[0].startsWith("ERCC")){
//                    if(tem[0].startsWith("TraesCS")){
                        if(!nameList.contains(tem[0])){
                            nameList.add(tem[0]);
                        }
                        int index=nameList.indexOf(tem[0]);
                        geneCount[index][fileList.indexOf(f.getName().replace("Count.txt", ""))]=Integer.parseInt(tem[1]);
                        allReadsNumber[fileList.indexOf(f.getName().replace("Count.txt", ""))]+=Integer.parseInt(tem[1]);
                    }
                }

            }
            catch (Exception ex) {
                System.out.println(tem[0]+"\t1234");
                ex.printStackTrace();

            }
        });
        String row = null;String [] tem =null;
        HashMap<String,String> geneNameAndLength = new HashMap();
        try{
//                BufferedReader br = IOUtils.getTextReader("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/ERCC92.gtf");
            BufferedReader br = IOUtils.getTextReader("/Users/xujun/Desktop/wheat/ERCC92/ERCC92.gtf");
            while((row = br.readLine()) != null){
                List<String> tList= PStringUtils.fastSplit(row);
                tem = tList.toArray(new String[tList.size()]);
                geneNameAndLength.put(tem[0], tem[4]);
            }

        }
        catch (Exception ex) {
            ex.printStackTrace();
        }
        for(int i=0;i<fList.size();i++){
            for(int j=0;j<this.geneNumber;j++){
                RPKM[j][i]=geneCount[j][i]*1000000000/allReadsNumber[i]/Integer.parseInt(geneNameAndLength.get(nameList.get(j)));
            }
        }
        String outputFileS = new File (this.outputDirS,"RowAndRPKM-countResult.txt").getAbsolutePath();
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            sb.append("Gene"+"\t");
            for(int i=0;i<fileList.size();i++){
                sb.append(fileList.get(i).replace("Count.txt", "")+"\t"+"RPKM"+"\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for(int i=0;i<geneCount.length;i++){
                sb = new StringBuilder();
                for(int j=0;j<fileList.size();j++){
                    if(j==0){
                        sb.append(nameList.get(i)+"\t");
                    }
                    sb.append(geneCount[i][j]+"\t"+RPKM[i][j]+"\t");
                }
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
    public void HTSeqCountMerge(){//顺序是HTSeq里面的顺序
        String inputDirS="/data1/home/junxu/wheat/SiPAS1.1/SiPAS/geneCount/";
        String outputDirS="/data1/home/junxu/wheat/SiPAS1.1/SiPAS/geneCount/";
        List <String> nameList=new ArrayList<>();
        HashSet nameSet = new HashSet();
        List<String> fileList=new ArrayList<>();
        String subFqDirS = new File (inputDirS).getAbsolutePath();
        File[] fs = new File(subFqDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
        List<File> fList = Arrays.asList(fs);
        Collections.sort(fList);
        for(int i=0;i<fList.size();i++){
            fileList.add(fList.get(i).getName());
        }
//    int[][][] geneCount = new int[chrs.length][][]; // chr, geneByChr, taxon
//    double [][][] RPKM=new double[chrs.length][][];
//    int [][] count=new int[92][fList.size()];
        StringBuilder wc = new StringBuilder();
        wc.append("wc -l ").append(fList.get(0));
        String command = wc.toString();
        System.out.println(command);
        try {
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                int length=temp.split(" ").length-1;
                geneNumber=Integer.valueOf(temp.replace(temp.split(" ")[length],"").replaceAll(" ", ""))-31;//单独跑的后果
            }
            p.waitFor();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(geneNumber);
        int [][] count=new int[107891][fList.size()];
        fList.stream().forEach(f -> {
            String temp=null;String[] tem = null;
            try{
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                int rowN=0;
                while((temp = br.readLine()) != null){
                    if(temp.contains("processed")) continue;
                    List<String> tList= PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
//                if(tem[0].startsWith("ERCC")){
                    if(!tem[0].startsWith("__")){
                        if(!nameSet.contains(tem[0])){
                            nameList.add(tem[0]);
                        }
                        nameSet.add(tem[0]);
                        count[rowN][fList.indexOf(f)]=Integer.parseInt(tem[1]);
                        rowN++;
                    }
                }

            }
            catch (Exception ex) {
                System.out.println(tem[0]+"\t1234");
                ex.printStackTrace();

            }
        });

        String outputFileS = new File (outputDirS,"countResult.txt").getAbsolutePath();
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outputFileS);
            sb.append("Gene"+"\t");
            for(int i=0;i<fileList.size();i++){
                sb.append(fileList.get(i).replace("Count.txt", "")+"\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for(int i=0;i<count.length;i++){
                sb = new StringBuilder();
                for(int j=0;j<fileList.size();j++){
                    if(j==0){
                        sb.append(nameList.get(i)+"\t");
                    }
                    sb.append(count[i][j]+"\t");
                }
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
    public void HTSeqCountSingle(){
        String inputDirS = this.inputDirS;
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s no ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat/rightchangewheat.gtf").append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File (this.outputDirS).getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
    public void HTSeqCountPair(){
        String inputDirS = "/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/Out/diffValumnSams-D";
//        String inputDirS = "/Users/xujun/Desktop/TEP/TEPOut/sams";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Aligned.out.sam");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            StringBuilder sb = new StringBuilder();
            sb.append("htseq-count").append(" -m intersection-nonempty -s reverse ");
            sb.append(f);
            sb.append(" /data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/ERCC92.gtf").append(" >> ");
            sb.append(f.getName().replace("Aligned.out.sam", "Count.txt"));
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File ("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC92/Out/HTSeqCount-D").getAbsolutePath());
                String []cmdarry ={"/bin/bash","-c",command};
                Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                p.waitFor();
            }
            catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished"+f);
        });
    }
    public void dataVolumn() {
        String inputDirS = new File ("/data2/junxu/SiPASResult/200415/SiPASU/withoutBandT/").getAbsolutePath();
//        String subFqDirS = "xujun/TEP/TEPOut/subFastqs";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".fq.gz");
        fs = IOUtils.listFilesStartsWith(fs, "SiPASURam");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            int a=fs[i].getName().split("_").length-1;
            nameSet.add(fs[i].getName().replace(fs[i].getName().split("_")[a],""));//存入的name是包含下划线的
        }
        int numCores = Runtime.getRuntime().availableProcessors();
        nameSet.parallelStream().forEach(f -> {
            try {
                File dir = new File(new File ("/data2/junxu/SiPASResult/200415/SiPAS/r2/SiPASUR/subFastqs/").getAbsolutePath());
                for( int i = 1 ;i<=12;i++){
                    StringBuilder sb = new StringBuilder();
                    sb.append("/data1/home/junxu/software/seqtk/seqtk sample -s100 ").append(inputDirS+"/"+f+"R1.fq.gz ");
                    sb.append(1000000*i+" ").append(" | gzip -c > ").append(" "+1000000*i+"_"+f+"R1.fq.gz");
                    String command = sb.toString();
                    System.out.println(command);
                    String []cmdarry ={"/bin/bash","-c",command};
                    Process p=Runtime.getRuntime().exec(cmdarry,null,dir);
                    p.waitFor();

                    StringBuilder sb1 = new StringBuilder();
                    sb1.append("/data1/home/junxu/software/seqtk/seqtk sample -s100 ").append(inputDirS+"/"+f+"R2.fq.gz ");
                    sb1.append(1000000*i+" ").append(" | gzip -c > ").append(" "+1000000*i+"_"+f+"R2.fq.gz");
                    String command1 = sb1.toString();
                    System.out.println(command1);
                    String []cmdarry1 ={"/bin/bash","-c",command1};
                    Process p1=Runtime.getRuntime().exec(cmdarry1,null,dir);
                    p1.waitFor();
                }
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
    }
    public void ployTLengthAndQUMI () {
//        String barcodeFileS = "/data1/home/junxu/wheat/RNA-seq20181030.txt";
        String inputDirS ="/data2/junxu/SiPASResult/200415/SiPASU/subFastqs/";
        String outputDirS ="/data2/junxu/SiPASResult/200415/SiPASU/TandQ/";
//        String outputDirS ="/Users/xujun/Desktop/test-data/ANNO_ANCBJ170529_PM-ANCBJ170529-13_2020-01-19/pmResult/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
//        fs = IOUtils.listFilesStartsWith(fs, "5000000");
//        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            int length=fs[i].getName().split("_").length;
            nameSet.add(fs[i].getName().replace(fs[i].getName().split("_")[length-1], ""));
        }
        nameSet.stream().forEach(p -> {
            String seq1 = null;String temp1=null;String phred1=null;
            String seq2 = null;String temp2=null;String phred2=null;
//            BufferedWriter bw = utils.IOUtils.getTextWriter(new File(outputDirS,p+"highQ.txt").getAbsolutePath());
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(inputDirS+p+"R1.fq.gz");
                BufferedReader br2 = IOUtils.getTextGzipReader(inputDirS+p+"R2.fq.gz");
                BufferedWriter bwL = IOUtils.getTextWriter(new File(outputDirS,p+"L.txt").getAbsolutePath());
                BufferedWriter bwQ = IOUtils.getTextWriter(new File(outputDirS,p+"Q.txt").getAbsolutePath());
                BufferedWriter bwtQ = IOUtils.getTextWriter(new File(outputDirS,p+"2Q.txt").getAbsolutePath());
                BufferedWriter bwAQ = IOUtils.getTextWriter(new File(outputDirS,p+"AQ.txt").getAbsolutePath());
                List L=new ArrayList();
                List Q=new ArrayList();
                List tQ=new ArrayList();
                List AQ=new ArrayList();
                while ((temp1 = br1.readLine()) != null){
                    int a=0;int allPhred1 = 0;int allPhred2 = 0;int lq=0;int q=0;
                    seq1=br1.readLine(); br1.readLine();phred1=br1.readLine();
                    br2.readLine();seq2=br2.readLine(); br2.readLine();phred2=br2.readLine();
                    for(int i=23;i<seq1.length();i++){//UMI(barcode/UMI已经去除):0
                        if(seq1.charAt(i)!='T'){
                            a++;
                        }else{
                            a=0;
                        }
                        if(a>2){
                            L.add(i-25);
                            for(int j=i-2;j<phred1.length();j++){
                                allPhred1+=(int)phred1.charAt(j)-33;
                            }
                            lq=152-i;
                            Q.add(allPhred1/lq);//UMI:129-i

                            break;
                        }
                    }
                    if(a<=2){
                        L.add(127);//UMI:150-8-10-5
//                        Q.add(0);
                        allPhred1=0;lq=0;
                    }
                    for(int j=0;j<phred2.length();j++){
                        allPhred2+=(int)phred2.charAt(j)-33;
                    }
                    tQ.add(allPhred2/150);
                    AQ.add((allPhred1+allPhred2)/(150+lq));
                }
                Collections.shuffle(L);
                Collections.shuffle(Q);
                Collections.shuffle(tQ);Collections.shuffle(AQ);
                int randomSeriesLength = 10000;
                List<Integer> randomSeriesL = L.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriesQ = Q.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriestQ = tQ.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriesAQ = AQ.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                    if (k % 25 == 0){
                        bwL.write(randomSeriesL.get(k)+ " ");bwL.newLine();
                        bwQ.write(randomSeriesQ.get(k)+ " ");bwQ.newLine();
                        bwtQ.write(randomSeriestQ.get(k)+ " ");bwtQ.newLine();
                        bwAQ.write(randomSeriesAQ.get(k)+ " ");bwAQ.newLine();
                    } else {
                        bwL.write(randomSeriesL.get(k)+ " ");
                        bwQ.write(randomSeriesQ.get(k)+ " ");
                        bwtQ.write(randomSeriestQ.get(k)+ " ");
                        bwAQ.write(randomSeriesAQ.get(k)+ " ");
                    }
                }
                br1.close();
                bwL.flush();bwL.close();
                bwQ.flush();bwQ.close();
                bwtQ.flush();bwtQ.close();
                bwAQ.flush();bwAQ.close();
//                bw.flush();bw.close();
            }
            catch (Exception ex) {
                System.out.println(seq1+"\t1234");
                ex.printStackTrace();

            }
        });
    }
    public void ployTLengthAndQ () {
//        String barcodeFileS = "/data1/home/junxu/wheat/RNA-seq20181030.txt";
//        String inputDirS ="/data2/junxu/SiPASResult/SiPASReverse/subFastqs/";
        String inputDirS ="/data2/junxu/SiPASResult/200415/SiPAS/subFastqs/";
//        String barcodeFileS = "/Users/xujun/Desktop/RNA-seq20181030.txt";
        String outputDirS ="/data2/junxu/SiPASResult/200415/SiPAS/TandQ/";
//        String outputDirS ="/data2/junxu/SiPASResult/SiPASReverse/TandQ/";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1.fq.gz");
//        fs = IOUtils.listFilesStartsWith(fs, "SiPASam");
//        fs = IOUtils.listFilesStartsWith(fs, "5000000");
//        List<File> fList = Arrays.asList(fs);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            int length=fs[i].getName().split("_").length;
            nameSet.add(fs[i].getName().replace(fs[i].getName().split("_")[length-1], ""));
        }
        nameSet.stream().forEach(p -> {
            String seq1 = null;String temp1=null;String phred1=null;
            String seq2 = null;String temp2=null;String phred2=null;
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(inputDirS+p+"R1.fq.gz");
                BufferedReader br2 = IOUtils.getTextGzipReader(inputDirS+p+"R2.fq.gz");
                BufferedWriter bwL = IOUtils.getTextWriter(new File(outputDirS,p+"L.txt").getAbsolutePath());
                BufferedWriter bwQ = IOUtils.getTextWriter(new File(outputDirS,p+"Q.txt").getAbsolutePath());
                BufferedWriter bwtQ = IOUtils.getTextWriter(new File(outputDirS,p+"2Q.txt").getAbsolutePath());
                BufferedWriter bwAQ = IOUtils.getTextWriter(new File(outputDirS,p+"AQ.txt").getAbsolutePath());
                List L=new ArrayList();
                List Q=new ArrayList();
                List tQ=new ArrayList();
                List AQ=new ArrayList();
                while ((temp1 = br1.readLine()) != null){
                    int a=0;int allPhred1 = 0;int allPhred2 = 0;int lq=0;
                    seq1=br1.readLine(); br1.readLine();phred1=br1.readLine();
                    br2.readLine();seq2=br2.readLine(); br2.readLine();phred2=br2.readLine();
                    for(int i=8;i<seq1.length();i++){//UMI(barcode/UMI已经去除):0 注意barcode是否去除
                        if(seq1.charAt(i)!='T'){
                            a++;
                        }else{
                            a=0;//这里如果等于0的话对ployT的定义更加严格 不等于0的话就宽松一点儿
                        }
                        if(a>2){
                            L.add(i-10);//根据是否含有barcode进行修改
                            for(int j=i-2;j<phred1.length();j++){
                                allPhred1+=(int)phred1.charAt(j)-33;
                            }
                            lq=152-i;
                            Q.add(allPhred1/lq);//UMI:134-i
                            break;
                        }
                    }
                    if(a<=2){
                        L.add(142);
//                        Q.add(0);
                        allPhred1=0;lq=0;
                    }
                    for(int j=0;j<phred2.length();j++){
                        allPhred2+=(int)phred2.charAt(j)-33;
                    }
                    tQ.add(allPhred2/150);//UMI:134-i
                    AQ.add((allPhred1+allPhred2)/(150+lq));
                }
                Collections.shuffle(L);
                Collections.shuffle(Q);
                Collections.shuffle(tQ);Collections.shuffle(AQ);
                int randomSeriesLength = 10000;
                List<Integer> randomSeriesL = L.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriesQ = Q.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriestQ = tQ.subList(0, randomSeriesLength+1);
                List<Integer> randomSeriesAQ = AQ.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                    if (k % 25 == 0){
                        bwL.write(randomSeriesL.get(k)+ " ");bwL.newLine();
                        bwQ.write(randomSeriesQ.get(k)+ " ");bwQ.newLine();
                        bwtQ.write(randomSeriestQ.get(k)+ " ");bwtQ.newLine();
                        bwAQ.write(randomSeriesAQ.get(k)+ " ");bwAQ.newLine();
                    } else {
                        bwL.write(randomSeriesL.get(k)+ " ");
                        bwQ.write(randomSeriesQ.get(k)+ " ");
                        bwtQ.write(randomSeriestQ.get(k)+ " ");
                        bwAQ.write(randomSeriesAQ.get(k)+ " ");
                    }
                }
                br1.close();
                bwL.flush();bwL.close();
                bwQ.flush();bwQ.close();
                bwtQ.flush();bwtQ.close();
                bwAQ.flush();bwAQ.close();
            }
            catch (Exception ex) {
                System.out.println(seq1+"\t1234");
                ex.printStackTrace();

            }
        });
    }
    public void averageQ () {
//        String outputDirS ="/data1/home/junxu/BRBUMI/UMItest/lengthPloyTAndQ";
        String outputDirS ="/Users/xujun/Desktop/test-data/ANNO_ANCBJ170529_PM-ANCBJ170529-13_2020-01-19/pmResult";
//        String inputDirS ="/data1/home/junxu/BRBUMI/UMItest/result";
        String inputDirS ="/Users/xujun/Desktop/test-data/ANNO_ANCBJ170529_PM-ANCBJ170529-13_2020-01-19/pmResult";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R2.fq.gz");
//        fs = IOUtils.listFilesStartsWith(fs, "6000000");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(p -> {
            String seq = null;
            String temp=null;
            String phred=null;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(p.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS,p.getName().replace(".fq.gz", "Q.txt")).getAbsolutePath());
                List Q=new ArrayList();
                while ((temp = br.readLine()) != null){
                    int allPhred=0;br.readLine();br.readLine();phred=br.readLine();
                    for(int j=0;j<phred.length();j++){
                        allPhred+=(int)phred.charAt(j)-33;
                    }
                    Q.add(allPhred/150);
                }
                Collections.shuffle(Q);
                int randomSeriesLength = 10000;
                List<Integer> randomSeries = Q.subList(0, randomSeriesLength+1);
                for (int k = 1;  k<randomSeriesLength+1; k++) {
                    if (k % 25 == 0){
                        bw.write(randomSeries.get(k)+ " ");bw.newLine();
                    } else {
                        bw.write(randomSeries.get(k)+ " ");
                    }
                }
                br.close();
                bw.flush();bw.close();
            }
            catch (Exception ex) {
                System.out.println(seq+"\t1234");
                ex.printStackTrace();

            }
        });
    }
    public void checkFq(){
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".fastq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.stream().forEach(f -> {
            try{
                BufferedReader br = null;
                if (f.getAbsolutePath().endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                System.out.println(f);
                String temp = null;
                int count =0;
                while((temp = br.readLine())!=null){
                    count++;
                    if(temp.startsWith("@")){
                        br.readLine();br.readLine();br.readLine();
                    }
                    else{
                        System.out.println(temp+"\t"+count);
                        break;
                    }
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }
    private void parseFq () {
        long startTimePoint = System.nanoTime();
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/splitByIndex";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            try {
                BufferedWriter[][] bws = null;
                bws = new BufferedWriter[indexList.size()][];
                BufferedWriter [] bw = new BufferedWriter[indexList.size()];
                for (int i = 0; i < indexList.size(); i++) {
                    bws[i]=new BufferedWriter[barcodeLists[i].size()];
                    for(int j=0;j<barcodeLists[i].size();j++){
                        String fileName=barcodeMethodMaps[i].get(barcodeLists[i].get(j));
                        bws[i][j]=IOUtils.getTextWriter(new File(outputDirS, fileName+".fq").getAbsolutePath());
                    }
                    if(barcodeLists[i].size()==0){
                        String fileName=barcodeMethodMaps[i].get(indexList.get(i));
                        bw[i]=IOUtils.getTextWriter(new File(outputDirS, fileName+".fq").getAbsolutePath());
                    }
                }
                BufferedReader br = null;
                if (f.getAbsolutePath().endsWith(".gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                String index1 = null ;String index2 = null ;
                String temp = null;
                String seq = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    index1=temp.split(" ")[1].split(":")[3].substring(1, 6);
                    if(indexList.contains(index1)){
                        int pos1 = indexList.indexOf(index1);
                        seq = br.readLine();
                        if(barcodeLists[pos1].size()!=0){
                            currentBarcode = seq.substring(1,8);
                            if(barcodeLists[pos1].contains(currentBarcode)){
                                int pos2 = barcodeLists[pos1].indexOf(currentBarcode);
                                bws[pos1][pos2].write(temp);bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(seq);bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(br.readLine());bws[pos1][pos2].newLine();
                                bws[pos1][pos2].write(br.readLine());bws[pos1][pos2].newLine();
                            }else{
                                br.readLine();br.readLine();
                            }

                        }else{
                            bw[pos1].write(temp);bw[pos1].newLine();
                            bw[pos1].write(seq);bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                        }
                    }else{
                        index2=temp.split(" ")[1].split(":")[3].substring(1);
                        if(indexList.contains(index2)){
                            seq=br.readLine();
                            int pos1 = indexList.indexOf(index2);
                            bw[pos1].write(temp);bw[pos1].newLine();
                            bw[pos1].write(seq);bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                            bw[pos1].write(br.readLine());bw[pos1].newLine();
                        }else{
                            br.readLine();br.readLine();br.readLine();
                        }
                    }
                }
                br.close();
                for(int i=0;i<indexList.size();i++){
                    for(int j=0;j<barcodeLists[i].size();j++){
                        bws[i][j].flush();bws[i][j].close();
                    }
                    if(barcodeLists[i].size()==0){
                        bw[i].flush();bw[i].close();
                    }
                }
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
    private void parseFqByIndex () {
        long startTimePoint = System.nanoTime();
        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq.txt");
//       RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030.txt");
        List<String> indexList = new ArrayList() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            barcodeStrain.put(t.getCell(i, 2), t.getCell(i, 0));
            indexList.add(t.getCell(i, 2));
        }
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/Novoseq";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, "R2_001.fq.gz");
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(f -> {
            try {
                BufferedWriter [] bw = new BufferedWriter[indexList.size()];
                for (int i = 0; i < indexList.size(); i++) {
//                   bw[i]=IOUtils.getTextWriter(new File(outputDirS, f.getName()+"_"+barcodeMethod.get(barcodeList.get(i))+".fq").getAbsolutePath());
                    bw[i]=IOUtils.getTextWriter(new File(outputDirS,barcodeStrain.get(indexList.get(i))+".fq").getAbsolutePath());
                }
                BufferedReader br = null;
//                BufferedReader br1 = null;
                if (f.getAbsolutePath().endsWith(".fq.gz")) {
                    br = IOUtils.getTextGzipReader(f.getAbsolutePath());
                }
                else {
                    br = IOUtils.getTextReader(f.getAbsolutePath());
                }
                int index1 = -1 ;String index2 = null ;
                String temp = null;String temp1 = null;
                String seq = null;String seq1 = null;
                String currentIndex = null;
                while((temp = br.readLine())!=null){
                    seq=br.readLine();//seq1=br1.readLine();temp1=br1.readLine();
                    currentIndex=temp.split(":")[9];
                    if(barcodeStrain.get(currentIndex)!=null){
                        index1=indexList.indexOf(currentIndex);
                        bw[index1].write(temp);bw[index1].newLine();
                        bw[index1].write(seq);bw[index1].newLine();
                        bw[index1].write(br.readLine());bw[index1].newLine();
                        bw[index1].write(br.readLine());bw[index1].newLine();
                    }else{
                        br.readLine();br.readLine();
                    }
                }
                for(int i=0;i<indexList.size();i++){
                    bw[i].flush();bw[i].close();
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
    private void parseFqByIndexAndBarcode () {
        long startTimePoint = System.nanoTime();
        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq20181030.txt");
//       RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/RNA-seq20181030.txt");
        List<String> indexList = new ArrayList() ;
        List<String> barcodeList = new ArrayList() ;
//        for (int i = 0; i < t.getRowNumber(); i++) {
//           barcodeMethod.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
//           barcodeList.add(t.getCell(i, 1).substring(1));
//        }
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0).split("-")[0]);
            barcodeStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            barcodeList.add(t.getCell(i, 1).substring(1));
        }
        String inputDirS = "/data1/home/junxu/wheat/unsplited";
        String outputDirS = "/data1/home/junxu/wheat/doubleAll";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, "_001.fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"-R1_001.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"-R2_001.fq.gz").getAbsolutePath();
                BufferedWriter [] bw = new BufferedWriter[barcodeList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[barcodeList.size()];
                for (int i = 0; i < barcodeList.size(); i++) {
                    bw[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R1.fq").getAbsolutePath());
                    bw1[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R2.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextReader(infile1);
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
    private void parseFqByBarcodePlateseq () {
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/RNA-seq20181030.txt");
//        RowTable<String> t = new RowTable<>("/data1/home/junxu/BRBUMI/UMItest/withoutUMI/IndexBarcode.txt");
//        RowTable<String> t = new RowTable<>("/data1/home/junxu/wheat/rnaseq20181204-ERCC/ERCC.txt");
        RowTable<String> t = new RowTable<>("/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/Rawdata/20191226ampm/IndexBarcode.txt");
        List<String> barcodeList = new ArrayList() ;
        List<String> strainList = new ArrayList() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            barcodeStrain.put(t.getCell(i, 2), t.getCell(i, 0)+"-"+t.getCell(i, 3));
            barcodeList.add(t.getCell(i, 2));
            strainList.add(t.getCell(i, 0)+"-"+t.getCell(i, 3));
        }
        int [] count =new int [strainList.size()];
//        String inputDirS = "/Users/xujun/Desktop/test-data/ANNO_ANCBJ170529_PM-ANCBJ170529-13_2020-01-19/Rawdata/20190427cs10am";
//        String outputDirS = "/Users/xujun/Desktop/test-data/ANNO_ANCBJ170529_PM-ANCBJ170529-13_2020-01-19/amResult";
//        String inputDirS = "/data1/home/junxu/BRBUMI/UMItest/withoutUMI/Rawdata/20191226ampm";
//        String outputDirS = "/data1/home/junxu/BRBUMI/UMItest/withoutUMI/result";
        String inputDirS = "/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/Rawdata/20191226ampm/";
        String outputDirS = "/Users/xujun/Desktop/test-data/ANNO_ANGBJ170529_PM-ANCBJ170529-11_2-2020-01-11/Rawdata/20191226ampm/";
        File[] fs = new File(inputDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "_R1_001.fastq.gz");
        fs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        List<File> fList = Arrays.asList(fs);
        nameSet.stream().forEach(p -> {
            try {
                String infile1 = new File (inputDirS, p+"_R1.fq.gz").getAbsolutePath();
                String infile2 = new File (inputDirS, p+"_R2.fq.gz").getAbsolutePath();
                BufferedWriter [] bwR1 = new BufferedWriter[barcodeList.size()];
                BufferedWriter [] bwR2 = new BufferedWriter[barcodeList.size()];
                BufferedWriter bw1 = IOUtils.getTextWriter(new File (outputDirS,"ratio.txt").getAbsolutePath());
                BufferedWriter bw2 = IOUtils.getTextWriter(new File (outputDirS,"error.txt").getAbsolutePath());
                for (int i = 0; i < barcodeList.size(); i++) {
                    bwR1[i]=IOUtils.getTextGzipWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R1.fq.gz").getAbsolutePath());
                    bwR2[i]=IOUtils.getTextGzipWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i))+"_R2.fq.gz").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);

                int index1 = -1 ;int index2 = -1 ;
                String temp = null;String temp1 = null;
                String seq = null;String seq1 = null;
                String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    seq=br.readLine();temp1=br1.readLine();seq1=br1.readLine();
                    currentBarcode=seq.substring(0,8);
                    if(barcodeStrain.get(currentBarcode)!=null){
                        index1=barcodeList.indexOf(currentBarcode);
                        index2=strainList.indexOf(barcodeStrain.get(currentBarcode));
                        count[index2]++;
//                        if(String.valueOf(barcodeMethod.get(currentBarcode)).contains("SiPAS")){
                        bwR1[index1].write(temp);bwR1[index1].newLine();
                        bwR1[index1].write(seq);bwR1[index1].newLine();
                        bwR1[index1].write(br.readLine());bwR1[index1].newLine();
                        bwR1[index1].write(br.readLine());bwR1[index1].newLine();
                        bwR2[index1].write(temp1);bwR2[index1].newLine();
                        bwR2[index1].write(seq1);bwR2[index1].newLine();
                        bwR2[index1].write(br1.readLine());bwR2[index1].newLine();
                        bwR2[index1].write(br1.readLine());bwR2[index1].newLine();

//                        }else{
//                            br1.readLine();br1.readLine();
//                        }
//                        br.readLine();br.readLine();
//                        br1.readLine();br1.readLine();
//                        bw[index1].write(temp);bw[index1].newLine();
//                        bw[index1].write(seq);bw[index1].newLine();
//                        bw[index1].write(br.readLine());bw[index1].newLine();
//                        bw[index1].write(br.readLine());bw[index1].newLine();
                    }else{
                        bw2.write(temp);bw2.newLine();
                        bw2.write(seq);bw2.newLine();
                        bw2.write(br.readLine());bw2.newLine();
                        bw2.write(br.readLine());bw2.newLine();
                        br1.readLine();br1.readLine();
                    }
                }
                for(int i=0;i<count.length;i++){
                    bw1.write(strainList.get(i)+"\t"+count[i]);
                    bw1.newLine();
                }
                for(int i=0;i<barcodeList.size();i++){
                    bwR1[i].flush();bwR1[i].close();
                    bwR2[i].flush();bwR2[i].close();
                }
                bw1.flush();bw1.close();
                bw2.flush();bw2.close();
                br.close();br1.close();
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
}

