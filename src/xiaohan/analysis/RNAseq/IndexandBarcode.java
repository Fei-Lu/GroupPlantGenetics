/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author yxh
 */
public class IndexandBarcode {
    HashMap barcodeStrain = new HashMap();
    HashMap indexStrain = new HashMap();
    
    public IndexandBarcode(){
        this.parseFqByIndex();
        //this.parseFqByIndexAndBarcode();
    }
    
    private void parseFqByIndexAndBarcode () {
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

    private void parseFqByIndex () {
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rnaseq/RNA-seq-coleoptile-20180315-1.txt");
        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/Book1.txt");

        List<String> IndexList = new ArrayList<String>() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            //barcodeStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0));
            IndexList.add(t.getCell(i, 1).substring(1));
        }
        String inputDirS = "/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/P101SC18112845-01-F004-WSW50-0129-weifen";
//        String inputDirS = "/data1/home/xiaohan/coleoptile/unsplitedtest";
        String outputDirS = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest";
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
                BufferedWriter [] bw = new BufferedWriter[IndexList.size()];
                BufferedWriter [] bw1 = new BufferedWriter[IndexList.size()];
                for (int i = 0; i < IndexList.size(); i++) {
                    bw[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(IndexList.get(i))+"_R1.fastq.fq").getAbsolutePath());
                    bw1[i]=IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(IndexList.get(i))+"_R2.fastq.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                int pos = -1 ;
                String temp = null;
                String seq = null;String index = null;
                String currentindex = null;
                //String currentBarcode = null;
                while((temp = br.readLine())!=null){
                    index=temp.split(":")[9].substring(1,6);
                    if(indexStrain.get(index)!=null){
                        currentindex = index;
                        if(barcodeStrain.get(currentindex)!=null){
                            pos=IndexList.indexOf(currentindex);
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
                for(int i=0;i<IndexList.size();i++){
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

    public static void main(String args[]) {
        new IndexandBarcode();
    }
}


