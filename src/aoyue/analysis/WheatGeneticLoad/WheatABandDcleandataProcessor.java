/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import format.table.RowTable;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class WheatABandDcleandataProcessor {

    public WheatABandDcleandataProcessor() {
        //this.getIDseqsubsamples();
        this.getBWAscript();
        
    }
    
    /**
     * Goal:生成BWA脚本文件;Map<String,List<T>> hashMap =new HashMap<String,List<T>>  List<T> list = new List<T>  hashMap.put("1",list)
     * 根据 /Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDlist.txt 进行AB 还是D材料的判断
     * 1.建立一个list类型的数组，list 大小是父目录的文件数，每个父目录包括8个文件名的数组元素；
     */
    public void getIDseqsubsamples(){
        String infileDirS = "/Volumes/LuLab4T_17";
        String outfileS16 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T16_S26.txt";
        String outfileS17 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T17_S40.txt";
        String outfileS18 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T18_S40.txt";
        String ploidyFileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDlist.txt";
        
        /******************** 建立ID和倍性的库 ****************************/
        HashMap<String,String> hmIDvsPloidy = new HashMap<>();
        RowTable<String> t = new RowTable<>(ploidyFileS);
        for(int i = 0 ; i < t.getRowNumber(); i++){
            String id = t.getCell(i, 0); String ploidy = t.getCell(i, 1);
            hmIDvsPloidy.put(id, ploidy);
        }
        
        /******************** 建立一个id对应的所有fq文件的数据库 *********************/
        try{
            File[] fs = new File(infileDirS).listFiles();
            
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS16);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS17);
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS18);
            
            bw.write("ID_BGI\tPloidy\tSeqSamples\n");
            for(int i =0 ; i < fs.length; i++){
                if(!(fs[i].getName().startsWith("A") || fs[i].getName().startsWith("W"))) continue;
                if(fs[i].isDirectory()){
                    String id = fs[i].getName();
                    File[] fsamples = new File(fs[i].getAbsolutePath()).listFiles();
                    for(int j = 0; j < fsamples.length; j++){
                        if(fsamples[j].isHidden()) continue;
                        if(fsamples[j].isDirectory()){
                            String sample = fsamples[j].getName();
                            bw.write( id + "\t" + hmIDvsPloidy.get(id) + "\t" + sample + "\n");
                        }
                    }
                }
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void getBWAscript(){
        /********** Path in script ***********/
        //bwa mem -t 8 -R '@RG\tID:AT08079A\tPL:illumina\tSM:AT08079A\tLB:AT08079A' /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 190520_I4_V300021171_L1_WHEtkyRAAAA-501_1.clean.fq.gz 190520_I4_V300021171_L1_WHEtkyRAAAA-501_2.clean.fq.gz | samtools view -S -b - > /data2/aoyue/vmap2_AB_project/002_bamfile/AB/190520_I4_V300021171_L1_WHEtkyRAAAA-501.pe.bam && echo '** bwa mapping done **' & 
        String refAB = "/data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz";
        String refD = "/data1/publicData/wheat/reference/v1.0/D/bwaLib/d_iwgscV1.fa.gz";
        String rawdataDirS = "/data2/aoyue/vmap2_AB_project/001_rawData/";
        String bampathAB = "/data2/aoyue/vmap2_AB_project/002_bamfile/AB/";
        String bampathD = "/data2/aoyue/vmap2_AB_project/002_bamfile/D/";
        try{
            /********** LuLab4T_16 ***********/
//            String scriptS16 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/4T_16/001_bwa.S16.sh";
//            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T16_S26.txt";
//            BufferedWriter bw = IOUtils.getTextWriter(scriptS16);
            /********** LuLab4T_17 ***********/
            String scriptS17 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/4T_17/001_bwa.S17.sh";
            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T17_S40.txt";
            BufferedWriter bw = IOUtils.getTextWriter(scriptS17);
            /********** LuLab4T_18 ***********/
//            String scriptS18 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/4T_18/001_bwa.S18.sh";
//            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T18_S40.txt";
//            BufferedWriter bw = IOUtils.getTextWriter(scriptS18);
            
            RowTable<String> t = new RowTable<>(db);
            for(int i =0; i <t.getRowNumber(); i++){
                String id = t.getCell(i, 0); String ploidy = t.getCell(i, 1); 
                String sample = t.getCell(i, 2);
                StringBuilder sb = new StringBuilder("bwa mem -t 8 -R ");
                sb.append("'@RG\\tID:" + id + "\\tPL:illumina\\tSM:" + id + "\\tLB:" + id + "' ");
                String bampath; String indexFileS;
                if(ploidy.equals("AB")){
                    indexFileS = refAB;
                    bampath = bampathAB;
                }
                else{
                    indexFileS = refD;
                    bampath = bampathD;
                }
                sb.append(indexFileS + " ");
                sb.append(rawdataDirS + id + "/" + sample + "/" + sample + "_1.clean.fq.gz" + " ");
                sb.append(rawdataDirS + id + "/" + sample + "/" + sample + "_2.clean.fq.gz");
                sb.append(" | ").append("samtools view -S -b - > ").append(new File(bampath,sample + ".pe.bam").
                        getAbsolutePath()).append(" && ");
                sb.append("echo '** bwa mapping done **'");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush(); bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
}
