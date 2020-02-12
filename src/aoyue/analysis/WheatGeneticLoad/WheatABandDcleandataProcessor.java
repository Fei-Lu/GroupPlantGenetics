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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class WheatABandDcleandataProcessor {

    public WheatABandDcleandataProcessor() {
        //this.getIDseqsubsamples();
        //this.getBWAscript();
        //this.getRmdupscript();
        //this.splitRmdupScript();
        //this.splitBwaScript();
        
        //this.getBWAscriptfromSuppleData();
        //this.getRmdupscriptfromSuppleData();
        //this.mergeBamfromSuppleData();
        
        //this.mergeBamScript();
        
        
        
        
        
        
        
    }
    
    /**
     * 将去处重复的bam文件进行合并
     */
    public void mergeBamScript(){
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_S106.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/003_mergeBam/listBams/"; //在HPC中/002_bamfile/ 文件夹下放置
        String scriptS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/003_mergeBam/001_script_mergeBam/001_allMergeBam.sh";
        
        String bamParentDirS = "/data1/home/aoyue/vmap2_AB_project/002_bamfile/";
        String mergeBamDirS = "/data1/home/aoyue/vmap2_AB_project/003_mergefile/";
        
        RowTable<String> t = new RowTable<>(dbfileS);
        Set<String> s = new HashSet<String>(t.getColumn(0));
        String[] IDs = s.toArray(new String [s.size()]);
        Arrays.sort(IDs);
        
        HashMap<String,String> hm = new HashMap<>(); //ID和倍性的映射关系
        
        /*************************得到id-samples列表，并写出文件到 outfileDirS 文件夹*****************************************/
        try{
            BufferedWriter[] bw = new BufferedWriter[IDs.length];
            for(int i=0;i<IDs.length;i++){
                String dbID = IDs[i];
                String outfileS = new File(outfileDirS,dbID + ".listbams.txt").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileS);
                BufferedReader br = IOUtils.getTextReader(dbfileS);
                String temp = br.readLine();
                while((temp=br.readLine()) != null){
                    String ID = PStringUtils.fastSplit(temp).get(0);
                    String ploidy = PStringUtils.fastSplit(temp).get(1);
                    String sample = PStringUtils.fastSplit(temp).get(2);
                    hm.put(ID, ploidy);
                    if(ID.equals(dbID)){
                        String bampathS;
                        if(ploidy.equals("AB")){
                            bampathS = bamParentDirS + "AB/" + sample + ".rmdup.bam";
                        }
                        else{
                            bampathS = bamParentDirS + "D/" + sample + ".rmdup.bam";
                        }
                        bw[i].write(bampathS);bw[i].newLine();
                    }
                }
                bw[i].flush(); bw[i].close();br.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        
        /*************************写出合并bam文件的脚本*****************************************/
        try{
            //samtools merge -f -b AT08079A.listbams.txt /data4/home/aoyue/vmap2/ab/003_mergefile/AB/AT08079A.rmdup.bam &
            BufferedWriter bw = IOUtils.getTextWriter(scriptS);
            for(int i=0;i<IDs.length;i++){
                String mergeBamS;
                if(hm.get(IDs[i]).equals("AB")){
                    mergeBamS = mergeBamDirS + "AB/" + IDs[i] + ".rmdup.bam";
                }
                else{
                    mergeBamS = mergeBamDirS + "D/" + IDs[i] + ".rmdup.bam";
                }
                String listbamfileS = bamParentDirS + "listBams/" + IDs[i] + ".listbams.txt";
                bw.write("samtools merge -f -b " + listbamfileS +" " + mergeBamS + " && samtools index " + mergeBamS + " &" );
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * 将补充测序的bam数据进行合并，文件名为测序ID名，文件内包含所有lanes的bam文件绝对路径;
     * ID_BGI	Ploidy	SeqSamples
    AT08079A	AB	190621_I5_V300018662_L2_WHEtkyRAAAA-501
    AT08079A	AB	190621_I5_V300018662_L2_WHEtkyRAAAA-502
    AT08079A	AB	190621_I5_V300018662_L2_WHEtkyRAAAA-503
    AT08079A	AB	190621_I5_V300018662_L2_WHEtkyRAAAA-504
    AT08079A	AB	190621_I5_V300018662_L2_WHEtkyRAAAA-505
     * 
     * 准备文件：*.listbams.txt
     * 
     */
    public void mergeBamfromSuppleData(){
        String dbfileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T15_S5.txt";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/003_mergeBam/listBams/"; //在集群上和rmdup合并前的文件在一个文件夹中放置，即/002_bamfile/ 文件夹下放置
        String scriptS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/003_mergeBam/001_script_mergeBam/001_mergeBam.sh";
        
        String bamParentDirS = "/data4/home/aoyue/vmap2/ab/002_bamfile/";
        String mergeBamDirS = "/data4/home/aoyue/vmap2/ab/003_mergefile/";
        
        RowTable<String> t = new RowTable<>(dbfileS);
        Set<String> s = new HashSet<String>(t.getColumn(0));
        String[] IDs = s.toArray(new String [s.size()]);
        Arrays.sort(IDs);
        
        HashMap<String,String> hm = new HashMap<>(); //ID和倍性的映射关系
        
        /*************************得到id-samples列表，并写出文件到 outfileDirS 文件夹*****************************************/
        try{
            BufferedWriter[] bw = new BufferedWriter[IDs.length];
            for(int i=0;i<IDs.length;i++){
                String dbID = IDs[i];
                String outfileS = new File(outfileDirS,dbID + "_1.listbams.txt").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileS);
                BufferedReader br = IOUtils.getTextReader(dbfileS);
                String temp = br.readLine();
                while((temp=br.readLine()) != null){
                    String ID = PStringUtils.fastSplit(temp).get(0);
                    String ploidy = PStringUtils.fastSplit(temp).get(1);
                    String sample = PStringUtils.fastSplit(temp).get(2);
                    hm.put(ID, ploidy);
                    if(ID.equals(dbID)){
                        String bampathS;
                        if(ploidy.equals("AB")){
                            bampathS = bamParentDirS + "AB/" + sample + ".rmdup.bam";
                        }
                        else{
                            bampathS = bamParentDirS + "D/" + sample + ".rmdup.bam";
                        }
                        bw[i].write(bampathS);bw[i].newLine();
                    }
                }
                bw[i].flush(); bw[i].close();br.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        
        /*************************写出合并bam文件的脚本*****************************************/
        try{
            //samtools merge -f -b AT08079A.listbams.txt /data4/home/aoyue/vmap2/ab/003_mergefile/AB/AT08079A.rmdup.bam &
            BufferedWriter bw = IOUtils.getTextWriter(scriptS);
            for(int i=0;i<IDs.length;i++){
                String mergeBamS;
                if(hm.get(IDs[i]).equals("AB")){
                    mergeBamS = mergeBamDirS + "AB/" + IDs[i] + "_1.rmdup.bam";
                }
                else{
                    mergeBamS = mergeBamDirS + "D/" + IDs[i] + "_1.rmdup.bam";
                }
                String listbamfileS = bamParentDirS + "listBams/" + IDs[i] + "_1.listbams.txt";
                bw.write("samtools merge -f -b " + listbamfileS +" " + mergeBamS + " && samtools index " + mergeBamS + " &" );
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Goal:将32个样品进行bam文件的排序去重复，并建立索引
     * 进行name排序，再加标签，然后坐标排序，最后去重复；
     */
    
    public void getRmdupscriptfromSuppleData(){
        String namesortpath=null; 
        String pebamPath=null;
        int core = 4; 
        //String memory = "4G";
        String memory = "8G";
        
        /********** Path in script ***********/
        String namesortpathAB = "/data4/home/aoyue/vmap2/ab/002_bamfile/AB/";
        String namesortpathD = "/data4/home/aoyue/vmap2/ab/002_bamfile/D/";
        String pebampathAB = "/data4/home/aoyue/vmap2/ab/002_bamfile/AB/";
        String pebampathD = "/data4/home/aoyue/vmap2/ab/002_bamfile/D/";
        
        try{
            String script = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/004_rmdup_suppleData.sh";
            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T15_S5.txt";
            BufferedWriter bw = IOUtils.getTextWriter(script);
            RowTable<String> t = new RowTable<>(db);
            for(int i =0; i <t.getRowNumber(); i++){
                String seqid = t.getCell(i, 0); String ploidy = t.getCell(i, 1); 
                String sample = t.getCell(i, 2);
                if(ploidy.equals("AB")){
                    namesortpath = namesortpathAB;
                    pebamPath = pebampathAB;
                }
                else{
                    namesortpath = namesortpathD;
                    pebamPath = pebampathD;
                }
                /****bampath********/
                String pefileS = new File(pebamPath,sample + ".pe.bam").getAbsolutePath();
                String namesorfileS = new File(namesortpath,sample + ".namesort.bam").getAbsolutePath();
                String fixmatefileS = new File(namesortpath,sample + ".fixmate.bam").getAbsolutePath();
                String fixmateposfileS = new File(namesortpath,sample + ".fixmate.pos.bam").getAbsolutePath();
                String rmdupfileS = new File(namesortpath,sample + ".rmdup.bam").getAbsolutePath();
                //samtools sort -n -m 20G -o /data2/aoyue/output/bamfile/GRC-L1.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/GRC-L1.pe.bam
                // && echo '** samtools namesort done **' &&
                //  "samtools fixmate -m /data2/aoyue/output/bamfile/GRC-L1.namesort.bam /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.pe.bam &&
                
                //  samtools sort -m 20G -o /data2/aoyue/output/bamfile/GRC-L1.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.namesort.bam &&
                //  samtools markdup -r /data2/aoyue/output/bamfile/GRC-L1.fixmate.pos.bam /data2/aoyue/output/bamfile/GRC-L1.rmdup.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam &&
                //  samtools index /data2/aoyue/output/bamfile/GRC-L1.rmdup.bam
                StringBuilder sb = new StringBuilder("samtools sort -n -m ");
                sb.append(memory).append(" -o ").append(namesorfileS).append(" -O bam -@ ").append(core).append(" ").append(pefileS);
                sb.append(" && echo '** samtools namesort done **' && ");
                sb.append("samtools fixmate -m ").append(namesorfileS).append(" ").append(fixmatefileS);
                sb.append(" && rm -f ").append(pefileS).append(" && ");
                
                sb.append("samtools sort -m ").append(memory).append(" -o ").append(fixmateposfileS).append(" -O bam -@ ").append(core).append(" ").append(fixmatefileS);
                sb.append(" && rm -f ").append(namesorfileS).append(" && ");
                sb.append("samtools markdup -r ").append(fixmateposfileS).append(" ").append(rmdupfileS);
                sb.append(" && rm -f ").append(fixmatefileS).append(" && ");
                sb.append("samtools index ").append(rmdupfileS).append(" &");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Goal:生成BWA脚本文件; 在204集群上运行剩余补测数据
     * 根据 /Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDlist.txt 进行AB 还是D材料的判断
     */
    public void getBWAscriptfromSuppleData(){
        /********** Path in script ***********/
        //bwa mem -t 8 -R '@RG\tID:AT08079A\tPL:illumina\tSM:AT08079A\tLB:AT08079A' /data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz 190520_I4_V300021171_L1_WHEtkyRAAAA-501_1.clean.fq.gz 190520_I4_V300021171_L1_WHEtkyRAAAA-501_2.clean.fq.gz | samtools view -S -b - > /data2/aoyue/vmap2_AB_project/002_bamfile/AB/190520_I4_V300021171_L1_WHEtkyRAAAA-501.pe.bam && echo '** bwa mapping done **' & 
        String refAB = "/data1/publicData/wheat/reference/v1.0/AB/bwaLib/ab_iwgscV1.fa.gz";
        String refD = "/data1/publicData/wheat/reference/v1.0/D/bwaLib/d_iwgscV1.fa.gz";
        String rawdataDirS = "/data4/home/aoyue/vmap2/ab/001_rawData/a_supplement_rawData/";
        String bampathAB = "/data4/home/aoyue/vmap2/ab/002_bamfile/AB/";
        String bampathD = "/data4/home/aoyue/vmap2/ab/002_bamfile/D/";
        try{
            String scriptS15 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/4T_15/001_bwa.S15.sh";
            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T15_S5.txt";
            BufferedWriter bw = IOUtils.getTextWriter(scriptS15);
            
            RowTable<String> t = new RowTable<>(db);
            for(int i =0; i <t.getRowNumber(); i++){
                String id = t.getCell(i, 0); String ploidy = t.getCell(i, 1); 
                String sample = t.getCell(i, 2);
                StringBuilder sb = new StringBuilder("bwa mem -t 4 -R ");
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
                sb.append("echo '** bwa mapping done **' &");
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
    
    public void splitBwaScript(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/bwa_20190705needRERUN.sh";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/splitScript";
        String shfileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/001_bwaScript/bwa_20190630needRERUN.sh";

        try{
        String[] outfileS= new String[20];
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter[] bw = new BufferedWriter[20];
        for(int i=0;i<outfileS.length;i++){
            String num = PStringUtils.getNDigitNumber(3, i+1);
            outfileS[i]=new File(outfileDirS,"bwa_20190706needRERUN_" + num + ".sh").getAbsolutePath();
            bw[i] = IOUtils.getTextWriter(outfileS[i]);
            String temp;
            for(int j=0; j<6;j++){
                if((temp = br.readLine()) != null){
                    bw[i].write(temp);
                    bw[i].newLine();
                }
            }
            bw[i].flush();bw[i].close();
        }
        br.close();
    }
    catch(Exception e){
        e.printStackTrace();
        System.exit(1);
    }

    try{
        File[] fs = new File(outfileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".sh");
        Arrays.sort(fs);
        BufferedWriter bw = IOUtils.getTextWriter(shfileS);
        for(int i=0; i<fs.length;i++){
            bw.write("sh " + fs[i].getName() + " &");
            bw.newLine();
        }
        bw.flush();bw.close();
    }
    catch(Exception e){
        e.printStackTrace();
        System.exit(1);
    }
            
        }
    
    
    
    /**
     * 将006_ABandD/000_cleandata/002_samtoolsRmdupScript/002_allRmdup_remove2test.sh中的926个命令拆分成多个脚本，在HPC上同时运行；
     * 先处理AB genome的脚本，合计639个文件， 160个sh脚本, 159（4CMD）+ 1（3CMD）；
     * 再处理D genme的脚本，合计287个文件，18个sh脚本，17（16CMD）+1（15CMD）；
     * 
     */
    public void splitRmdupScript(){
        String infileABS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/002_allRmdup_AB_20190709.sh";
        //String infileDS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/003_rmdup_D_remove2test.sh";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/splitScript_rmdup_20190709/";
        String shfileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/sh_splitScript.sh";
        
        /********** Path in script ***********/
        try{
            String[] outfileAB= new String[80];
            BufferedReader br = IOUtils.getTextReader(infileABS);
            BufferedWriter[] bw = new BufferedWriter[80];
            for(int i=0;i<outfileAB.length;i++){
                String num = PStringUtils.getNDigitNumber(3, i+1);
                outfileAB[i]=new File(outfileDirS,"rmdup_AB_" + num + ".sh").getAbsolutePath();
                bw[i] = IOUtils.getTextWriter(outfileAB[i]);
                String temp;
                for(int j=0; j<8;j++){
                    if((temp = br.readLine()) != null){
                        bw[i].write(temp);
                        bw[i].newLine();
                    }
                }
                bw[i].flush();bw[i].close();
            }
            br.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
//        try{
//            String[] outfileD= new String[18];
//            BufferedReader br = IOUtils.getTextReader(infileDS);
//            BufferedWriter[] bw = new BufferedWriter[18];
//            for(int i=0;i<outfileD.length;i++){
//                String num = PStringUtils.getNDigitNumber(3, i+1);
//                outfileD[i]=new File(outfileDirS,"rmdup_D_" + num + ".sh").getAbsolutePath();
//                bw[i] = IOUtils.getTextWriter(outfileD[i]);
//                String temp;
//                for(int j=0; j<16;j++){
//                    if((temp = br.readLine()) != null){
//                        bw[i].write(temp);
//                        bw[i].newLine();
//                    }
//                }
//                bw[i].flush();bw[i].close();
//            }
//            br.close();
//        }
//        catch(Exception e){
//            e.printStackTrace();
//            System.exit(1);
//        }
        
        try{
            File[] fs = new File(outfileDirS).listFiles();
            fs = IOUtils.listFilesEndsWith(fs, ".sh");
            Arrays.sort(fs);
            BufferedWriter bw = IOUtils.getTextWriter(shfileS);
            for(int i=0; i<fs.length;i++){
                bw.write("sh " + fs[i].getName() + " &");
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Goal:将928个样品进行bam文件的排序去重复，并建立索引
     * 分为3步：1.进行name排序，再加标签，然后坐标排序，最后去重复；
     */
    
    public void getRmdupscript(){
        String namesortpath=null; 
        String pebamPath=null;
        int core = 4; 
        String memory = "4G";
        
        /********** Path in script ***********/
        String namesortpathAB = "/data1/home/aoyue/vmap2_AB_project/002_bamfile/AB/";
        String namesortpathD = "/data1/home/aoyue/vmap2_AB_project/002_bamfile/D/";
        String pebampathAB = "/data2/aoyue/vmap2_AB_project/002_bamfile/AB/";
        String pebampathD = "/data2/aoyue/vmap2_AB_project/002_bamfile/D/";
        
        try{
            String script = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/002_samtoolsRmdupScript/001_allRmdup_20190709.sh";
            String db = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_S106.txt";
            BufferedWriter bw = IOUtils.getTextWriter(script);
            RowTable<String> t = new RowTable<>(db);
            for(int i =0; i <t.getRowNumber(); i++){
                String seqid = t.getCell(i, 0); String ploidy = t.getCell(i, 1); 
                String sample = t.getCell(i, 2);
                if(ploidy.equals("AB")){
                    namesortpath = namesortpathAB;
                    pebamPath = pebampathAB;
                }
                else{
                    namesortpath = namesortpathD;
                    pebamPath = pebampathD;
                }
                /****bampath********/
                String pefileS = new File(pebamPath,sample + ".pe.bam").getAbsolutePath();
                String namesorfileS = new File(namesortpath,sample + ".namesort.bam").getAbsolutePath();
                String fixmatefileS = new File(namesortpath,sample + ".fixmate.bam").getAbsolutePath();
                String fixmateposfileS = new File(namesortpath,sample + ".fixmate.pos.bam").getAbsolutePath();
                String rmdupfileS = new File(namesortpath,sample + ".rmdup.bam").getAbsolutePath();
                //samtools sort -n -m 20G -o /data2/aoyue/output/bamfile/GRC-L1.namesort.bam -O bam -@ 1 /data2/aoyue/output/bamfile/GRC-L1.pe.bam
                // && echo '** samtools namesort done **' &&
                //  "samtools fixmate -m /data2/aoyue/output/bamfile/GRC-L1.namesort.bam /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.pe.bam &&
                
                //  samtools sort -m 20G -o /data2/aoyue/output/bamfile/GRC-L1.fixmate.pos.bam -O bam -@ 1 /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.namesort.bam &&
                //  samtools markdup -r /data2/aoyue/output/bamfile/GRC-L1.fixmate.pos.bam /data2/aoyue/output/bamfile/GRC-L1.rmdup.bam 
                // && rm -f /data2/aoyue/output/bamfile/GRC-L1.fixmate.bam &&
                //  samtools index /data2/aoyue/output/bamfile/GRC-L1.rmdup.bam
                StringBuilder sb = new StringBuilder("samtools sort -n -m ");
                sb.append(memory).append(" -o ").append(namesorfileS).append(" -O bam -@ ").append(core).append(" ").append(pefileS);
                sb.append(" && echo '** samtools namesort done **' && ");
                sb.append("samtools fixmate -m ").append(namesorfileS).append(" ").append(fixmatefileS);
                sb.append(" && ");
                
                sb.append("samtools sort -m ").append(memory).append(" -o ").append(fixmateposfileS).append(" -O bam -@ ").append(core).append(" ").append(fixmatefileS);
                sb.append(" && rm -f ").append(namesorfileS).append(" && ");
                sb.append("samtools markdup -r ").append(fixmateposfileS).append(" ").append(rmdupfileS);
                sb.append(" && rm -f ").append(fixmatefileS).append(" && ");
                sb.append("samtools index ").append(rmdupfileS);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /**
     * Goal:生成BWA脚本文件;Map<String,List<T>> hashMap =new HashMap<String,List<T>>  List<T> list = new List<T>  hashMap.put("1",list)
     * 根据 /Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDlist.txt 进行AB 还是D材料的判断
     * 1.建立一个list类型的数组，list 大小是父目录的文件数，每个父目录包括8个文件名的数组元素；
     */
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
    
    public void getIDseqsubsamples(){
        //String infileDirS = "/Volumes/LuLab4T_17";
        String infileDirS = "/data4/home/aoyue/vmap2/ab/001_rawData/a_supplement_rawData";
        String outfileS15 = "/data4/home/aoyue/vmap2/ab/001_rawData/ABandDsubfileList_4T15_S5.txt";
        String outfileS16 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T16_S26.txt";
        String outfileS17 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T17_S40.txt";
        String outfileS18 = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDsubfileList_4T18_S40.txt";
        
        //String ploidyFileS = "/Users/Aoyue/project/wheatVMapII/006_ABandD/000_cleandata/000_dataCheck/ABandDlist.txt";
        String ploidyFileS = "/data4/home/aoyue/vmap2/ab/001_rawData/ABandDlist.txt";
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
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS17);
            //BufferedWriter bw = IOUtils.getTextWriter(outfileS18);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS15);
            
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
}
