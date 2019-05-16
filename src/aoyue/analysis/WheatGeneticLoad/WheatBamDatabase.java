/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class WheatBamDatabase {

    public WheatBamDatabase() {
//        this.VMapIAgenome(); //deprecated this method
//        this.VMapIABgenome(); //deprecated this method
//        this.VMapIABDgenome(); //deprecated this method
//        this.VMapIDgenome(); //deprecated this method
        
        this.VMapIgenomeMethod2();
        this.mergeVMapIbamDB(); 
        this.wheatJiaoABDgenome();
        //this.renameJiaoABD();
        this.mergeVMapIJiaobamDB();
        
        
    }
    
    public void mergeVMapIJiaobamDB(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String infileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/003_step3_VMapIandJiao/AllWheatBamDatabase.txt";
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            br = IOUtils.getTextReader(infileS);
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();

            br=IOUtils.getTextReader(infileS1);
            temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                bw.write(temp); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void renameJiaoABD(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Jiao_ABD/Jiao_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/renameJiao.sh";
        /**
         * SeqID	Accessions	DatabaseID
            AFG-L1	PI 80741	TW0055
            AFG-L2	PI 127098	TW0056
            AFG-L3	PI 366569	TW0057
         */
        try{
            ///data3/wgs/bam/ABD
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            String temp = br.readLine(); // header
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append("mv /data3/wgs/bam/ABD/").append(seqID).append(".rmdup.bam /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam").
                        append(" && ").append("mv /data3/wgs/bam/ABD/").append(seqID).append(".rmdup.bam.bai /data3/wgs/bam/ABD/").append(databaseID).append(".rmdup.bam.bai");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
            
            
            
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    

    /**
     * 修改焦老师的文件名，并且进行数据库的建立，再与VMapI进行合并；
     */
    public void wheatJiaoABDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/Jiao_ABD/Jiao_ABD_SeqIDAccession.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/002_step2_jiao/Jiao_ABD_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String taxa = PStringUtils.fastSplit(temp).get(3);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = PStringUtils.fastSplit(temp).get(2);
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append("ABD\t").append("/data3/wgs/bam/ABD/" + databaseID + ".rmdup.bam").append("\t").
                        append("350\t").append("Illumina HiSeq X\t").append("10X\tJiaoLab");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    /**
     * 通过同时读取VMapI_A.txt VMapI_AB.txt VMapI_ABD.txt VMapI_D.txt 4个文件，进行数据库的生成，一步到位。
     */
    public void VMapIgenomeMethod2(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2/";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, "VMapI_");
        for (int i=0; i<fs.length;i++){
            System.out.println(fs[i].getName());
        }
        System.out.println(fs.length);
        
        List<File> fList = Arrays.asList(fs);
        fList.parallelStream().forEach(file -> {
            String genome = file.getName().split("I_")[1].split(".t")[0];
            System.out.println(genome + "   genome is being executed");
            String infileS = file.getAbsolutePath();
            String outfileS = new File(outfileDirS,file.getName().split(".txt")[0] + "_bam_temp.txt" ).getAbsolutePath();
            String outfileS1 = new File(outfileDirS,file.getName().split(".txt")[0] + "_bam.txt" ).getAbsolutePath();
            try{
                BufferedReader br = IOUtils.getTextReader(infileS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                bw.write("DataBaseID\tTaxa\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
                String temp = br.readLine(); //表头
                while((temp = br.readLine()) != null){
                    String taxa = PStringUtils.fastSplit(temp).get(2);
                    String seqID = PStringUtils.fastSplit(temp).get(0);
                    String accessions = PStringUtils.fastSplit(temp).get(1);
                    String databaseID = null;
                    if(seqID.startsWith("B0")){
                        databaseID = seqID.replaceFirst("B0", "B00");
                    }
                    else if (seqID.startsWith("B1")) {
                        databaseID = seqID.replaceFirst("B1", "B01");
                    }
                    else if(seqID.startsWith("A0")){
                        databaseID = seqID.replaceFirst("A0", "A00");
                    }
                    else if(seqID.startsWith("TW0")){
                        databaseID = seqID.replaceFirst("TW0", "TW00");
                    }
                    else if(seqID.startsWith("D0")){
                        databaseID = seqID.replaceFirst("D0", "D00");
                    }
                    StringBuilder sb = new StringBuilder();
                    sb.append(databaseID).append("\t").append(taxa).append("\t").append(accessions).append("\t").append(genome).append("\t").append("/data3/wgs/bam/" + genome + "/" + 
                            databaseID + ".rmdup.bam").append("\t").
                            append("350\t").append("NovaSeq 6000\t").append("3X\tLuLab");
                    bw.write(sb.toString()); bw.newLine();
                }
                br.close();bw.flush();bw.close();
            }
            catch(Exception e){
                e.printStackTrace();
                System.exit(1);
            }
            
            RowTable t = new RowTable(outfileS);
            t.sortAsText(0);
            t.writeTextTable(outfileS1, IOFileFormat.Text);

            new File(outfileS).delete();
        });
        
    }
    
    /**
     * 将生成的A AB ABD 和 D bam文件数据库合并起来，成为一个文件。
     * 1.建立表头，写入头文件；
     * 2.建立数组，将文件名字按顺序读入；
     * 3.进行for循环，合并文件。
     */
    public void mergeVMapIbamDB(){
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1_method2";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/001_mergeVMapI/VMapI_allBam.txt";
        String[] genomes = {"A", "AB", "ABD", "D"};
        try{
            //只读入表头
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine());bw.newLine();
            br.close();
            //按文件顺序合并
            for(int i =0; i < genomes.length; i++){
                String infileS = new File(infileDirS,"VMapI_" + genomes[i] + "_bam.txt").getAbsolutePath();
                br = IOUtils.getTextReader(infileS);
                String temp = br.readLine(); //表头
                while((temp = br.readLine()) != null){
                    bw.write(temp); bw.newLine();
                }
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
    }
    
    public void VMapIDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_D.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_D_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_D_bam.txt";

        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("D0", "D00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("D\t").append("/data3/wgs/bam/D/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void VMapIABDgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_ABD.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_ABD_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_ABD_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("TW0", "TW00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("ABD\t").append("/data3/wgs/bam/ABD/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    public void VMapIABgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_AB.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_AB_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_AB_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID;
                if(seqID.startsWith("B0")){
                    databaseID = seqID.replaceFirst("B0", "B00");
                }
                else{
                    databaseID = seqID.replaceFirst("B1", "B01");
                }
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("AB\t").append("/data3/wgs/bam/AB/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
    
    /**
     * 1.get the name of each SeqID, and then rename it!
     * 2.Make table like DataBaseID	SeqID	Accessions	Genome-type	Bam-Path	Insert-size	Sequencing-platform	Coverage
     * 
     * 将文本读进去，建立HashMap
     */
    public void VMapIAgenome(){
        String infileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/source/VMapI_A.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam_temp.txt";
        String outfileS1 = "/Users/Aoyue/project/wheatVMapII/003_dataAnalysis/001_bamDatabase/000_step1/VMapI_A_bam.txt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("DataBaseID\tSeqID\tAccessions\tGenome-type\tBam-Path\tInsert-size(bp)\tSequencing-platform\tCoverage\tDataSource"); bw.newLine();
            String temp = br.readLine(); //表头
            while((temp = br.readLine()) != null){
                String seqID = PStringUtils.fastSplit(temp).get(0);
                String accessions = PStringUtils.fastSplit(temp).get(1);
                String databaseID = seqID.replaceFirst("A0", "A00");
                StringBuilder sb = new StringBuilder();
                sb.append(databaseID).append("\t").append(seqID).append("\t").append(accessions).append("\t").append("A\t").append("/data3/wgs/bam/A/").append("\t").
                        append("350\t").append("NovaSeq 6000\t").append("3X\tVMapI");
                bw.write(sb.toString()); bw.newLine();
            }
            br.close();bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        RowTable t = new RowTable(outfileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS1, IOFileFormat.Text);
        
        new File(outfileS).delete();
    }
}
