/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import pgl.infra.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class PopGenPara {

    public PopGenPara() {
        //this.getGroups();
        //this.mkFstCommands();
        //this.fstTable();
        //this.mergeChrFst();
        //this.mkPiCommand();
        //this.averagePi();
        //this.mkPiCommandbasedWindow();
        //this.mergeChrPi();
        //this.mkTajimaDCommandbasedWindow();
        //this.mergeChrTajimaD();
        this.mergePi10KwindowBased();
        
        
    }
    
 
    
    public PopGenPara(String b) {
        //this.fstTable(b);
        
    }
    
    public void mergePi10KwindowBased(){
        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/003_mergePi10KwindowBased/Merged_pi10K.txt";
//        File[] fs = new File(infileDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".txt");
//        for(int i =0; i < fs.length; i ++){
//            System.out.print(fs[i] + "\n");
//        }
        String[] infilename = {"Teosinte_pi10K.txt","Tropical_subtropical_pi10K.txt","Mixed_pi10K.txt",
            "China_specific_pi10K.txt","Non_stiff_stalk_pi10K.txt","Stiff_stalk_pi10K.txt"};
        String[] group = {"Teosinte","Tropical-subtropical","Mixed","China specific","Non-stiff stalk","Stiff stalk"};
        try{
            BufferedReader br = IOUtils.getTextReader(infileDirS + "China_specific_pi10K.txt");
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write(br.readLine() + "\tGeneticGruop\n" );
            for(int i = 0; i < group.length; i++){
                String infileS = infileDirS + infilename[i];
                br = IOUtils.getTextReader(infileS);
                String temp = br.readLine(); //read header
                while((temp = br.readLine()) != null){
                    bw.write(temp + "\t" + group[i] + "\n");
                }
                br.close();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }

        
        
        
    }
    
    public void mergeChrTajimaD(){
        
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/China_specific";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/China_specific_Tajima.D10K.txt";

//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Mixed";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Mixed_Tajima.D10K.txt";

//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Non_stiff_stalk";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Non_stiff_stalk_Tajima.D10K.txt";

//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Stiff_stalk";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Stiff_stalk_Tajima.D10K.txt";

        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Teosinte";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Teosinte_Tajima.D10K.txt";
        

        //String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Tropical_subtropical";
        //String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/000_tajimaD10KwindowBased/Tropical_subtropical_Tajima.D10K.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".Tajima.D");
        Arrays.sort(fs);
        try{
            double size = 0;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tBIN_START\tN_SNPS\tTajimaD");bw.newLine();
            for(int i=0;i<fs.length;i++){ 
                String nameprefix = fs[i].getName().split("chr")[0]; //China_specific_hmp321_agpv4_chr1.windowed.pi
                String chr = String.valueOf(i+1);
                String infileS = new File(infileDirS,nameprefix + "chr" + chr + ".Tajima.D").getAbsolutePath();
                size = new File(infileDirS,nameprefix + "chr" + chr + ".Tajima.D").length() + size;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = br.readLine();
                while((temp = br.readLine()) != null){
                        bw.write(temp);bw.newLine();
                }
                br.close();
            }
            size = size/1024/1024;
            DecimalFormat df = new DecimalFormat ("0.00");
            
            System.out.println(df.format(size));
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        
    }
    
    public void mkTajimaDCommandbasedWindow(){
        /*测试数据路径*/
        String groupFileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/000_groups";
        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/oriData/piTestData";
        String ScriptDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/testScriptBased10Kwindow";
        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/003_tajimaD/testTajimaD_Result_Based10Kwindow";
        /*HPC路径*/
//        String groupFileDirS = "/data2/aoyue/popGene/002_parameters/000_groups";
//        String infileDirS = "/data2/aoyue/maizeData/hmp321_agp4/";
//        String ScriptDirS = "/data2/aoyue/popGene/002_parameters/003_tajimaD/000_scriptBased10Kwindow";
//        String outfileDirS = "/data2/aoyue/popGene/002_parameters/003_tajimaD/000_tajimaDBased10Kwindow";
        
        File[] fs = new File(groupFileDirS).listFiles();
        for(int i=0; i<fs.length; i++){
            if(fs[i].isHidden()) fs[i].delete();
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(groupFileDirS).listFiles();
        List<String> scriptList = new ArrayList();
        try{
            for(int i=0; i<fs.length; i++){
                String scriptS = new File (ScriptDirS,fs[i].getName().replaceFirst(".txt", "")+"based10Kwindow.sh").getAbsolutePath();
                scriptList.add(scriptS);
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                File outfileS = new File(outfileDirS,fs[i].getName().replaceFirst(".txt", ""));
                outfileS.mkdirs();
                for(int j=0; j<10;j++){
                    String infileS = new File(infileDirS,"hmp321_agpv4_chr" + (j+1) + ".vcf.gz").getAbsolutePath();
                //vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --keep China_specific.txt --TajimaD 10000 --out hmp321_agpv4_chr1_china_specific
                    StringBuilder sb = new StringBuilder("vcftools --gzvcf ");
                    sb.append(infileS).append(" --keep ").append(fs[i]).append(" --TajimaD 10000 --out ").append(outfileS.getAbsolutePath()).append("/")
                            .append(fs[i].getName().replaceFirst(".txt", "")).append("_hmp321_agpv4_chr").append(j+1).append(" &");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();bw.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            String out = new File(ScriptDirS,"sh_tajimaDBased10Kwindow.sh").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(out);
            for(int i =0; i<scriptList.size();i++){
                bw.write("sh " + scriptList.get(i) + " &");
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mergeChrPi(){
        
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/China_specific";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/China_specific_pi10K.txt";
        
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Mixed";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Mixed_pi10K.txt";
        
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Non_stiff_stalk";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Non_stiff_stalk_pi10K.txt";

//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Stiff_stalk";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Stiff_stalk_pi10K.txt";

//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Teosinte";
//        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Teosinte_pi10K.txt";

        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Tropical_subtropical";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_pi10KwindowBased/Tropical_subtropical_pi10K.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".windowed.pi");
        Arrays.sort(fs);
        try{
            double size = 0;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tPI");bw.newLine();
            for(int i=0;i<fs.length;i++){ 
                String nameprefix = fs[i].getName().split("chr")[0]; //China_specific_hmp321_agpv4_chr1.windowed.pi
                String chr = String.valueOf(i+1);
                String infileS = new File(infileDirS,nameprefix + "chr" + chr + ".windowed.pi").getAbsolutePath();
                size = new File(infileDirS,nameprefix + "chr" + chr + ".windowed.pi").length() + size;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = br.readLine();
                while((temp = br.readLine()) != null){
                        bw.write(temp);bw.newLine();
                }
                br.close();
            }
            size = size/1024/1024;
            DecimalFormat df = new DecimalFormat ("0.00");
            
            System.out.println(df.format(size));
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        
    }
    
    public void mkPiCommandbasedWindow(){
        /*测试数据路径*/
//        String groupFileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/000_groups";
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/oriData/piTestData";
//        String ScriptDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/testScriptBased10Kwindow";
//        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/testPi_Result_Based10Kwindow";
        /*HPC路径*/
        String groupFileDirS = "/data2/aoyue/popGene/002_parameters/000_groups";
        String infileDirS = "/data2/aoyue/maizeData/hmp321_agp4/";
        String ScriptDirS = "/data2/aoyue/popGene/002_parameters/002_pi/scriptBased10Kwindow";
        String outfileDirS = "/data2/aoyue/popGene/002_parameters/002_pi/Pi_Result_Based10Kwindow";
        
        File[] fs = new File(groupFileDirS).listFiles();
        for(int i=0; i<fs.length; i++){
            if(fs[i].isHidden()) fs[i].delete();
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(groupFileDirS).listFiles();
        List<String> scriptList = new ArrayList();
        try{
            for(int i=0; i<fs.length; i++){
                String scriptS = new File (ScriptDirS,fs[i].getName().replaceFirst(".txt", "")+"based10Kwindow.sh").getAbsolutePath();
                scriptList.add(scriptS);
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                File outfileS = new File(outfileDirS,fs[i].getName().replaceFirst(".txt", ""));
                outfileS.mkdirs();
                for(int j=0; j<10;j++){
                    String infileS = new File(infileDirS,"hmp321_agpv4_chr" + (j+1) + ".vcf.gz").getAbsolutePath();
                //vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --keep China_specific.txt --window-pi 10000 --out hmp321_agpv4_chr1_china_specific
                    StringBuilder sb = new StringBuilder("vcftools --gzvcf ");
                    sb.append(infileS).append(" --keep ").append(fs[i]).append(" --window-pi 10000 --out ").append(outfileS.getAbsolutePath()).append("/")
                            .append(fs[i].getName().replaceFirst(".txt", "")).append("_hmp321_agpv4_chr").append(j+1).append(" &");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();bw.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            String out = new File(ScriptDirS,"sh_piBased10Kwindow.sh").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(out);
            for(int i =0; i<scriptList.size();i++){
                bw.write("sh " + scriptList.get(i) + " &");
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /*
    pi值有 -nan 0.02892等
    double[] pi 是放每个亚群中，每条染色体的pi平均值，TDoubleArrayList piList是放每条染色体的每个点的pi值。
    一个亚群一个亚群进行计算
    */
    public void averagePi () {
        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/000_piBasedSNP/";
        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/001_piTable/";
        
        //String group = "China_specific";
        //String group = "Stiff_stalk";
        //String group = "Non_stiff_stalk";
        //String group = "Teosinte";
        //String group = "Tropical_subtropical";
        String group = "Mixed";
        
        String outfileS = new File(outfileDirS,group + ".txt").getAbsolutePath();
        File[] fs = new File(infileDirS,group).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "gz");
        try{
            BufferedWriter bw =IOUtils.getTextWriter(outfileS);
            bw.write("group\tchr1\tchr2\tchr3\tchr4\tchr5\tchr6\tchr7\tchr8\tchr9\tchr10\ttotal");bw.newLine();
            bw.write(group);bw.write("\t");
            double[] pi = new double[fs.length];
            double sum = 0;
            StringBuilder sb = new StringBuilder();
            for(int i=0; i<fs.length; i++){ //China_specific_hmp321_agpv4_chr3.sites.pi.gz
                String infileS = new File(infileDirS,group).getAbsolutePath() + "/" + fs[i].getName().split("chr")[0] + "chr" + String.valueOf(i+1) + ".sites.pi.gz";
                BufferedReader br = IOUtils.getTextGzipReader(infileS);
                String temp = br.readLine();
                TDoubleArrayList piList = new TDoubleArrayList();
                while ((temp = br.readLine()) != null) {
                    String[] tem = temp.split("\t");
                    if (tem[2].startsWith("-n") || tem[2].startsWith("n")) continue;
                    piList.add(Double.valueOf(tem[2]));
                }
                br.close();
                DescriptiveStatistics d = new DescriptiveStatistics(piList.toArray());
                pi[i] = d.getMean();
                sum = sum+pi[i];
                sb.append(String.valueOf(pi[i])).append("\t");
                System.out.println("chr " + (i+1) + " is processing");
            }
            bw.write(sb.toString());bw.write(String.valueOf(sum));bw.newLine();
            bw.flush();bw.close();
            System.out.println(group + "\tfinished statistic");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    /*****
     1.进入group文件，遍历6个gruop
     * 2.对每个group进行的每条染色体进行pi计算；
     * 一共生成6*10=60个命令，每个group写入一个脚本；6个脚本.
     ****/
    public void mkPiCommand(){
        /*测试数据路径*/
//        String groupFileDirS = "/Users/Aoyue/Documents/Data/project/maizeGeneticLoad/005_popGen/002_parameters/000_groups/";
//        String infileDirS = "/Users/Aoyue/Documents/Data/project/maizeGeneticLoad/oriData/piTestData/";
//        String ScriptDirS = "/Users/Aoyue/Documents/Data/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/testScriptBasedSNP";
//        String outfileDirS = "/Users/Aoyue/Documents/Data/project/maizeGeneticLoad/005_popGen/002_parameters/002_pi/000_piBasedSNP";
        /*HPC路径*/
        String groupFileDirS = "/data2/aoyue/popGene/002_parameters/000_groups";
        String infileDirS = "/data2/aoyue/maizeData/hmp321_agp4/";
        String ScriptDirS = "/data2/aoyue/popGene/002_parameters/002_pi/scriptBasedSNP";
        String outfileDirS = "/data2/aoyue/popGene/002_parameters/002_pi/000_piBasedSNP";
        
        
        File[] fs = new File(groupFileDirS).listFiles();
        for(int i=0; i<fs.length; i++){
            if(fs[i].isHidden()) fs[i].delete();
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(groupFileDirS).listFiles();
        List<String> scriptList = new ArrayList();
        try{
            for(int i=0; i<fs.length; i++){
                String scriptS = new File (ScriptDirS,fs[i].getName().replaceFirst(".txt", "")+".sh").getAbsolutePath();
                scriptList.add(scriptS);
                BufferedWriter bw = IOUtils.getTextWriter(scriptS);
                File outfileS = new File(outfileDirS,fs[i].getName().replaceFirst(".txt", ""));
                outfileS.mkdirs();
                for(int j=0; j<10;j++){
                    String infileS = new File(infileDirS,"hmp321_agpv4_chr" + (j+1) + ".vcf.gz").getAbsolutePath();
                //vcftools --gzvcf hmp321_agpv4_chr1.vcf.gz --keep China_specific.txt --site-pi --out hmp321_agpv4_chr1_china_specific
                    StringBuilder sb = new StringBuilder("vcftools --gzvcf ");
                    sb.append(infileS).append(" --keep ").append(fs[i]).append(" --site-pi --out ").append(outfileS.getAbsolutePath()).append("/")
                            .append(fs[i].getName().replaceFirst(".txt", "")).append("_hmp321_agpv4_chr").append(j+1).append(" &");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();bw.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        
        try{
            String out = new File(ScriptDirS,"sh_pi.sh").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(out);
            for(int i =0; i<scriptList.size();i++){
                bw.write("sh " + scriptList.get(i) + " &");
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public void mergeChrFst(){
//        String infileDirS = "/Users/Aoyue/Documents/teoVSstiff";
//        String outfileS = "/Users/Aoyue/Documents/teoVSstiff.txt";
        
//        String infileDirS = "/Users/Aoyue/Documents/000_Non_stiff_stalkVSChina_specific";
//        String outfileS = "/Users/Aoyue/Documents/000_Non_stiff_stalkVSChina_specific.txt";
        
//        String infileDirS = "/Users/Aoyue/Documents/001_Non_stiff_stalkVSStiff_stalk";
//        String outfileS = "/Users/Aoyue/Documents/001_Non_stiff_stalkVSStiff_stalk.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/002_Non_stiff_stalkVSTropical_subtropical";
//        String outfileS = "/Users/Aoyue/Documents/002_Non_stiff_stalkVSTropical_subtropical.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/003_Stiff_stalkVSChina_specific";
//        String outfileS = "/Users/Aoyue/Documents/003_Stiff_stalkVSChina_specific.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/004_TeosinteVSChina_specific";
//        String outfileS = "/Users/Aoyue/Documents/004_TeosinteVSChina_specific.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/005_TeosinteVSNon_stiff_stalk";
//        String outfileS = "/Users/Aoyue/Documents/005_TeosinteVSNon_stiff_stalk.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/006_TeosinteVSStiff_stalk";
//        String outfileS = "/Users/Aoyue/Documents/006_TeosinteVSStiff_stalk.txt";
//        
//        String infileDirS = "/Users/Aoyue/Documents/007_TeosinteVSTropical_subtropical";
//        String outfileS = "/Users/Aoyue/Documents/007_TeosinteVSTropical_subtropical.txt";

//        String infileDirS = "/Users/Aoyue/Documents/008_Tropical_subtropicalVSChina_specific";
//        String outfileS = "/Users/Aoyue/Documents/008_Tropical_subtropicalVSChina_specific.txt";

        String infileDirS = "/Users/Aoyue/Documents/009_Tropical_subtropicalVSStiff_stalk";
        String outfileS = "/Users/Aoyue/Documents/009_Tropical_subtropicalVSStiff_stalk.txt";
        
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".weir.fst");
        Arrays.sort(fs);
        try{
            double size = 0;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tWEIGHTED_FST\tMEAN_FST");bw.newLine();
            for(int i=0;i<fs.length;i++){
                String nameprefix = fs[i].getName().split("chr")[0];
                String chr = PStringUtils.getNDigitNumber(3, i+1);
                String infileS = new File(infileDirS,nameprefix + "chr" + chr + "_fst.txt.windowed.weir.fst").getAbsolutePath();
                size = new File(infileDirS,nameprefix + "chr" + chr + "_fst.txt.windowed.weir.fst").length() + size;
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = br.readLine();
                while((temp = br.readLine()) != null){
                        bw.write(temp);bw.newLine();
                }
                br.close();
            }
            size = size/1024/1024;
            DecimalFormat df = new DecimalFormat ("0.00");
            
            System.out.println(df.format(size));
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        
    }
    
    public void fstTable (String b) {
        int chr = Integer.valueOf(b);
        /** 本机运行测试 **/
//        String groupDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/000_groups";
//        String fstDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/001_Fst/TestFst";
//        String outfileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/001_Fst/testFstTable.txt";
        /** HPC运行 **/
        String groupDirS = "/data2/aoyue/popGene/002_parameters/000_groups/";
        String fstDirS = "/data2/aoyue/popGene/002_parameters/001_Fst/";
        String outfileS = "/data2/aoyue/popGene/002_parameters/fstTable_chr" + chr +".txt";
        
        File[] fs = new File (groupDirS).listFiles();
        for(int i=0; i < fs.length; i++){
            if(fs[i].isHidden()) fs[i].delete();
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File (groupDirS).listFiles();
        String[] groups = new String[fs.length];
        for (int i = 0; i < groups.length; i++) {
            groups[i] = fs[i].getName().replaceFirst(".txt", "");
        }
        Arrays.sort(groups);
        double[][] fstValues = new double[groups.length][groups.length]; //建立二维数组
        for (int i = 0; i < groups.length-1; i++) {
            for (int j = i + 1; j < groups.length; j++) {
                String infileS = groups[i]+"VS"+groups[j]+"_chr"+PStringUtils.getNDigitNumber(3, chr)+"_fst.txt.weir.fst";
                infileS = new File (fstDirS, infileS).getAbsolutePath();
                if (!(new File (infileS).exists())) { //如果China_specificVSMixed_chr010_fst.txt.weir.fst 不存在,则是相反的顺序。
                    infileS = groups[j]+"VS"+groups[i]+"_chr"+PStringUtils.getNDigitNumber(3, chr)+"_fst.txt.weir.fst";
                    infileS = new File (fstDirS, infileS).getAbsolutePath();
                }
                try {
                    BufferedReader br = IOUtils.getTextReader(infileS);
                    /**
                     * CHROM	POS	WEIR_AND_COCKERHAM_FST
                        10	186842	0.0310016
                        10	1094857	-0.00302266
                        10	1304813	nan
                        10	1516109	0.314593
                        10	1583996	nan
                     */
                    String temp = br.readLine();
                    TDoubleArrayList vList = new TDoubleArrayList();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = PStringUtils.fastSplit(temp);
                        String t = l.get(2); //第二列为fst的值
                        if (t.startsWith("-n") || t.startsWith("n")) continue; //忽略那些没有n的位点
                        double v = Double.valueOf(t);
                        if (v < 0) v = 0;
                        vList.add(v);
                    }
                    br.close();
                    double[] v = vList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(v);
                    fstValues[i][j] = d.getMean();
                    System.out.println(infileS);
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        /*********** 建立相互关系的表格 *****************/
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder("Pupoluation");
            for (int i = 0; i < groups.length; i++) {
                sb.append("\t").append(groups[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < groups.length; i++) {
                sb = new StringBuilder(groups[i]);
                for (int j = 0; j < groups.length; j++) {
                    if (fstValues[i][j] == 0) {
                        if (i == j) {
                            sb.append("\tNA");
                        }
                        else {
                            sb.append("\t");
                        }
                    }
                    else {
                        sb.append("\t").append(fstValues[i][j]);
                    }
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
    
    /**
     * 这个方法有一个亮点，就是将写出的几十个脚本，统一写在一个脚本中，执行一下。
     */
    
    public void mkFstCommands () {
        /******************** 输入文件的准备 ***************************/
        /***测试数据输入*****/
        //String genotypeDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData/fstTestData";
        //String groupDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/000_groups";
        //String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/TestFst";
        // String shScriptDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/TestScript";
        
        //String outputDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/001_Fst/testFst_based10k";
        //String shScriptDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/001_Fst/testScript_based10k";

        /***HPC数据输入*****/
        String genotypeDirS = "/data2/aoyue/maizeData/hmp321_agp4";
        String groupDirS = "/data2/aoyue/popGene/002_parameters/000_groups";
        String outputDirS = "/data2/aoyue/popGene/002_parameters/001_Fst";
        String shScriptDirS = "/data2/aoyue/popGene/002_parameters/001_scripts";
        
        File[] genotypeFiles = new File(genotypeDirS).listFiles();
        genotypeFiles = IOUtils.listFilesEndsWith(genotypeFiles, "vcf.gz");
        File[] groupFiles = new File(groupDirS).listFiles();
        for(int i=0; i < groupFiles.length; i++){
            if(groupFiles[i].isHidden()) groupFiles[i].delete();
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        groupFiles = new File(groupDirS).listFiles();
        new File(outputDirS).mkdir();
        new File (shScriptDirS).mkdir();
        ArrayList<String> perlList = new ArrayList(); //在循环外建立perlList集合， 每个集合包含多个字符串，一个字符串代表一个文件。
        for (int i = 0; i < groupFiles.length-1; i++) {
            String group1 = groupFiles[i].getName().replace(".txt", "");
            for (int j = i+1; j < groupFiles.length; j++) {
                String group2 = groupFiles[j].getName().replace(".txt", "");
                try {
                    String outfileS = new File (shScriptDirS, group1+"VS"+group2+".sh").getAbsolutePath(); //写入的是命令，每个文件包含n条染色体的命令。
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    //vcftools --gzvcf test.vcf.gz --weir-fst-pop ../groups/Teosinte.txt --weir-fst-pop ../groups/Stiff_stalk.txt --weir-fst-pop ../groups/Non_stiff_stalk.txt --out out.txt
                    for (int k = 0; k < genotypeFiles.length; k++) {
                        int chr = Integer.valueOf(genotypeFiles[k].getName().replaceFirst("hmp321_agpv4_chr", "").replaceFirst(".vcf.gz", ""));
                        String output = new File (outputDirS, group1+"VS"+group2+"_chr"+PStringUtils.getNDigitNumber(3, chr)+"_fst.txt").getAbsolutePath();
                        StringBuilder sb = new StringBuilder("vcftools --gzvcf ");
                        sb.append(genotypeFiles[k].getAbsolutePath()).append(" --weir-fst-pop ").append(groupFiles[i]).append(" --weir-fst-pop ").append(groupFiles[j]);
                        sb.append(" --fst-window-size 10000").append(" --fst-window-step 2000");
                        //sb.append(" --fst-window-size 10000");
                        sb.append(" --out ").append(output).append(" &");
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    perlList.add(outfileS); //加入的是一个文件路径字符串
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
        
        /******** 将15个组合的脚本命令统一执行，有一个统一执行的脚本 sh_fst.sh *******/
        try {
            String outfileS = new File (shScriptDirS, "sh_fst.sh").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < perlList.size(); i++) {
                StringBuilder sb = new StringBuilder("sh ");
                sb.append(perlList.get(i))
                        .append(" &");;
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
    
    /*********** 将 Teosinte	Tropical-subtropical	Non-stiff stalk	Stiff stalk	Mixed	China specific 写入6个文件中，只有taxa名字 **********/
    public void getGroups(){
        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/000_group/geneticGroup.manual.txt";
        String outfileDirS = "/Users/Aoyue/Documents/maizeGeneticLoad/005_popGen/002_parameters/000_groups/";
        /**
         * Taxa	PC1	PC2	GroupIndex	GeneticGroup	CAU	CIMMYT-BGI	German	HapMap2	HapMap2_extra	282-2x	282-4x
         * 30-1	-0.032416	0.001835	3	Stiff stalk	1	0	0	0	0	0	0
         * 8-64	0.002665	0.053714	2	Non-stiff stalk	1	0	0	0	0	0	0
        */
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(4);
        Set<String> s = new HashSet(l);
        String[] groups = s.toArray(new String[s.size()]);
        Arrays.sort(groups);
        try{
            for(int i=0; i < groups.length; i++){
                String group = groups[i].replaceAll("-", "_").replaceAll(" ", "_");
                String outfileS = new File(outfileDirS,group+ ".txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                BufferedReader br = IOUtils.getTextReader(infileS);
                String temp = br.readLine();
                while((temp = br.readLine()) != null){
                    l = PStringUtils.fastSplit(temp);
                    String key = l.get(4);
                    String taxa = l.get(0);
                    if(key.equals(groups[i])){
                        bw.write(taxa);
                        bw.newLine();
                    }
                }
                bw.flush();bw.close();br.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    
}
