package daxing;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import pgl.app.hapScanner.HapScanner;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
        runHapScnner(args[0]);
    }

    public static void runHapScnner(String paramerterFileDir){
        List<File> files = IOTool.getFileListInDir(paramerterFileDir);
        HapScanner hapScanner;
        for (int i = 0; i < files.size(); i++) {
            hapScanner = new HapScanner(files.get(i).getAbsolutePath());
            long start = System.nanoTime();
            System.out.println(files.get(i).getName() + " completed in "+ Benchmark.getTimeSpanHours(start)+ " hours");
        }
    }

    public static void posAlleleFile(String posAlleleDir, String outDir_posAllele, String outDir_pos){
        List<File> files = IOTool.getFileListInDir(posAlleleDir);
        String[] outNames_posAllele = files.stream().map(File::getName).map(s -> s.replaceAll(".pos",".posAllele.txt")).toArray(String[]::new);
        String[] outNames_pos = files.stream().map(File::getName).map(s -> s.replaceAll(".pos",".pos.txt")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->{
            try (BufferedWriter bw_posAllele = IOTool.getWriter(new File(outDir_posAllele, outNames_posAllele[e]));
                 BufferedWriter bw_pos = IOTool.getWriter(new File(outDir_pos, outNames_pos[e]));
                 BufferedReader br = IOTool.getReader(files.get(e).getAbsolutePath())) {
                String header = "Chr\tPos\tRef\tAlt";
                bw_posAllele.write(header);
                bw_posAllele.newLine();
                String line;
                List<String> temp;
                StringBuilder sb =new StringBuilder();
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    sb.append(temp.get(0)).append("\t").append(temp.get(1));
                    bw_pos.write(sb.toString());
                    bw_pos.newLine();
                    bw_posAllele.write(line);
                    bw_posAllele.newLine();
                }
                bw_pos.flush();
                bw_posAllele.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void prepareParameterFile(String taxaRefBamFileDir, String posAlleleDir, String posDir, int[] chrID,
                                            String samtoolsPath, String threadsNum,
                                            String prepareParameterOutDir){
        List<File> taxaRefBamFiles = IOTool.getFileListInDirEndsWith(taxaRefBamFileDir, ".txt");
        List<File> posAlleleFiles = IOTool.getVisibleFileListInDir(posAlleleDir);
        List<File> posFiles =IOTool.getVisibleFileListInDir(posDir);
        // 490,000,000
        StringBuilder sb = new StringBuilder();
        BufferedWriter bw;
        double errorRate = 0.05;
        try {
            for (int i = 0; i < chrID.length; i++) {
                sb.setLength(0);
                sb.append("#App:\tHapScanner\n#Author:\tFei Lu\n#Email:\tflu@genetics.ac.cn; dr.lufei@gmail.com\n#Homepage" +
                        ":\thttps://plantgeneticslab.github.io/home/\n\n#HapScanner is used to perform genotyping of diploid " +
                        "species from whole genome sequenceing data, based on an existing genetic variation library.\n#To run" +
                        " and pipeline, the machine should have both Java 8 and samtools installed. The lib directory should " +
                        "stay with TIGER.jar in the same folder.\n#Command line example. java -Xmx100g -jar TIGER.jar -a " +
                        "HapScanner -p parameter_hapscanner.txt > log.txt &\n#To specify options, please edit the the " +
                        "parameters below. Also, please keep the order of parameters.\n\n#Parameter 1: The taxaRefBam file " +
                        "containing information of taxon and its corresponding refernece genome and bam files. The bam file " +
                        "should have .bai file in the same folder\n");
                sb.append(taxaRefBamFiles.get(i)).append("\n\n");
                sb.append("#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF " +
                        "format). The positions come from genetic variation library.\n");
                sb.append(
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C" +
                                ".\n#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n");
                sb.append(posAlleleFiles.get(i)).append("\n\n");
                sb.append("#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from" +
                        " haplotype library, which is used in mpileup.\n");
                sb.append(posFiles.get(i)).append("\n\n");
                sb.append("#Parameter 4: Chromosome or region on which genotyping will be performed (e.g. chromosome 1 is" +
                        " designated as 1. Region 1bp to 100000bp on chromosome 1 is 1:1,100000)\n");
                sb.append(chrID[i]).append(":1,490000000\n\n");
                sb.append("#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping " +
                        "are more likely to be genotyped as homozygote when the combined error rate is high.\n");
                sb.append(errorRate).append("\n\n");
                sb.append("#Parameter 6: The path of samtools\n");
                sb.append(samtoolsPath).append("\n\n");
                sb.append("#Parameter 7: Number of threads\n");
                sb.append(threadsNum).append("\n\n");
                sb.append("#Parameter 8: The directory of output\n");
                sb.append(prepareParameterOutDir);
                bw = IOTool.getWriter(new File(prepareParameterOutDir,
                        "chr"+PStringUtils.getNDigitNumber(3,chrID[i])+"_parameterFile.txt"));
                bw.write(sb.toString());
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void taxaRefBam(String gbsSampleDir, String refDir, String namParentsFile, String taxaRefBamOutDir){
        List<File> gbsSampleFiles = IOTool.getFileListInDirEndsWith(gbsSampleDir, ".bam");
        List<File> refChrsPath = IOTool.getFileListInDirEndsWith(refDir, ".fa.gz");
        Map<String, String> namParentsBamMap = RowTableTool.getMap(namParentsFile, 0, 5);
        String[] taxaRefBamOutFileNames = refChrsPath.stream().map(File::getName).map(s -> s.substring(0,6)).map(s -> s+"_taxaRefBam.txt").toArray(String[]::new);
        BufferedWriter bw;
        try {
            StringBuilder sb = new StringBuilder();
            List<String> keyList = new ArrayList<>(namParentsBamMap.keySet());
            Collections.sort(keyList);
            for (int i = 1; i < refChrsPath.size()-2; i++) {
                bw = IOTool.getWriter(new File(taxaRefBamOutDir, taxaRefBamOutFileNames[i]));
                sb.append("Taxa\tReference\tBams(A list of bams of the taxon, seperated by the delimiter of Tab)");
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < keyList.size(); j++) {
                    sb.setLength(0);
                    sb.append(keyList.get(j)).append("\t");
                    sb.append(refChrsPath.get(i)).append("\t");
                    sb.append(namParentsBamMap.get(keyList.get(j)));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                for (int j = 0; j < gbsSampleFiles.size(); j++) {
                    String sampleName = gbsSampleFiles.get(j).getName().substring(0,7);
                    sb.setLength(0);
                    sb.append(sampleName).append("\t");
                    sb.append(refChrsPath.get(i)).append("\t");
                    sb.append(gbsSampleFiles.get(j));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}