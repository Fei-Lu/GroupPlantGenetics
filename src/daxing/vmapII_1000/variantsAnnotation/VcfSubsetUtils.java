package daxing.vmapII_1000.variantsAnnotation;

import com.google.common.collect.Table;
import daxing.common.ArrayTool;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

/**
 * Utils of sampling methods for different purposes
 */
public class VcfSubsetUtils {

    /**
     * Sampling a fixed number of variants based on the distribution of variants on chromosomes
     * @param chrSNPNumberFile
     * Chr	SnpNumber
     * 1	9314623
     * 2	2120245
     * @param variantsNum
     */
    public static void sampleFixedNumVariants(String chrSNPNumberFile, String vcfInputDir, int variantsNum,
                                             String outDir, String variantsNumSummaryFile){
        System.out.println("start sample fixed number variants ...");
        long start0 = System.nanoTime();
        RowTableTool<String> table = new RowTableTool<>(chrSNPNumberFile);
        int[] snpNumArray=table.getColumnAsIntArray(1);
        double[] snpFreArray= ArrayTool.getElementPercent(snpNumArray);
        int[] subsetSnpNumArray=new int[snpFreArray.length];
        for (int i = 0; i < subsetSnpNumArray.length; i++) {
            subsetSnpNumArray[i]=(int)(variantsNum*snpFreArray[i]);
        }
        int[] intervalArray=new int[snpFreArray.length];
        for (int i = 0; i < intervalArray.length; i++) {
            intervalArray[i]=snpNumArray[i]/subsetSnpNumArray[i];
        }
        List<File> fileList = IOTool.getFileListInDirEndsWith(vcfInputDir, ".vcf.gz");
        String[] outNames=
                fileList.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","."+variantsNum/1000+"k.vcf.gz")).toArray(String[]::new);
        ConcurrentHashMap<Integer, Integer> chrSnpNumMap = new ConcurrentHashMap<>();
        IntStream.range(0, snpFreArray.length).parallel().forEach(e->{
            long start= System.nanoTime();
            int[] indexArray=IntStream.iterate(1, d -> d+intervalArray[e]).limit(subsetSnpNumArray[e]).toArray();
            TIntArrayList indexList=new TIntArrayList(indexArray);
            try (BufferedReader br = IOTool.getReader(fileList.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                int count=0;
                int cnt=0;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    count++;
                    if (!indexList.contains(count)) continue;
                    cnt++;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
                chrSnpNumMap.put(e+1, cnt);
                System.out.println(outNames[e]+" complicated in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
        writeVariantsNum(variantsNumSummaryFile, chrSnpNumMap);
        System.out.println("sample fixed number variants complicated in "+ Benchmark.getTimeSpanHours(start0)+ " hours");
    }

    private static void writeVariantsNum(String variantsNumSummaryFile, Map<Integer, Integer> chrSnpNumMap) {
        System.out.println("start writing variants number summary file...");
        long start = System.nanoTime();
        try (BufferedWriter bw = IOTool.getWriter(variantsNumSummaryFile)) {
            bw.write("Chr\tVariantsNum");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            Set<Integer> keySet=chrSnpNumMap.keySet();
            List<Integer> keyList= new ArrayList<>(keySet);
            Collections.sort(keyList);
            for (Integer i : keyList){
                sb.setLength(0);
                sb.append(i).append("\t").append(chrSnpNumMap.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("writing variants number summary file complicated in "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
    }

//    /**
//     * Extract the variants in gene region
//     * @param pgfFileGeneHC
//     * @param vcfInputDir
//     * @param outDir
//     */
//    public static void extractVariantsInGene(String pgfFileGeneHC, String vcfInputDir, String outDir,
//                                             String variantsNumSummaryFile){
//        System.out.println("start extract variants in gene region ...");
//        long start0 = System.nanoTime();
//        List<File> files = IOTool.getFileListInDirEndsWith(vcfInputDir, ".vcf.gz");
//        String[] outNames= files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".gene.vcf.gz")).toArray(String[]::new);
//        PGF pgf = new PGF(pgfFileGeneHC);
//        pgf.sortGeneByGeneRange();
//        ConcurrentHashMap<Integer, Integer> chrSnpNumMap = new ConcurrentHashMap<>();
//        IntStream.range(0, files.size()).parallel().forEach(e->{
//            long start= System.nanoTime();
//            try (BufferedReader br = IOTool.getReader(files.get(e));
//                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
//                String line, subLine;
//                List<String> temp;
//                int chrID, pos;
//                while ((line=br.readLine()).startsWith("##")){
//                    bw.write(line);
//                    bw.newLine();
//                }
//                bw.write(line);
//                bw.newLine();
//                int cnt=0;
//                while((line=br.readLine())!=null){
//                    subLine = line.substring(0,20);
//                    temp = PStringUtils.fastSplit(subLine);
//                    chrID = Integer.parseInt(temp.get(0));
//                    pos = Integer.parseInt(temp.get(1));
//                    int geneIndex = pgf.getGeneIndex(chrID, pos);
//                    if (geneIndex < 0) continue;
//                    cnt++;
//                    bw.write(line);
//                    bw.newLine();
//                }
//                bw.flush();
//                chrSnpNumMap.put(e+1, cnt);
//                System.out.println(outNames[e]+" complicated in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
//            } catch (IOException ioException) {
//                ioException.printStackTrace();
//            }
//        });
//        writeVariantsNum(variantsNumSummaryFile, chrSnpNumMap);
//        System.out.println("extract variants in gene region complicated in "+ Benchmark.getTimeSpanHours(start0)+ " hours");
//    }

    public static void extractVariantsWithAncestralState(String ancestralDir, String vcfDir, String outDir, String variantsNumSummaryFile){
        System.out.println("start extract variants with ancestral state ...");
        long start0 = System.nanoTime();
        List<File> ancestralFiles = IOTool.getFileListInDirEndsWith(ancestralDir, ".txt.gz");
        List<File> vcfFiles = IOTool.getFileListInDirEndsWith(vcfDir, ".vcf.gz");
        String[] outNames =
                vcfFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".ancestral.vcf.gz")).toArray(String[]::new);
        Map<Integer, Integer> chrSnpNumMap = new HashMap<>();
        for (int i = 0; i < vcfDir.length(); i++) {
            long start= System.nanoTime();
            Table<String, String, String> table = RowTableTool.getTable(ancestralFiles.get(i).toString(), 2);
            try (BufferedReader br = IOTool.getReader(vcfFiles.get(i));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[i]))) {
                String line, subLine;
                List<String> temp;
                int cnt=0;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    subLine = line.substring(0, 20);
                    temp = PStringUtils.fastSplit(subLine);
                    if (!table.contains(temp.get(0), temp.get(1))) continue;
                    cnt++;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
                chrSnpNumMap.put(i+1, cnt);
                System.out.println(outNames[i]+" complicated in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        }
//        ConcurrentHashMap<Integer, Integer> chrSnpNumMap = new ConcurrentHashMap<>();
//        IntStream.range(0, vcfFiles.size()).forEach(e->{
//            long start= System.nanoTime();
//            Table<String, String, String> table = RowTableTool.getTable(ancestralFiles.get(e).toString(), 2);
//            try (BufferedReader br = IOTool.getReader(vcfFiles.get(e));
//                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
//                String line, subLine;
//                List<String> temp;
//                int cnt=0;
//                while ((line=br.readLine()).startsWith("##")){
//                    bw.write(line);
//                    bw.newLine();
//                }
//                bw.write(line);
//                bw.newLine();
//                while ((line=br.readLine())!=null){
//                    subLine = line.substring(0, 20);
//                    temp = PStringUtils.fastSplit(subLine);
//                    if (!table.contains(temp.get(0), temp.get(1))) continue;
//                    cnt++;
//                    bw.write(line);
//                    bw.newLine();
//                }
//                bw.flush();
//                chrSnpNumMap.put(e+1, cnt);
//                System.out.println(outNames[e]+" complicated in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
//            } catch (IOException ioException) {
//                ioException.printStackTrace();
//            }
//        });
        writeVariantsNum(variantsNumSummaryFile, chrSnpNumMap);
        System.out.println("extract variants with ancestral state complicated in "+ Benchmark.getTimeSpanHours(start0)+ " hours");
    }

//    public static void go(String chrSNPNumFile, String vcfInputDir, int variantsNum, String pgfFileGeneHC,
//                          String ancestraDir, String outDir){
//        String[] subDirArray={"001_vmap2_"+variantsNum/1000+"k","002_vmap2_gene","003_vmap2_ancestral"};
//        String[] snpNumSummaryFileArray = {subDirArray[0]+".summary.txt", subDirArray[1]+".summary.txt",
//                subDirArray[2]+".summary.txt"};
//        for (String str: subDirArray){
//            new File(outDir, str).mkdirs();
//        }
//        VcfSubsetUtils.sampleFixedNumVariants(chrSNPNumFile, vcfInputDir, variantsNum, new File(outDir, subDirArray[0]).getAbsolutePath(),
//                new File(outDir, snpNumSummaryFileArray[0]).getAbsolutePath());
//        VcfSubsetUtils.extractVariantsInGene(pgfFileGeneHC, vcfInputDir, new File(outDir, subDirArray[1]).getAbsolutePath(),
//                new File(outDir, snpNumSummaryFileArray[1]).getAbsolutePath());
//        VcfSubsetUtils.extractVariantsWithAncestralState(ancestraDir, vcfInputDir, new File(outDir, subDirArray[2]).getAbsolutePath(),
//                new File(outDir, snpNumSummaryFileArray[2]).getAbsolutePath());
//    }
}
