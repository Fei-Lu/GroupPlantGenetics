package daxing.hybrid.dection;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Utils {



    public static void getMissingRateByTaxon(String inputDir, String taxaNameFile, String outDir){
        List<File> files= IOTool.getVisibleFileListInDir(inputDir);
        String[] outFileNames= files.stream().map(file ->
                file.getName().replaceAll("vcf","missingRate.byTaxon")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->getMissingRateByTaxon(files.get(e), taxaNameFile, new File(outDir,
                outFileNames[e])));
    }

    private static void getMissingRateByTaxon(File inputVCFFile, String taxaNameFile,File outFile){
        GenotypeGrid genotypeGrid=new GenotypeGrid(inputVCFFile.getAbsolutePath(), GenoIOFormat.VCF);
        List<String> taxa= RowTableTool.getColumnList(taxaNameFile, 0);
        int[] index=new int[taxa.size()];
        for (int i = 0; i < taxa.size(); i++) {
            index[i]=genotypeGrid.getTaxonIndex(taxa.get(i));
        }
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxa\tMissingRateByTaxon");
            bw.newLine();
            StringBuilder stringBuilder=new StringBuilder();
            for (int i = 0; i < index.length; i++) {
                double missingNum=genotypeGrid.getMissingNumberByTaxon(index[i]);
                double total=genotypeGrid.getSiteNumber();
                stringBuilder.setLength(0);
                stringBuilder.append(taxa.get(i)).append("\t").append(missingNum/total);
                bw.write(stringBuilder.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getMissingRateBySite(String inputVCFDir, String outDir){
        List<File> files= IOTool.getVisibleFileListInDir(inputVCFDir);
        String[] outFileNames= files.stream().map(file ->
                file.getName().replaceAll("vcf","missingRate.bySites")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->getMissingRateBySite(files.get(e), new File(outDir,
                outFileNames[e])));
    }

    private static void getMissingRateBySite(File inputVCFFile, File outFile){
        long start=System.nanoTime();
        BufferedWriter bw = IOTool.getWriter(outFile);
        try {
            bw.write("Chr\tPos\tMissingRate");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            double missingNum, taxaNum;
            NumberFormat numberFormat=NumberFormat.getInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            GenotypeGrid genotypeGrid=new GenotypeGrid(inputVCFFile.getAbsolutePath(),  GenoIOFormat.VCF);
            for (int i = 0; i < genotypeGrid.getSiteNumber(); i++) {
                missingNum=genotypeGrid.getMissingNumberBySite(i);
                taxaNum=genotypeGrid.getTaxaNumber();
                sb.setLength(0);
                sb.append(genotypeGrid.getChromosome(i)).append("\t");
                sb.append(genotypeGrid.getPosition(i)).append("\t");
                sb.append(numberFormat.format(missingNum/taxaNum));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(inputVCFFile.getName() + " completed in "+ Benchmark.getTimeSpanSeconds(start) + " " +
                "s");
    }

    /**
     * taxon1 (sorted by taxon name) vs all other
     * @param inputDir vcfDir
     * @param taxonNameFile Specified taxon list
     * @param outDir
     */
    public static void getSeparationSiteBetweenParents(String inputDir, String taxonNameFile,String outDir,
                                                       int depthThreshold){
        List<File> fileList= IOTool.getFileListInDir(inputDir);
        Predicate<File> p0= file -> file.getName().startsWith("00", 4);
        Predicate<File> p43=file -> file.getName().startsWith("43", 4);
        Predicate<File> p44=file -> file.getName().startsWith("44", 4);
        Predicate<File> p=p0.or(p43).or(p44);
        List<File> files=fileList.stream().filter(p.negate()).collect(Collectors.toList());
        List<String> taxaList=getTaxaList(taxonNameFile);
        int[][] snpCount= new int[files.size()][];
        for (int i = 0; i < files.size(); i++) {
            snpCount[i]=countSNP(files.get(i).getAbsolutePath(),taxaList,depthThreshold);
        }
        try (BufferedWriter bw = IOTool.getWriter(new File(outDir, "snpCount.GBS.txt"))) {
            int[] chrIDs= RefV1Utils.getChrIDs();
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < taxaList.size(); i++) {
                sb.append(taxaList.get(i)).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write("chr\tTaxa\t"+sb.toString());
            bw.newLine();
            for (int i = 0; i < snpCount.length; i++) {
                sb.setLength(0);
                sb.append(chrIDs[i]).append("\t").append(taxaList.get(0)).append("\t");
                for (int j = 0; j < snpCount[i].length; j++) {
                    sb.append(snpCount[i][j]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static List<String> getTaxaList(String taxonNameFile){
        List<String> res=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(taxonNameFile)) {
            String line;
            List<String> temp;
//            br.readLine();
//            StringBuilder sb=new StringBuilder();
            int count=0;
            while ((line=br.readLine())!=null){
                count++;
                temp= PStringUtils.fastSplit(line);
//                sb.setLength(0);
//                sb.append(temp.get(0)).append("_").append(temp.get(1)).append("_");
//                sb.append(temp.get(2)).append("_").append(temp.get(3)).append("_");
//                sb.append(temp.get(4));
                if (count > 30) break;
                res.add(temp.get(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        Collections.sort(res);
        return res;
    }

    private static List<String> getTaxaList_vmap2(String taxonNameFile){
        List<String> res=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(taxonNameFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                res.add(temp.get(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        Collections.sort(res);
        return res;
    }

    /**
     *
     * @param vcfFile
     * @param taxonList
     * @return
     */
    private static int[] countSNP(String vcfFile, List<String> taxonList, int deptThreshold){
        long start=System.nanoTime();
        int[] res=new int[taxonList.size()];
        try (BufferedReader br = IOTool.getReader(vcfFile)) {
            String line;
            List<String> temp, tem, te;
            int depth;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            int index= temp.indexOf(taxonList.get(0));
            int[] taxaIndex=new int[taxonList.size()];
            for (int i = 0; i < taxonList.size(); i++) {
                taxaIndex[i]=temp.indexOf(taxonList.get(i));
            }
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                if (temp.get(index).startsWith("./.")) continue;
                if (temp.get(index).startsWith("0/1")) continue;
                for (int i = 0; i < taxaIndex.length; i++) {
                    if (temp.get(taxaIndex[i]).startsWith("./.")) continue;
                    if (temp.get(taxaIndex[i]).startsWith("0/1")) continue;
                    if (temp.get(index).startsWith("0/0") && temp.get(taxaIndex[i]).startsWith("0/0")) continue;
                    if (temp.get(index).startsWith("1/1") && temp.get(taxaIndex[i]).startsWith("1/1")) continue;
//                    res[i]++;
                    tem=PStringUtils.fastSplit(temp.get(taxaIndex[i]), ":");
                    te=PStringUtils.fastSplit(tem.get(0), ",");
                    depth=Integer.parseInt(te.get(0))+ Integer.parseInt(te.get(1));
                    if (depth < deptThreshold) continue;
                    res[i]++;
                }
            }
            System.out.println(vcfFile + " has completed in "+ Benchmark.getTimeSpanSeconds(start)+ " s");
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }
}
