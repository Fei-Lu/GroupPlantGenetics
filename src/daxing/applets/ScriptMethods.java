package daxing.applets;

import daxing.common.*;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class ScriptMethods {

    public static void getTopRows(File inputFile, int n, File outputFile){
        try{
            long start=System.nanoTime();
            BufferedReader br;
            BufferedWriter bw;
            if (inputFile.getName().endsWith("gz")){
                br=IOUtils.getTextGzipReader(inputFile.getAbsolutePath());
            }else {
                br=IOUtils.getTextReader(inputFile.getAbsolutePath());
            }
            if (outputFile.getName().endsWith(".gz")){
                bw=IOUtils.getTextGzipWriter(outputFile.getAbsolutePath());
            }
            else {
                bw=IOUtils.getTextWriter(outputFile.getAbsolutePath());
            }
            String line;
            int count=1;
            while ((line=br.readLine())!=null){
                if (count > n) break;
                bw.write(line);
                bw.newLine();
                count++;
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(outputFile.getName()+" is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void getTopRowsFromDir(String inputDir, int n, String outputDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] names= files.stream().map(File::getName).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e-> ScriptMethods.getTopRows(files.get(e), n, new File(outputDir, names[e])));
    }

    /**
     *
     * @param inputFile 有header
     * @param rate
     * @param subsetFile 有header
     */
    public static void getSubsetFromFile(String inputFile, double rate, String subsetFile){
        long start=System.nanoTime();
        BufferedReader br=null;
        BufferedWriter bw = null;
        try {
            if (inputFile.endsWith("gz")){
                br=IOUtils.getTextGzipReader(inputFile);
                bw=IOUtils.getTextGzipWriter(subsetFile);
            }else {
                br=IOUtils.getTextReader(inputFile);
                bw=IOUtils.getTextWriter(subsetFile);
            }
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            double r=-1d;
            int count=0;
            int total=0;
            while ((line=br.readLine())!=null){
                total++;
                r=Math.random();
                if (r>rate) continue;
                count++;
                bw.write(line);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("samping "+NumberTool.parse(count)+"("+NumberTool.parse(total)+") row from "
                    +new File(inputFile).getName()+" into "+new File(subsetFile).getName()+" in "
                    +Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void getSubsetFromDir(String inputDir, double rate, String outDir){
        long start=System.nanoTime();
        File[] input=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        String[] files= Arrays.stream(input).filter(p.negate()).map(File::getAbsolutePath).toArray(String[]::new);
        String[] filesName=Arrays.stream(input).filter(p.negate()).map(File::getName)
                .map(str->str.replaceAll("vcf", "subset.vcf")).toArray(String[]::new);
        IntStream.range(0, files.length).forEach(e->ScriptMethods.getSubsetFromFile(files[e], rate, new File(outDir,
                filesName[e]).getAbsolutePath()));
        System.out.println(outDir+" subset is completed in "+Benchmark.getTimeSpanHours(start)+" hours");
    }

    public static void getSubsetFromDir(String inputDir, String outDir){
        ScriptMethods.getSubsetFromDir(inputDir, 0.001, outDir);
    }

    public static void calculateLD(String ingputDir, int binWidth_kb, int threshForDistance_Mb, String outDir){
        File[] files= IOUtils.listRecursiveFiles(new File(ingputDir));
        Predicate<File> h= File::isHidden;
        File[] f=Arrays.stream(files).filter(h.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f).map(File::getName).toArray(String[]::new);
        IntStream.range(0, f.length).parallel().forEach(e->{
            calculateLD(f[e], binWidth_kb, threshForDistance_Mb, new File(outDir, outNames[e]));
        });
    }

    public static void calculateLD(File ingputFile, int binWidth_kb, int threshForDistance_Mb, File outFile){
        try (BufferedReader br = IOUtils.getTextReader(ingputFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> temp;
            int distance=Integer.MIN_VALUE;
            double r2=Double.MIN_VALUE;
            int thresh=binWidth_kb*1000;
            int limit=threshForDistance_Mb*1000000;
            int kb=binWidth_kb;
            DescriptiveStatistics r2Stats=new DescriptiveStatistics();
            br.readLine();
            bw.write("window_kb\tnumberInWindow\tmeanOfR2\n");
            StringBuilder sb;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                distance=Integer.parseInt(temp.get(0));
                if (distance > limit) break;
                r2=Double.parseDouble(temp.get(1));
                if (distance > thresh){
                    sb=new StringBuilder();
                    sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
                    bw.write(sb.toString());
                    bw.newLine();
                    thresh+=binWidth_kb*1000;
                    kb+=binWidth_kb;
                    r2Stats=new DescriptiveStatistics();
                }
                r2Stats.addValue(r2);
            }
            sb=new StringBuilder();
            sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void getAncestral(String inputDir, String ancestralDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName).map(str->str.replaceAll("\\.txt",".ancestral.txt"))
                .toArray(String[]::new);
        String[] subgenome=files.stream().map(File::getName).map(str->str.substring(0,1)).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->getAncestral(files.get(e),new File(ancestralDir, outNames[e]),
                subgenome[e]));
    }

    public static void getAncestral(File inputFile, File ancestralFile, String subgenome){
        String[] acgt={"A","C","G","T"};
        try (BufferedReader bufferedReader = IOTool.getReader(inputFile);
             BufferedWriter bufferedWriter=IOTool.getTextWriter(ancestralFile)) {
            String line;
            List<String> temp, secer, hv;
            int chrID, pos, indexOfSecer, indexOfHv;
            bufferedWriter.write("chr\tpos\tancestral\n");
            StringBuilder sb;
            while ((line=bufferedReader.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chrID= RefV1Utils.getChrID(temp.get(0).substring(3,4)+subgenome, Integer.parseInt(temp.get(2)));
                pos=RefV1Utils.getPosOnChrID(temp.get(0).substring(3,4)+subgenome, Integer.parseInt(temp.get(2)));
                secer=PStringUtils.fastSplit(temp.get(4),",");
                hv=PStringUtils.fastSplit(temp.get(5), ",");
                indexOfSecer=secer.indexOf("1");
                indexOfHv=hv.indexOf("1");
                if (indexOfHv!=indexOfSecer) continue;
                sb=new StringBuilder();
                sb.append(chrID).append("\t").append(pos).append("\t").append(acgt[indexOfHv]);
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void retainVmapIExonGeno(String pgfFile, String geneHCFile, String genoInputDir, String outDir){
        PGF pgf=new PGF(pgfFile);
        System.out.println(pgf.getGeneNumber());
        Set<String> hcGenes=RowTableTool.getColumnSet(geneHCFile,0);
        Predicate<PGF.Gene> hcGenePredict=gene -> hcGenes.contains(gene.getGeneName());
        pgf.removeIf(hcGenePredict.negate());
        System.out.println(pgf.getGeneNumber());
        List<File> genoFiles=IOUtils.getVisibleFileListInDir(genoInputDir);
        try {
            BufferedReader bufferedReader;
            BufferedWriter bufferedWriter;
            String line, outFileName, header;
            List<String> temp;
            int chrID, vcfPos, geneIndex;
            for (int i = 0; i < genoFiles.size(); i++) {
                outFileName=genoFiles.get(i).getName().replaceAll("geno","exon.geon");
                bufferedReader=IOTool.getReader(genoFiles.get(i));
                bufferedWriter=IOTool.getTextGzipWriter(new File(outDir, outFileName));
                header=bufferedReader.readLine();
                bufferedWriter.write(header);
                bufferedWriter.newLine();
                while ((line=bufferedReader.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID= RefV1Utils.getChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    vcfPos=RefV1Utils.getPosOnChrID(temp.get(0), Integer.parseInt(temp.get(1)));
                    geneIndex=pgf.getGeneIndex(chrID, vcfPos);
                    if (geneIndex<0) continue;
                    if (!pgf.isWithinThisGeneExon(geneIndex, chrID, vcfPos)) continue;
                    bufferedWriter.write(line);
                    bufferedWriter.newLine();
                }
                bufferedReader.close();
                bufferedWriter.flush();
                bufferedWriter.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
