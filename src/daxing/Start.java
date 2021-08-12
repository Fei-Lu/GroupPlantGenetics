package daxing;

import daxing.common.ChrRange;
import daxing.common.IOTool;
import daxing.load.ancestralSite.SNPAnnotation;
import daxing.load.complementary.ComplementaryGo;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {



////        IOTool.viewHeader("/Users/xudaxing/Desktop/test/003_parameterFile/017_WheatVMap2_GermplasmInfo.txt");
//
        String exonSNPAnnoDir=args[0];
        String exonVCFDir=args[1];
        String taxa_InfoDBFile=args[2];
        String triadFile=args[3];
        String nonoverlapGeneFile=args[4];
        int blockGeneNum=Integer.parseInt(args[5]);
        String pgfFile=args[6];
        String outDir=args[7];
//////
//        String exonSNPAnnoDir="/Users/xudaxing/Desktop/test/001_geneSiteaAnnoDB";
//        String exonVCFDir="/Users/xudaxing/Desktop/test/002_geneVCF";
//        String taxa_InfoDBFile="/Users/xudaxing/Desktop/test/003_parameterFile/017_WheatVMap2_GermplasmInfo.txt";
//        String triadFile="/Users/xudaxing/Desktop/test/003_parameterFile/triadGenes1.1.txt";
//        String nonoverlapGeneFile="/Users/xudaxing/Desktop/test/003_parameterFile/wheat_v1.1_nonoverlap.txt";
//        int blockGeneNum=Integer.parseInt("20");
//        String pgfFile="/Users/xudaxing/Desktop/test/003_parameterFile/wheat_v1.1_Lulab.pgf";
//        String outDir="/Users/xudaxing/Desktop/test/004_outDir";
        SNPAnnotation.MethodCallDeleterious methodCallDeleterious= SNPAnnotation.MethodCallDeleterious.SIFT_GERP;
        ComplementaryGo.go(exonSNPAnnoDir,exonVCFDir, taxa_InfoDBFile, triadFile, nonoverlapGeneFile, blockGeneNum, pgfFile, outDir,
                methodCallDeleterious);

//        String inputFile_bin=args[0];
//        String positionsDir_42=args[1];
//        String maxIndexPerChrFile=args[2];
//        String alleleAgeSH=args[3];
//        AlleleAge.buildSH_alleleAge(inputFile_bin, positionsDir_42, maxIndexPerChrFile, alleleAgeSH);

//        ScriptMethods.splitSh(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));

//        ScriptMethods.splitSh(args[0], args[1], Integer.parseInt(args[2]), Integer.parseInt(args[3]));
    }



    public static void recombination(String inputFile1, String inputFile2, String outFile){
        try (BufferedReader br1 = IOTool.getReader(inputFile1);
             BufferedReader br2 = IOTool.getReader(inputFile2);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)");
            bw.newLine();
            String line1, line2;
            List<String> temp1, temp2;
            br1.readLine();
            br2.readLine();
            String refChr;
            int refPos, start, end;
            int chrID, pos;
            double geneticsPos, recombinationRate;
            List<ChrRange> chrRangeList = new ArrayList<>();
            DoubleList recombinationRateList = new DoubleArrayList();
            ChrRange chrRange;
            StringBuilder sb = new StringBuilder();
            while ((line2= br2.readLine())!=null){
                temp2 = PStringUtils.fastSplit(line2);
                refChr = temp2.get(0).substring(3,5);
                start = Integer.parseInt(temp2.get(1));
                end = Integer.parseInt(temp2.get(2));
                recombinationRate = Double.parseDouble(temp2.get(4));
                chrRange = new ChrRange(refChr, start, end+1);
                chrRangeList.add(chrRange);
                recombinationRateList.add(recombinationRate);
            }
            while((line1= br1.readLine())!=null){
                temp1 = PStringUtils.fastSplit(line1);
                refChr = temp1.get(1).substring(3,5);
                refPos = Integer.parseInt(temp1.get(2));
                geneticsPos= Double.parseDouble(temp1.get(3));
                chrRange = new ChrRange(refChr, refPos, refPos+1);
                int hit = Collections.binarySearch(chrRangeList, chrRange);
                int index =  hit < 0 ? -hit-2 : hit;
                recombinationRate = recombinationRateList.getDouble(index);
                chrID = RefV1Utils.getChrID(refChr, refPos);
                pos = RefV1Utils.getPosOnChrID(refChr, refPos);
                sb.setLength(0);
                sb.append(chrID).append("\t").append(pos).append("\t").append(recombinationRate).append("\t");
                sb.append(geneticsPos);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void biAllele(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, "vcf.gz");
        String[] outNames=
                files.stream().map(File::getName).map(s -> s.replaceAll("2.0.vcf.gz","2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, subLine;
                List<String> temp;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write("#CHROM\t"+line.substring(5));
                bw.newLine();
                int cnt=0;
                Set<String> altSet=new HashSet<>();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0,40);
                    temp= PStringUtils.fastSplit(subLine);
                    if (temp.get(4).length() > 1){
                        altSet.add(temp.get(4));
                    }
                    if (temp.get(4).equals("A") || temp.get(4).equals("T") ||
                            temp.get(4).equals("C") || temp.get(4).equals("G") ||
                            temp.get(4).equals("I") || temp.get(4).equals("D")){
                        bw.write(line);
                        bw.newLine();
                        cnt++;
                    }else {
                        altSet.add(temp.get(4));
                    }
                }
                bw.flush();
//                altSet.stream().forEach(System.out::println);
                System.out.println(new File(outDir, outNames[e]).getName()+  " have "+cnt+" biVariants");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void correctMaf(String vmap2Dir, String outDir){
        List<File> files =IOTool.getFileListInDirEndsWith(vmap2Dir, ".vcf.gz");
        String[] outNames=
                files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp,tem,te;
                double maf;
                StringBuilder sb = new StringBuilder();
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    tem=PStringUtils.fastSplit(temp.get(7), ";");
                    te=PStringUtils.fastSplit(tem.get(6), "=");
                    maf = Double.parseDouble(te.get(1));
                    maf = maf > 0.5 ? 1-maf : maf;
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0, 7))).append("\t");
                    sb.append(String.join(";", tem.subList(0, 6))).append(";");
                    sb.append(te.get(0)).append("=").append(maf).append("\t");
                    sb.append(String.join("\t", temp.subList(8, temp.size())));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }


    public static void extract(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            String line, subLine, ref, alt;
            List<String> temp, tem, te, t;
            String[] acgt={"A","C","G","T"};
            int[] countACGT;
            double maf;
            String chr;
            int refIndex, altIndex, minorCount, majorCount, chrID, pos, refPos, majorIndex, minorIndex;
            int refAlleleCount, altAlleleCount;
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                while ((line = br.readLine()).startsWith("##")) continue;
                StringBuilder sb= new StringBuilder();
                while ((line=br.readLine())!=null){
                    subLine= line.substring(0,150);
                    temp=PStringUtils.fastSplit(subLine);
                    ref = temp.get(3);
                    alt = temp.get(4);
                    tem= PStringUtils.fastSplit(temp.get(7), ";");
                    te=PStringUtils.fastSplit(tem.get(6), "=");
                    maf = Double.parseDouble(te.get(1));
                    te = PStringUtils.fastSplit(tem.get(4), "=");
                    t = PStringUtils.fastSplit(te.get(1),  ",");
                    refAlleleCount = Integer.parseInt(t.get(0)) * 2 + Integer.parseInt(t.get(1));
                    altAlleleCount = Integer.parseInt(t.get(2)) * 2 + Integer.parseInt(t.get(1));
                    refIndex = Arrays.binarySearch(acgt, ref);
                    altIndex = Arrays.binarySearch(acgt, alt);
                    minorCount = (int)Math.round(199*maf);
                    majorCount = 199-minorCount;
                    minorIndex = refAlleleCount > altAlleleCount ? altIndex : refIndex;
                    majorIndex = refAlleleCount > altAlleleCount ? refIndex : altIndex;
                    countACGT = new int[acgt.length];
                    countACGT[minorIndex]=minorCount;
                    countACGT[majorIndex]=majorCount;
                    sb.setLength(0);
                    chrID = Integer.parseInt(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    chr= RefV1Utils.getChromosome(chrID, pos);
                    refPos = RefV1Utils.getPosOnChromosome(chrID, pos);
                    sb.append("chr").append(chr).append("\t").append(refPos-1).append("\t").append(refPos).append("\t");
                    for (int i = 0; i < countACGT.length; i++) {
                        sb.append(countACGT[i]).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }


}