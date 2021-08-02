package daxing;

import daxing.common.IOTool;
import daxing.common.VCF;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {
//        String vcfInputDir=args[0];
//        String pgfFile=args[1];
//        String nonOverlapGeneFile=args[2];
//        String ancestralDir=args[3];
//        String outDir=args[4];
//        GeneSiteAnnoDB.extractGeneSiteInfoForAnnoDB(vcfInputDir, pgfFile, nonOverlapGeneFile, ancestralDir, outDir);
        ChrSNPAnnoDB.splitSNPAnnoDB("/Users/xudaxing/Documents/deleteriousMutation/001_analysis/005_vmap2_1000/003_Gene_Variant_genotype/002_variantsAnnotation/004_GeneSiteAnno/geneSiteAnno.txt.gz",
                "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/005_vmap2_1000/003_Gene_Variant_genotype/002_variantsAnnotation/004_GeneSiteAnno/004_byChrID_parsimony_leftJoinAoyue_vepSnpeff");
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

    public static void filterMaf(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");
        String[] outNames= files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".maf0.1.vcf.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e -> {
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    if (Double.parseDouble(VCF.calculateMaf(line)) < 0.1) continue;
                    bw.write(line);
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