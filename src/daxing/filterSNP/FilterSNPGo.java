package daxing.filterSNP;

import daxing.common.CollectionTool;
import pgl.infra.pos.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class FilterSNPGo {

    /**
     * 根据从VCF得到的DepthDB, 以top density 0.75为标准过滤VCF，输出chr pos
     * @param depthsDBDir chr pos depth sd
     * @param outPutDir chr pos
     * @param topDensity 0.75
     */
    public static void getTopDensity(String depthsDBDir, String outPutDir, double topDensity){
        File[] input= IOUtils.listRecursiveFiles(new File(depthsDBDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).sorted().toArray(File[]::new);
        String[] outFileName= Arrays.stream(files).map(File::getName).map(str->str.replaceAll("depth.txt.gz$", "chrpos.txt.gz")).sorted().toArray(String[]::new);
        IntStream.range(0, files.length).parallel().forEach(e->{
            Cells.getTopCellDensity(files[e], new File(outPutDir, outFileName[e]), topDensity);
        });
    }

    public static void getFilteredVCF(String vcfDir, String chrPosDir, String filteredVCFOutDir){
        File[] vcfFilesInput=IOUtils.listRecursiveFiles(new File(vcfDir));
        File[] chrPosFilesInput=IOUtils.listRecursiveFiles(new File(chrPosDir));
        Predicate<File> p=File::isHidden;
        File[] vcfFiles=Arrays.stream(vcfFilesInput).filter(p.negate()).sorted().toArray(File[]::new);
        File[] chrPosFiles=Arrays.stream(chrPosFilesInput).filter(p.negate()).sorted().toArray(File[]::new);
        if (vcfFiles.length!=chrPosFiles.length){
            System.out.println("please check vcfDir and chrPosDir, they must have equal number files");
            System.exit(1);
        }
        String[] filteredVCFOutName=Arrays.stream(chrPosFiles).map(File::getName)
                .map(str->str.replaceAll("chrpos.txt.gz$", "filtered0.75.vcf.gz")).sorted().toArray(String[]::new);
        IntStream.range(0, vcfFiles.length).parallel().forEach(e->{
            FilterSNPGo.getFilteredVCFfile(vcfFiles[e], chrPosFiles[e], new File(filteredVCFOutDir, filteredVCFOutName[e]));
        });
    }

    private static void getFilteredVCFfile(File vcfFile, File chrPosFile, File filteredVCFOutFile){
        String line=null;
        try(BufferedReader br1=IOUtils.getTextReader(vcfFile.getAbsolutePath());
            BufferedReader br2=IOUtils.getTextGzipReader(chrPosFile.getAbsolutePath());
            BufferedWriter bw=IOUtils.getTextGzipWriter(filteredVCFOutFile.getAbsolutePath())){
            long start=System.nanoTime();
            ChrPos[] chrPosArray=br2.lines().skip(1).map(PStringUtils::fastSplit)
                    .map(e->new ChrPos(Short.parseShort(e.get(0)), Integer.parseInt(e.get(1)))).sorted().toArray(ChrPos[]::new);
            List<String> lineList;
            short chr;
            int pos;
            ChrPos chrPos;
            int index;
            StringBuilder sb=new StringBuilder();
            while ((line=br1.readLine()).startsWith("##")){
                sb.append(line).append("\n");
            }
            sb.append(line).append("\n");
            bw.write(sb.toString());
            while ((line=br1.readLine())!=null){
                lineList=PStringUtils.fastSplit(line.substring(0, 20));
                chr= Short.parseShort(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                index=Arrays.binarySearch(chrPosArray, new ChrPos(chr, pos));
                if(index<0) continue;
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
            System.out.println(filteredVCFOutFile.getName()+" is complicated in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
            System.out.println(vcfFile.getName()+"\t"+line);
        }
    }

    /**
     *  最多只保留2个Alt allele
     * @param abdVCFDir1
     * @param abVCFDir
     * @param dVCFDir
     * @param outFile Chr Pos	Ref	Alt 有header
     */
    public static void mergePosList (String abdVCFDir1, String abVCFDir, String dVCFDir, String outFile) {
        File[] inputABD = IOUtils.listRecursiveFiles(new File(abdVCFDir1));
        File[] inputAB = IOUtils.listRecursiveFiles(new File(abVCFDir));
        File[] inputD = IOUtils.listRecursiveFiles(new File(dVCFDir));
        Predicate<File> p = File::isHidden;
        File[] abdFiles = Arrays.stream(inputABD).filter(p.negate()).toArray(File[]::new);
        File[] abFiles = Arrays.stream(inputAB).filter(p.negate()).toArray(File[]::new);
        File[] dFiles = Arrays.stream(inputD).filter(p.negate()).toArray(File[]::new);
        List<File> files = new ArrayList<>();
        String[] filesName = Arrays.stream(abdFiles).map(File::getName).map(str -> str.replaceAll(".ABDgenome.filtered0.75.vcf$", "_PosAllele.txt.gz")).toArray(String[]::new);
        files.addAll(CollectionTool.changeToList(abdFiles));
        files.addAll(CollectionTool.changeToList(abFiles));
        files.addAll(CollectionTool.changeToList(dFiles));
        Comparator<File> comparator=Comparator.comparing(f->f.getName());
        Collections.sort(files, comparator);
        int[] aa = IntStream.iterate(0, n -> n + 2).limit(42).toArray();
        Arrays.stream(aa).forEach(e -> FilterSNPGo.mergePosList(files.get(e), files.get(e + 1), new File(outFile, filesName[e / 2])));
    }

    /**
     * 最多只保留2个Alt allele
     * @param vcfInputFile1
     * @param vcfInputFile2
     * @param outFile chr pos  ref alt
     */
    public static void mergePosList (File vcfInputFile1, File vcfInputFile2, File outFile) {
        long start=System.nanoTime();
        String inFileS1 = vcfInputFile1.getAbsolutePath();
        String inFileS2 = vcfInputFile2.getAbsolutePath();
        String outfileS = outFile.getAbsolutePath();
        String[] alleles = {"A", "C", "G", "T", "D", "I"};
        Arrays.sort(alleles);
        double[] fre = new double[6];
        try {
            int chr1 = Integer.MIN_VALUE;
            int chr2 = Integer.MIN_VALUE;
            int taxaNum1 = Integer.MIN_VALUE;
            int taxaNum2 = Integer.MIN_VALUE;
            TIntArrayList posList1 = new TIntArrayList();
            List<String> referList1 = new ArrayList<>();
            List<String> altList1 = new ArrayList<>();
            List<String> altDepthList1 = new ArrayList<>();
            TIntArrayList posList2 = new TIntArrayList();
            List<String> referList2 = new ArrayList<>();
            List<String> altList2 = new ArrayList<>();
            List<String> altDepthList2 = new ArrayList<>();
            BufferedReader br = IOUtils.getTextReader(inFileS1);
            String temp = null;
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum1 = temp.split("\t").length-9;
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 100);
                tem = temp.split("\t");
                chr1 = Integer.parseInt(tem[0]);
                posList1.add(Integer.parseInt(tem[1]));
                referList1.add(tem[3]);
                altList1.add(tem[4]);
                altDepthList1.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            br.close();
            br = IOUtils.getTextReader(inFileS2);
            temp = null;
            while ((temp = br.readLine()).startsWith("##")) {};
            taxaNum2 = temp.split("\t").length-9;
            double weight1 = (double)taxaNum1/(taxaNum1+taxaNum2);
            double weight2 = (double)taxaNum2/(taxaNum1+taxaNum2);
            tem = null;
            while ((temp = br.readLine()) != null) {
                temp = temp.substring(0, 100);
                tem = temp.split("\t");
                chr2 = Integer.parseInt(tem[0]);
                posList2.add(Integer.parseInt(tem[1]));
                referList2.add(tem[3]);
                altList2.add(tem[4]);
                altDepthList2.add(tem[7].split(";")[1].replace("AD=", ""));
            }
            br.close();
            if (chr1 != chr2) {
                System.out.println("Wrong input files! Program quits.");
                System.exit(0);
            }
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            bw.write("Chr\tPos\tRef\tAlt");
            bw.newLine();
            TIntHashSet mergedPosSet = new TIntHashSet(posList1);
            mergedPosSet.addAll(posList2);
            int[] mergedPos = mergedPosSet.toArray();
            Arrays.sort(mergedPos);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < mergedPos.length; i++) {
                sb = new StringBuilder();
                int index1 = posList1.binarySearch(mergedPos[i]);
                int index2 = posList2.binarySearch(mergedPos[i]);
                if (index1 < 0 && index2 > -1) {
                    if (altList2.get(index2).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList2.get(index2)).append("\t").append(referList2.get(index2)).append("\t").append(altList2.get(index2));
                }
                else if (index1 > -1 && index2 < 0) {
                    if (altList1.get(index1).length() > 3) continue;
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t").append(altList1.get(index1));
                }
                else  {
                    for (int j = 0; j < fre.length; j++) {
                        fre[j] = -1;
                    }
                    tem = altList1.get(index1).split(",");
                    String[] fretem = altDepthList1.get(index1).split(",");
                    double[] depth = new double[fretem.length];
                    double[] fre1 = new double[depth.length];
                    double sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre1[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        fre[index] = fre1[j+1]*weight1;
                    }

                    tem = altList2.get(index2).split(",");
                    fretem = altDepthList2.get(index2).split(",");
                    depth = new double[fretem.length];
                    double[] fre2 = new double[depth.length];
                    sum = 0;
                    for (int j = 0; j < depth.length; j++) {
                        depth[j] = Integer.parseInt(fretem[j]);
                        sum+=depth[j];
                    }
                    for (int j = 0; j < depth.length; j++) {
                        fre2[j] = depth[j]/sum;
                    }
                    for (int j = 0; j < tem.length; j++) {
                        int index = Arrays.binarySearch(alleles, tem[j]);
                        if (fre[index] < 0) {
                            fre[index] = fre2[j+1]*weight2;
                        }
                        else {
                            fre[index] = fre[index]+fre2[j+1]*weight2;
                        }
                    }

                    int[] indices = PArrayUtils.getIndicesByDescendingValue(fre);
                    sb.append(chr1).append("\t").append(posList1.get(index1)).append("\t").append(referList1.get(index1)).append("\t");
                    for (int j = 0; j < 2; j++) {
                        if (fre[indices[j]] > 0) {
                            sb.append(alleles[indices[j]]).append(",");
                        }
                        else {
                            break;
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            System.out.println(outFile.getName()+" is completed in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

//    public static void main(String[] args) {
//        FilterSNPGo.getTopDensity(args[0], args[1], Double.parseDouble(args[2]));
//        FilterSNPGo.getFilteredVCF(args[0], args[1], args[2]);
//        FilterSNPGo.mergePosList(args[0], args[1], args[2], args[3]);
//    }
}
