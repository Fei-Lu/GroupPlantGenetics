package daxing.common;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import org.apache.commons.math3.util.CombinatoricsUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.text.NumberFormat;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Daxing Xu
 */
public class VCF {

    private String meta;
    private List<String> header;
    private List<List<String>> data;

    public VCF(String inputFile){
        this.initilize(inputFile);
    }

    public VCF(File inputFile){
        this.initilize(inputFile.getAbsolutePath());
    }

    public VCF(String meta, List<String> header, List<List<String>> data){
        this.meta=meta;
        this.header=header;
        this.data=data;
    }

    private void initilize(String inputFile){
        BufferedReader br;
        try{
            if (inputFile.endsWith("gz")){
                br= IOUtils.getTextGzipReader(inputFile);
            }else {
                br=IOUtils.getTextReader(inputFile);
            }
            StringBuilder sb=new StringBuilder();
            String temp;
            List<List<String>> lists=new ArrayList<>(1000);
            while ((temp=br.readLine()).startsWith("##")){
                sb.append(temp).append("\n");
            }
            this.meta =sb.toString();
            this.header=PStringUtils.fastSplit(temp);
            while ((temp=br.readLine())!=null){
                lists.add(PStringUtils.fastSplit(temp));
            }
            this.data=lists;
            br.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 根据VCF目录和染色体编号，返回对应染色体编号的各个VCF文件输入路径
     * @param vcfInputDir VCF目录
     * @param chrArray 染色体编号
     * @return 染色体编号对应的所有VCF文件输入路径
     */
    public static List<String> getAllVcfInputPath(String vcfInputDir, int[] chrArray){
        int temp=(int)Arrays.stream(chrArray).distinct().count();
        if(chrArray.length>temp){
            System.out.println("please check your input array, it contains duplicate value");
            System.exit(1);
        }
        File[] files=IOUtils.listRecursiveFiles(new File(vcfInputDir));
        List<Integer> chrList=Arrays.stream(chrArray).boxed().collect(Collectors.toList());
        return Arrays.stream(files).filter(e->chrList.contains(StringTool.getNumFromString(e.getName()))).map(File::getAbsolutePath).collect(Collectors.toList());
    }

    /**
     * merge chr001, chr002, chr003, ... to chr.Asubgenome.vcf, chr.Bsubgenome.vcf, chr.Dsubgenome.vcf
     * @param inputVcfDir inputVcfDir
     * @param outDir outDir
     */
    public static void mergeVCFtoLineage(String inputVcfDir, String outDir){
        File[] files=new File(inputVcfDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
        assert files != null;
        File[] f=Arrays.stream(files).filter(p).sorted().toArray(File[]::new);
        int[] chrIDArray=Arrays.stream(f).filter(p).map(File::getName).map(str->str.substring(3,6))
                .mapToInt(Integer::parseInt).toArray();
        int[][] lineage=new int[3][];
        lineage[0]=WheatLineage.valueOf("A").getChrID();
        lineage[1]=WheatLineage.valueOf("B").getChrID();
        lineage[2]=WheatLineage.valueOf("D").getChrID();
        TIntArrayList[] indexArray=new TIntArrayList[3];
        for (int i = 0; i < 3; i++) {
            indexArray[i]=new TIntArrayList();
        }
        int index;
        for (int i = 0; i < lineage.length; i++) {
            for (int j = 0; j < lineage[i].length; j++) {
                index=Arrays.binarySearch(chrIDArray, lineage[i][j]);
                if (index<0) continue;
                indexArray[i].add(index);
            }
        }
        VCF[] vcfArray=new VCF[3];
        for (int i = 0; i < indexArray.length; i++) {
            if (indexArray[i].size()==0) continue;
            vcfArray[i]=new VCF(f[indexArray[i].get(0)]);
            for (int j = 1; j < indexArray[i].size(); j++) {
                if (indexArray[i].size()==1) continue;
                vcfArray[i].addVCF(new VCF(f[indexArray[i].get(j)]));
            }
        }
        String[] outNames={"chr.Asubgenome.vcf", "chr.Bsubgenome.vcf", "chr.Dsubgenome.vcf"};
        for (int i = 0; i < vcfArray.length; i++) {
            if (vcfArray[i]==null) continue;
            vcfArray[i].changeToRefChr();
            vcfArray[i].write(outDir, outNames[i]);
        }
    }

    /**
     * merge chr001 chr002 chr003 chr004, ... to chr1A chr1B chr1D, ...
     * suitable for small VCF file
     * @param inputVcfDir inputVcfDir
     * @param outDir outDir
     */
    public static void mergeVCFtoChr(String inputVcfDir, String outDir){
        File[] files=new File(inputVcfDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
        assert files != null;
        File[] f=Arrays.stream(files).filter(p).sorted().toArray(File[]::new);
        int[] chrIDArray=Arrays.stream(f).filter(p).map(File::getName).map(str->str.substring(3,6))
                .mapToInt(Integer::parseInt).sorted().toArray();
        TIntArrayList indexList=new TIntArrayList();
        for (int i = 0; i < chrIDArray.length-1; i++) {
            if (NumberTool.isOdd(chrIDArray[i])){
                if (chrIDArray[i]+1!=chrIDArray[i+1]) continue;
                indexList.add(i);
                indexList.add(i+1);
                i=i+1;
            }
        }
        VCF vcf;
        String[] outNameArray= Arrays.stream(f).map(File::getName).map(str->str.substring(6)).toArray(String[]::new);
        String outName;
        for (int i = 0; i < indexList.size(); i=i+2) {
            vcf=new VCF(f[indexList.get(i)]);
            vcf.addVCF(new VCF(f[indexList.get(i+1)]));
            vcf.changeToRefChr();
            outName=vcf.getData().get(0).get(0);
            vcf.write(outDir, "chr"+outName+outNameArray[indexList.get(i)]);
        }
    }

    /**
     *
     * @param inputVcfDir inputVcfDir
     * @param outDir outDir
     */
    public static void fastMergeVCFtoLineage(String inputVcfDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow()+ " start");
        File[] files=new File(inputVcfDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
        assert files != null;
        File[] f=Arrays.stream(files).filter(p).sorted().toArray(File[]::new);
        TIntArrayList[] abd=new TIntArrayList[3];
        abd[0]=new TIntArrayList(WheatLineage.valueOf("A").getChrID());
        abd[1]=new TIntArrayList(WheatLineage.valueOf("B").getChrID());
        abd[2]=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> ap=fa->abd[0].contains(Integer.parseInt(fa.getName().substring(3,6)));
        Predicate<File> bp=fa->abd[1].contains(Integer.parseInt(fa.getName().substring(3,6)));
        Predicate<File> dp=fa->abd[2].contains(Integer.parseInt(fa.getName().substring(3,6)));
        File[][] abd_lineageFile=new File[3][];
        abd_lineageFile[0]= Arrays.stream(f).filter(ap).toArray(File[]::new);
        abd_lineageFile[1]= Arrays.stream(f).filter(bp).toArray(File[]::new);
        abd_lineageFile[2]= Arrays.stream(f).filter(dp).toArray(File[]::new);
        String[] outNames={"chr.Asubgenome.vcf", "chr.Bsubgenome.vcf", "chr.Dsubgenome.vcf"};
        BufferedWriter[] bws=new BufferedWriter[3];
        for (int i = 0; i < bws.length; i++) {
            bws[i]=IOUtils.getTextWriter(new File(outDir, outNames[i]).getAbsolutePath());
        }
        try {
            BufferedReader br;
            boolean ifFirst=true;
            StringBuilder sb;
            for (int i = 0; i < abd_lineageFile.length; i++) {
                for (int j = 0; j < abd_lineageFile[i].length; j++) {
                    sb=new StringBuilder(1000);
                    br=IOUtils.getTextGzipReader(abd_lineageFile[i][j].getAbsolutePath());
                    String line;
                    while ((line=br.readLine()).startsWith("##")){
                        sb.append(line);
                        sb.append("\n");
                    }
                    sb.append(line);
                    sb.append("\n");
                    if (ifFirst){
                        bws[i].write(sb.toString());
                        ifFirst=false;
                    }
                    while ((line=br.readLine())!=null){
                        bws[i].write(line);
                        bws[i].newLine();
                    }
                    br.close();
                }
                bws[i].flush();
                bws[i].close();
                ifFirst=true;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    public static void fastMergeVCFtoAB(String a_lineageVCF, String b_lineageVCF, String out_AB_lineageFile){
        try (BufferedReader bufferedReader= IOUtils.getTextReader(a_lineageVCF);
             BufferedReader bufferedReader1=IOUtils.getTextReader(b_lineageVCF);
             BufferedWriter bufferedWriter=IOUtils.getTextWriter(out_AB_lineageFile)) {
            String line;
            while ((line=bufferedReader.readLine()).startsWith("##")){
                bufferedWriter.write(line);
                bufferedWriter.newLine();
            }
            bufferedWriter.write(line);
            bufferedWriter.newLine();
            while ((line=bufferedReader.readLine())!=null){
                bufferedWriter.write(line);
                bufferedWriter.newLine();
            }
            while ((line=bufferedReader1.readLine()).startsWith("##")){}
            while ((line=bufferedReader1.readLine())!=null){
                bufferedWriter.write(line);
                bufferedWriter.newLine();
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * only support even, between chr001 and chr042, suitable for large VCF file
     * @param inputVcfDir inputVcfDir
     * @param outDir outDir
     */
    public static void fastMergeVCFtoChr(String inputVcfDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow()+ "start");
        long start= System.nanoTime();
        List<File> files=IOUtils.getFileListInDirEndsWith(inputVcfDir, "gz");
        Predicate<File> hidden=File::isHidden;
        File[] f=files.stream().filter(hidden.negate()).sorted().toArray(File[]::new);
        List<String> outChrs=WheatLineage.abdLineage();
        try{
            BufferedReader br1, br2;
            BufferedWriter bw;
            List<String> temp;
            String line;
            int vcfChr, vcfPos;
            String refChr;
            int refPos, chrNum;
            for (int i = 0; i < f.length; i=i+2) {
                br1=IOTool.getReader(f[i].getAbsolutePath());
                br2=IOTool.getReader(f[i+1].getAbsolutePath());
                chrNum=Integer.parseInt(f[i].getName().substring(3,6));
                bw= IOTool.getWriter(new File(outDir, "chr"+RefV1Utils.getChromosome(chrNum, 0)+".vcf.gz").getAbsolutePath());
                while ((line=br1.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br1.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    vcfChr=Integer.parseInt(temp.get(0));
                    vcfPos=Integer.parseInt(temp.get(1));
                    refChr= RefV1Utils.getChromosome(vcfChr, vcfPos);
                    refPos= RefV1Utils.getPosOnChromosome(vcfChr, vcfPos);
                    temp.set(0, refChr);
                    temp.set(1, String.valueOf(refPos));
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                }
                br1.close();
                while ((line=br2.readLine()).startsWith("##")){}
                while ((line=br2.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    vcfChr=Integer.parseInt(temp.get(0));
                    vcfPos=Integer.parseInt(temp.get(1));
                    refChr= RefV1Utils.getChromosome(vcfChr, vcfPos);
                    refPos= RefV1Utils.getPosOnChromosome(vcfChr, vcfPos);
                    temp.set(0, refChr);
                    temp.set(1, String.valueOf(refPos));
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                }
                br2.close();
                bw.flush();
                bw.close();
                System.out.println(new File(outDir, "chr"+outChrs.get(i/2)+".vcf").getName()+" completed in "
                        +Benchmark.getTimeSpanMinutes(start)+" minutes");
            }
            System.out.println(" all chromosomes completed in "+Benchmark.getTimeSpanHours(start)+" hours");
            System.out.println(DateTime.getDateTimeOfNow()+" end");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * split A.vcf B.vcf D.vcf to chr001.vcf chr002.vcf et.al.
     * @param subgenomeDir subgenomeDir
     * @param outDir outDir
     */
    public static void splitSubgenome(String subgenomeDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(subgenomeDir);
        IntStream.range(0, files.size()).forEach(e->splitSubgenome(files.get(e), outDir));
    }

    public static void splitSubgenome(File subgenomeFile, String outDir){
        Set<String> set=RowTableTool.getColumnSet(subgenomeFile.getAbsolutePath(), 0);
        List<String> chrs=new ArrayList<>(set);
        chrs.sort(Comparator.comparing(Integer::parseInt));
        Map<String, BufferedWriter> chrBufferedWriter=new HashMap<>();
        BufferedWriter bw;
        for (String s : chrs) {
            int chr = Integer.parseInt(s);
            bw = IOTool.getWriter(new File(outDir, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".txt"));
            chrBufferedWriter.put(s, bw);
        }
        try (BufferedReader bufferedReader = IOTool.getReader(subgenomeFile)) {
            String line;
            List<String> temp;
            String header=bufferedReader.readLine();
            for(Map.Entry<String, BufferedWriter> entry: chrBufferedWriter.entrySet()){
                entry.getValue().write(header);
                entry.getValue().newLine();
            }
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                bw=chrBufferedWriter.get(temp.get(0));
                bw.write(line);
                bw.newLine();
            }
            for(Map.Entry<String, BufferedWriter> entry: chrBufferedWriter.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Genotype for 0/0, 0/1, 1/1, or 0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles
     * @param vcfLine vcfLine
     * @return maf or -1, if 3 or more alt alleles exist
     */
    public static String calculateMaf(String vcfLine){
        double[] alleleFrequency=new double[3]; // 0 1 2
        String[] genotypeArray={"0/0", "0/1", "0/2", "1/1", "1/2", "2/2"};
        long[] genotypeCount=new long[genotypeArray.length];
        List<String> vcfLineList=PStringUtils.fastSplit(vcfLine);
        Predicate<String> p=str->str.startsWith("./.");
        List<String> vcfl=vcfLineList.stream().skip(9).filter(p.negate()).map(str->str.substring(0, 3)).collect(Collectors.toList());
        Map<String, Long> map=vcfl.stream().collect(Collectors.groupingBy(s-> s, Collectors.counting()));
        List<String> keyList=new ArrayList<>(map.keySet());
        Collections.sort(keyList);
        int index;
        for (String s : keyList) {
            index = Arrays.binarySearch(genotypeArray, s);
            if (index < 0) {
                System.out.println(s + " genotype found");
                return "-1";
//                System.out.println("Only supports two or three alleles, program quit");
//                System.exit(1);
            }
            genotypeCount[index] = map.get(s);
        }
        double sum=Arrays.stream(genotypeCount).sum();
        alleleFrequency[0]=(genotypeCount[0]*2+genotypeCount[1]+genotypeCount[2])/(sum*2);
        alleleFrequency[1]=(genotypeCount[1]+genotypeCount[3]*2+genotypeCount[4])/(sum*2);
        alleleFrequency[2]=(genotypeCount[2]+genotypeCount[4]+genotypeCount[5]*2)/(sum*2);
        Arrays.sort(alleleFrequency);
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setGroupingUsed(false);
        numberFormat.setMaximumFractionDigits(5);
        return numberFormat.format(alleleFrequency[1]);
    }

    /**
     * 位点为二等位时，计算Daf, 不包括三等位、四等位等
     * @param maf maf
     * @param majorAllele majorAllele
     * @param minorAllele minorAllele
     * @param ancestralAllele ancestralAllele
     * @return daf or -1 if two or more alt alleles exist
     */
    public static double calculateDaf(double maf, String majorAllele, String minorAllele, String ancestralAllele){
        double daf=-1;
        if (majorAllele.equals(ancestralAllele)){
            daf=maf;
        }else if (minorAllele.equals(ancestralAllele)){
            daf=1-maf;
        }
        return daf;
    }

    /**
     *
     * @param vcfDir vcfDir
     * @param subsetFileDir subsetFileDir
     * @param rate 0.001
     * @param numThreads 36
     */
    public static void getSubSetVcfFromDir(String vcfDir, String subsetFileDir, double rate, int numThreads){
        long totalStart=System.nanoTime();
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] f2=Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("vcf.*$", "subset"+rate+".vcf"))
                .toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, numThreads);
        for (int[] index : indices) {
            Integer[] subLibIndices = new Integer[index[1] - index[0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = index[0] + j;
            }
            List<Integer> integerList = Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e -> {
                BufferedReader br;
                BufferedWriter bw;
                long start = System.nanoTime();
                try {
                    if (f1[e].getName().endsWith("vcf")) {
                        br = IOUtils.getTextReader(f1[e].getAbsolutePath());
                        bw = IOUtils.getTextWriter(new File(subsetFileDir, f2[e]).getAbsolutePath());
                    } else {
                        br = IOUtils.getTextGzipReader(f1[e].getAbsolutePath());
                        bw = IOUtils.getTextWriter(new File(subsetFileDir, f2[e]).getAbsolutePath());
                    }
                    String line;
                    double r;
                    int count = 0;
                    int total = 0;
                    while ((line = br.readLine()).startsWith("##")) {
                        bw.write(line);
                        bw.newLine();
                    }
                    bw.write(line);
                    bw.newLine();
                    while ((line = br.readLine()) != null) {
                        total++;
                        r = Math.random();
                        if (r > rate) continue;
                        count++;
                        bw.write(line);
                        bw.newLine();
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("samping " + count + "(" + total + ") row from "
                            + f1[e].getName() + " into " + new File(subsetFileDir, f2[e]).getName() + " in "
                            + Benchmark.getTimeSpanMinutes(start) + " minutes");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            });
        }
        System.out.println("samping "+vcfDir+" into "+subsetFileDir+" is completed in "+Benchmark.getTimeSpanHours(totalStart)+" hours");
    }

    public static double calculateR2(List<String> vcfLine1, List<String> vcfLine2){
        String[] genotype1=vcfLine1.stream().skip(9).map(str->str.substring(0, 3)).toArray(String[]::new);
        String[] genotype2=vcfLine2.stream().skip(9).map(str->str.substring(0, 3)).toArray(String[]::new);
        return calculateR2(genotype1, genotype2);
    }

    /**
     * calculate r2 from genotype
     * note: only support binary allele, and genotype must be one of the following, "0/0", "1/1", "0/1", "./."
     * @param genotype1 genotype1
     * @param genotype2 genotype2
     * @return r2
     */
    public static double calculateR2(String[] genotype1, String[] genotype2){
        if (genotype1.length != genotype2.length){
            System.out.println("please check your input genotype array, its length is not same ");
            System.exit(1);
        }
        String[] genoypes={"0/0", "1/1", "0/1", "./."};
        double[] values={0, 1};
        Map<String, Double> genotypeValue=new HashMap<>();
        for (int i = 0; i < values.length; i++) {
            genotypeValue.put(genoypes[i], values[i]);
        }
        TDoubleArrayList array1=new TDoubleArrayList();
        TDoubleArrayList array2=new TDoubleArrayList();
        for (int i = 0; i < genotype1.length; i++) {
            if (genotype1[i].equals("./.")) continue;
            if (genotype2[i].equals("./.")) continue;
            if (genotype1[i].equals("0/1") && genotype2[i].equals("0/1")) continue;
            if (genotype1[i].equals("0/1")){
                array1.add(0);
                array1.add(1);
                array2.add(genotypeValue.get(genotype2[i]));
                array2.add(genotypeValue.get(genotype2[i]));
                continue;
            }
            if (genotype2[i].equals("0/1")){
                array2.add(0);
                array2.add(1);
                array1.add(genotypeValue.get(genotype1[i]));
                array1.add(genotypeValue.get(genotype1[i]));
                continue;
            }
            array1.add(genotypeValue.get(genotype1[i]));
            array1.add(genotypeValue.get(genotype1[i]));
            array2.add(genotypeValue.get(genotype2[i]));
            array2.add(genotypeValue.get(genotype2[i]));
        }
        if (array1.size()<20) return Double.NaN;
        PearsonsCorrelation pearsonsCorrelation=new PearsonsCorrelation();
        double r=pearsonsCorrelation.correlation(array1.toArray(), array2.toArray());
        return Math.pow(r, 2);
    }

    /**
     * recode 1A, 1B, 1D, 2A, 2B, 2D, ... to 1, 2, 3, 4, 5, 6, ...
     * @param vcfFile vcfFile
     * @param outFile outFile
     */
    public static void recode(File vcfFile, File outFile){
        int[] integers= IntStream.range(1,8).toArray();
        String[] abd={"A","B","D"};
        Map<String, Integer> chrIDmap=new HashMap<>();
        int count=0;
        for (int integer : integers) {
            for (String s : abd) {
                count++;
                chrIDmap.put(integer + s, count);
            }
        }
        try (BufferedReader br = IOUtils.getTextReader(vcfFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            bw.write(line);
            bw.newLine();
            List<String> lineList;
            String key;
            Integer value;
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                key=lineList.get(0);
                value=chrIDmap.get(key);
                lineList.set(0, String.valueOf(value));
                bw.write(String.join("\t", lineList));
                bw.newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * parallel
     * @param vcfDir vcfDir
     * @param outDir outDir
     * @param numThreads numThreads
     */
    public static void recode(String vcfDir, String outDir, int numThreads){
        long start=System.nanoTime();
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] outNames= Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("vcf$", "recoded.vcf"))
                .toArray(String[]::new);
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, numThreads);
        for (int[] ints : indices) {
            Integer[] subLibIndices = new Integer[ints[1] - ints[0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = ints[0] + j;
            }
            List<Integer> integerList = Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(index -> {
                        long start1 = System.nanoTime();
                        VCF.recode(f1[index], new File(outDir, outNames[index]));
                        System.out.println(new File(outDir, outNames[index]).getName() + " completed in " + Benchmark.getTimeSpanMinutes(start1) + " minutes");
                    });
        }
        System.out.println("completed in "+Benchmark.getTimeSpanHours(start)+" hours");
    }

    /**
     *
     * @param vcfInputFile vcfInputFile
     * @param columnIndex columnIndex
     * @return chr pos columnIndex, no header
     */
    public static Table<Integer, Integer, String> getTable(String vcfInputFile, int columnIndex){
        BufferedReader bufferedReader;
        Table<Integer, Integer, String> table= HashBasedTable.create();
        if (vcfInputFile.endsWith("gz")){
            bufferedReader= IOUtils.getTextGzipReader(vcfInputFile);
        }else {
            bufferedReader=IOUtils.getTextReader(vcfInputFile);
        }
        try {
            String line;
            List<String> temp;
            String[] tem;
            while ((line=bufferedReader.readLine()).startsWith("##")){}
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                int chr=Integer.parseInt(temp.get(0));
                int pos=Integer.parseInt(temp.get(1));
                tem=StringUtils.split(temp.get(columnIndex), ",");
                if (tem.length>1) continue;
                table.put(chr, pos, tem[0]);
            }
            bufferedReader.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        return table;
    }

    public static List<String> getColumnList(File vcfFile, int columnIndex){
        BufferedReader bufferedReader;
        if (vcfFile.toString().endsWith("gz")){
            bufferedReader=IOUtils.getTextGzipReader(vcfFile.getAbsolutePath());
        }else {
            bufferedReader=IOUtils.getTextReader(vcfFile.getAbsolutePath());
        }
        String line;
        List<String> temp;
        List<String> res = new ArrayList<>(1000);
        try {
            while ((line= bufferedReader.readLine()).startsWith("##")){}
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                res.add(temp.get(columnIndex));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    public static List<String>[] getColumnList(String vcfDir, int columnIndex){
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> predicate=File::isHidden;
        File[] f=Arrays.stream(files).filter(predicate.negate()).toArray(File[]::new);
        List<String>[] res=new List[f.length];
        List<String> pos;
        for (int i = 0; i < f.length; i++) {
            pos=VCF.getColumnList(f[i], columnIndex);
            res[i]=pos;
        }
        return res;
    }

    public static double getHeterozygousProportionBySite(List<String> temp){
        List<String> genotypeList=temp.subList(9, temp.size());
        double heteroCount=0, totalCount=0;
        for (String s : genotypeList) {
            if (s.startsWith("./.")) continue;
            totalCount++;
            if (s.startsWith("0/0")) continue;
            if (s.startsWith("1/1")) continue;
            heteroCount++;
        }
        return heteroCount/totalCount;
    }

    /**
     *
     * @return 返回VCF文件所有taxa的总数
     */
    public int getNumberOfTaxa(){
        return header.size()-9;
    }

    public List<List<String>> getData(){return this.data;}

    public List<String> getHeader(){
        return header;
    }

    /**
     *
     * @param columnIndex columnIndex
     * @return chr pos columnIndex, no header
     */
    public Table<Integer, Integer, String> getTable(int columnIndex){
        Table<Integer, Integer, String> table= HashBasedTable.create();
        List<List<String>> data=this.getData();
        String target;
        for (List<String> datum : data) {
            int chr = Integer.parseInt(datum.get(0));
            int pos = Integer.parseInt(datum.get(1));
            target = datum.get(columnIndex);
            table.put(chr, pos, target);
        }
        return table;
    }

    /**
     *
     * @param columnIndex columnIndex
     * @return the specified column
     */
    public List<String> getColumnList(int columnIndex){
        List<String> res = new ArrayList<>();
        for (int i = 0; i < data.size(); i++) {
            res.add(data.get(i).get(columnIndex));
        }
        return res;
    }

    /**
     * Removes all of the rows of this rowTable that satisfy the given predicate
     * @param p 针对行进行过滤的函数式接口
     * @return true, if any row were removed
     */
    public boolean removeIf(Predicate<List<String>> p){
        List<List<String>> data=this.getData();
        return data.removeIf(p);
    }

    public VCF subsetVCF(Set<String> subsetTaxonSet){
        List<String> headList=this.getHeader();
        List<String> subsetTaxonList=new ArrayList<>(subsetTaxonSet);
        TIntArrayList subsetTaxonIndexList=new TIntArrayList();
        int index;
        for (String s : subsetTaxonList) {
            index = headList.indexOf(s);
            subsetTaxonIndexList.add(index);
        }
        subsetTaxonIndexList.add(IntStream.range(0, 9).toArray());
        subsetTaxonIndexList.sort();
        List<String> subHeader=new ArrayList<>();
        List<List<String>> subData=new ArrayList<>();
        for (int i = 0; i < headList.size(); i++) {
            if (!subsetTaxonIndexList.contains(i)) continue;
            subHeader.add(headList.get(i));
        }
        List<String> lineList;
        List<String> subLineList;
        Predicate<String> p=str->str.startsWith("./.");
        Predicate<String> only0_0=str->str.equals("0/0");
        Predicate<String> pp=only0_0.or(str->str.equals("1/1"));
        for (int i = 0; i < this.getData().size(); i++) {
            lineList=this.getData().get(i);
            subLineList=new ArrayList<>();
            for (int j = 0; j < lineList.size(); j++) {
                if (!subsetTaxonIndexList.contains(j)) continue;
                subLineList.add(lineList.get(j));
            }
            if(subLineList.stream().skip(9).filter(p.negate()).map(str->str.substring(0, 3)).allMatch(pp)) continue;
            subData.add(subLineList);
        }
        return new VCF(this.meta, subHeader, subData);
    }

    /**
     *
     * @param heterozygousProportionBySite heterozygousProportionBySite
     */
    public void filterOnHeterozygousProportionBySite(double heterozygousProportionBySite){
        List<List<String>> data=this.data;
        List<List<String>> newdata=new ArrayList<>();
        for (List<String> datum : data) {
            if (VCF.getHeterozygousProportionBySite(datum) > heterozygousProportionBySite) continue;
            newdata.add(datum);
        }
        this.data=newdata;
    }


    /**
     * 根据VCF文件的CHR和POS对该对象排序
     */
    public void sort(){
        if (StringTool.isNumeric(data.get(0).get(0))){
            Comparator<List<String>> c=Comparator.comparing(e->Integer.parseInt(e.get(0)));
            c=c.thenComparing(e->Integer.parseInt(e.get(1)));
            this.sort(c);
        }else {
            Comparator<List<String>> c=Comparator.comparing(e->e.get(0));
            c=c.thenComparing(e->Integer.parseInt(e.get(1)));
            this.sort(c);
        }
    }

    /**
     *
     * @param comparator comparator
     */
    public void sort(Comparator<List<String>> comparator){
        data.sort(comparator);
    }

    public void removeVarianceWithHighGenotypeMissingRate(double genotypeMissingRate){
        double threshold= this.getNumberOfTaxa()*(1-genotypeMissingRate);
        int[] genotyedTaxonNum=this.getGenotypedTaxonNum();
        int[] reservedLine=IntStream.range(0,genotyedTaxonNum.length).filter(e->genotyedTaxonNum[e]>=threshold).toArray();
        List<List<String>> l=new ArrayList<>(reservedLine.length);
        for (int value : reservedLine) {
            l.add(data.get(value));
        }
        this.data=l;
    }

    /**
     *
     * @return 返回所有位点的基因型缺失比例
     */
    public double[] getGenotypeMissingRate(){
        double[] genotypeMissingRate= new double[data.size()];
        for(int i=0,size=data.size();i<size;i++){
            String info=data.get(i).stream().limit(8).collect(Collectors.toList()).get(7);
            String ns=PStringUtils.fastSplit(info, ";").get(2);
            double nsValue=Integer.parseInt(PStringUtils.fastSplit(ns, "=").get(1));
            genotypeMissingRate[i]=(this.getNumberOfTaxa()-nsValue)/this.getNumberOfTaxa();
        }
        return genotypeMissingRate;
    }

    /**
     * 统计每个variant中具有基因型的taxon数目
     * @return  the taxon number which having genotype
     */
    public int[] getGenotypedTaxonNum(){
        int[] numberOfSamplesWithAllele =new int[data.size()];
        for(int i=0,size=data.size();i<size;i++){
            String info=data.get(i).stream().limit(8).collect(Collectors.toList()).get(7);
            String ns=PStringUtils.fastSplit(info, ";").get(2);
            int nsValue=Integer.parseInt(PStringUtils.fastSplit(ns, "=").get(1));
            numberOfSamplesWithAllele[i]=nsValue;
        }
        return numberOfSamplesWithAllele;
    }

    /**
     * 将指定的VCF与这个VCF进行融合（两个VCF的meta与header必须相同）
     * @param vcf 指定的VCF对象
     * @return 当VCF对象域data在调用后发生变化时返回true
     */
    public boolean addVCF(VCF vcf){
        return this.data.addAll(vcf.data);
    }

    /**
     * change vcfChr and vcfPos to refChr and refPos
     */
    public void changeToRefChr(){
        List<List<String>> data=this.data;
        List<List<String>> res=new ArrayList<>();
        for (int i = 0; i < data.size(); i++) {
            res.add(new ArrayList<>());
        }
        List<String> line;
        short vcfChr;
        int vcfPos;
        String refChr;
        int refPos;
        for (int i = 0; i < data.size(); i++) {
            line=data.get(i);
            vcfChr=Short.parseShort(line.get(0));
            vcfPos=Integer.parseInt(line.get(1));
            refChr= RefV1Utils.getChromosome(vcfChr, vcfPos);
            refPos= RefV1Utils.getPosOnChromosome(vcfChr,vcfPos);
            line.set(0, refChr);
            line.set(1, String.valueOf(refPos));
            res.get(i).addAll(line);
        }
        this.data=res;
    }

    /**
     * change refChr and reffPos to vcfChr and vcfPos
     */
    public void changeToVCFChr(){
        List<List<String>> data=this.data;
        List<List<String>> res=new ArrayList<>();
        for (int i = 0; i < data.size(); i++) {
            res.add(new ArrayList<>());
        }
        List<String> line;
        String refChr;
        int refPos;
        int vcfChr;
        int vcfPos;
        for (int i = 0; i < data.size(); i++) {
            line=data.get(i);
            refChr=line.get(0);
            refPos=Integer.parseInt(line.get(1));
            vcfChr=RefV1Utils.getChrID(refChr, refPos);
            vcfPos=RefV1Utils.getPosOnChrID(refChr, refPos);
            line.set(0, String.valueOf(vcfChr));
            line.set(1, String.valueOf(vcfPos));
            res.get(i).addAll(line);
        }
        this.data=res;
    }

    public void writeR2(String outFile){
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write("Chr1\tPos1\tChr2\tPos2\tr2\n");
            List<List<String>> data=this.getData();
            double r2;
            StringBuilder sb;
            Iterator<int[]> iterator = CombinatoricsUtils.combinationsIterator(data.size(), 2);
            int[] combinationIndex;
            while (iterator.hasNext()) {
                combinationIndex = iterator.next();
                r2=VCF.calculateR2(data.get(combinationIndex[0]), data.get(combinationIndex[1]));
                sb=new StringBuilder(30);
                sb.append(data.get(combinationIndex[0]).get(0)).append("\t");
                sb.append(data.get(combinationIndex[0]).get(1)).append("\t");
                sb.append(data.get(combinationIndex[1]).get(0)).append("\t");
                sb.append(data.get(combinationIndex[1]).get(1)).append("\t");
                sb.append(r2);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }


    }

    public void write(String outFile){
        this.sort();
        BufferedWriter bw;
        if (outFile.endsWith(".vcf")){
            bw=IOUtils.getTextWriter(outFile);
        }else {
            bw=IOUtils.getTextGzipWriter(outFile);
        }
        try{
            bw.write(this.meta);
            StringBuilder sb=new StringBuilder();
            for (String s : header) {
                sb.append(s);
                sb.append("\t");
            }
            sb.deleteCharAt(sb.length()-1).append("\n");
            bw.write(sb.toString());
            sb=new StringBuilder();
            bw.write(sb.toString());
            for (List<String> datum : data) {
                for (int j = 0; j < datum.size(); j++) {
                    sb = sb.append(datum.get(j)).append("\t");
                }
                sb.deleteCharAt(sb.length() - 1).append("\n");
                bw.write(sb.toString());
                sb = new StringBuilder();
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void write(String outputDir, String fileName){
       this.write(outputDir+"/"+fileName);
    }

    /**
     * 将VCF按照"chr000.vcf, chr001.vcf"的形式进行输出
     * @param outputDir 输出目录
     */
    public void writeVcfToSplitedChrID(String outputDir){
        this.changeToVCFChr();
        this.sort();
        List<Integer> chrList=data.stream().flatMap(e->e.stream().limit(1)).mapToInt(Integer::valueOf).boxed()
                              .distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
        Map<Integer, BufferedWriter> strToBufferedWriterMap=new HashMap<>();
        Integer key;
        BufferedWriter value;
        for (Integer integer : chrList) {
            key = integer;
            value = IOUtils.getTextWriter(outputDir + "/chr" + PStringUtils.getNDigitNumber(3, key) + ".vcf");
            strToBufferedWriterMap.put(key, value);
        }
        StringBuilder sb=new StringBuilder();
        sb.append(meta);
        for (String s : header) {
            sb.append(s);
            sb.append("\t");
        }
        sb.deleteCharAt(sb.length()-1).append("\n");
        try{
            for(Map.Entry<Integer, BufferedWriter> entry:strToBufferedWriterMap.entrySet()){
                entry.getValue().write(sb.toString());
            }
            for (List<String> datum : data) {
                for (int j = 0; j < datum.size(); j++) {
                    key = Integer.valueOf(datum.get(0));
                    value = strToBufferedWriterMap.get(key);
                    value.write(String.join("\t", datum));
                    value.newLine();
                    break;
                }
            }
            for(Map.Entry<Integer, BufferedWriter> entry:strToBufferedWriterMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 将VCF按照"chr1A, chr1B.vcf"的形式进行输出
     * @param outputDir outputDir
     */
    public void writeVcfToSplitedChr(String outputDir){
        if (StringTool.isNumeric(this.data.get(0).get(0))){
            this.changeToRefChr();
        }
        this.sort();
        List<String> chrList=data.stream().flatMap(e->e.stream().limit(1)).distinct().sorted().collect(Collectors.toList());
        Map<String, BufferedWriter> strToBufferedWriterMap=new HashMap<>();
        String key;
        BufferedWriter value;
        for (String s : chrList) {
            key = s;
            value = IOUtils.getTextWriter(outputDir + "/chr" + key + ".vcf");
            strToBufferedWriterMap.put(key, value);
        }
        StringBuilder sb=new StringBuilder();
        sb.append(meta);
        for (String s : header) {
            sb.append(s);
            sb.append("\t");
        }
        sb.deleteCharAt(sb.length()-1).append("\n");
        try{
            for(Map.Entry<String, BufferedWriter> entry:strToBufferedWriterMap.entrySet()){
                entry.getValue().write(sb.toString());
            }
            for (List<String> datum : data) {
                for (int j = 0; j < datum.size(); j++) {
                    key = datum.get(0);
                    value = strToBufferedWriterMap.get(key);
                    value.write(String.join("\t", datum));
                    value.newLine();
                    break;
                }
            }
            for(Map.Entry<String, BufferedWriter> entry:strToBufferedWriterMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void midFormatToVCF(String inputDir, String outDir){
        List<File> fileList= IOUtils.getVisibleFileListInDir(inputDir);
        String[] outFiles= fileList.stream().map(File::getName).map(f->f.replaceAll(".mid.gz", ".vcf.gz")).toArray(String[]::new);
        IntStream.range(0, fileList.size()).parallel().forEach(e->midFormatToVCF(fileList.get(e), new File(outDir,
                outFiles[e])));
    }

    private static void midFormatToVCF(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, genotypeTemp;
            line=br.readLine();
            temp=PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(8, temp.size());
            StringBuilder sb=new StringBuilder();
            sb.setLength(0);
            sb.append("##fileformat=VCFv4.2\n##fileDate=20210122\n");
            sb.append("##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
            sb.append("##INFO=<ID=GN,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            sb.append("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">\n");
            sb.append("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            sb.append("##Species=Rice\n");
            sb.append("##ReferenceGenome=IRGSP 4.0\n");
            sb.append("##Reference=\"A map of rice genome variation reveals the origin of cultivated rice\"");
            bw.write(sb.toString());
            bw.newLine();
            sb.setLength(0);
            sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
            sb.append(String.join("\t", taxonNames));
            bw.write(sb.toString());
            bw.newLine();
            String chr, position, refBase, altBase;
            int genotypedNum, alleleCount;
            double maf;
            NumberFormat numberFormt = NumberFormat.getInstance();
            numberFormt.setMaximumFractionDigits(5);
            numberFormt.setGroupingUsed(false);
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chr=temp.get(0);
                position=temp.get(1);
                refBase=temp.get(2);
                altBase=temp.get(3);
                genotypedNum=Integer.parseInt(temp.get(4));
                maf=Double.parseDouble(temp.get(5));
                alleleCount=Integer.parseInt(temp.get(7));
                sb.setLength(0);
                sb.append(chr).append("\t").append(position).append("\t").append(chr).append("-").append(position).append("\t");
                sb.append(refBase).append("\t").append(altBase).append("\t");
                sb.append(".").append("\t").append(".").append("\t");
                sb.append("GN=").append(genotypedNum).append(";");
                sb.append("AC=").append(alleleCount).append(";");
                sb.append("MAF=").append(numberFormt.format(maf)).append("\t");
                sb.append("GT").append("\t");
                genotypeTemp=temp.subList(8, temp.size());
                for (String genotype: genotypeTemp){
                    if (genotype.equals(refBase)){
                        sb.append("0/0").append("\t");
                    }else if (genotype.equals(altBase)){
                        sb.append("1/1").append("\t");
                    }else {
                        sb.append("./.").append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param taxonListFile containing a list of individuals to include in subsequent analysis
     * @param vcfDir input VCF dir
     * @param outDir out dir
     */
    public static void extractPopulation(String taxonListFile, String vcfDir, String outDir){
        try {
            List<String> taxonList= Files.readAllLines(Paths.get(taxonListFile));
            List<File> files= IOUtils.getVisibleFileListInDir(vcfDir);
            String[] outFileNames=files.stream().map(File::getName).map(str->
                    str.replaceAll(".vcf.gz", ".subset.vcf.gz")).toArray(String[]::new);
            IntStream.range(0, files.size()).parallel().forEach(e->extractPopulation(taxonList, files.get(e),
                    new File(outDir, outFileNames[e])));
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void extractPopulation(List<String> taxonList, File vcfFile, File outFile){
        try (BufferedReader br = IOTool.getReader(vcfFile);
             BufferedWriter bw =IOTool.getWriter(outFile)){
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            StringBuilder sb=new StringBuilder();
            temp=PStringUtils.fastSplit(line);
            TIntArrayList taxonIndexList=new TIntArrayList();
            for (String taxon:taxonList){
                taxonIndexList.add(temp.indexOf(taxon));
            }
            taxonIndexList.sort();
            sb.append(String.join("\t", temp.subList(0,9))).append("\t");
            for (int i = 0; i < taxonIndexList.size(); i++) {
                sb.append(temp.get(taxonIndexList.get(i))).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            String genotype;
            Set<String> genotypeSet;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                genotypeSet=new HashSet<>();
                for (int i = 0; i < taxonIndexList.size(); i++) {
                    genotype=temp.get(taxonIndexList.get(i));
                    if (genotype.equals("./.")) continue;
                    genotypeSet.add(genotype);
                }
                if (genotypeSet.size()<2) continue;
                sb.setLength(0);
                sb.append(String.join("\t", temp.subList(0, 9))).append("\t");
                for (int i = 0; i < taxonIndexList.size(); i++) {
                    sb.append(temp.get(taxonIndexList.get(i))).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
