package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import utils.*;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class VCF {

    private String meta;
    private List<String> header;
    private List<List<String>> data;
    public static Map<Integer,String> chrToChrMap =VCF.getchrToChrMap();

    public VCF(String inputFile){
        this.initilize(inputFile);
    }

    public VCF(File inputFile){
        this.initilize(inputFile.getAbsolutePath());
    }

    private void initilize(String inputFile){
        BufferedReader br;
        try{
            if (inputFile.endsWith(".vcf")){
                br=IOUtils.getTextReader(inputFile);
            }else {
                br= IOUtils.getTextGzipReader(inputFile);
            }
            StringBuilder sb=new StringBuilder();
            String temp;
            List<List<String>> lists=new ArrayList<>();
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

    private static Map<Integer, String> getchrToChrMap(){
        Map<Integer,String> chrToChrMap=new HashMap<>();
        List<Integer> numOfChr=IntStream.range(1,43).boxed().collect(Collectors.toList());
        List<Integer> int1_7= IntStream.range(1,8).boxed().collect(Collectors.toList());
        List<Integer> chrList=new ArrayList<>();
        for(int i=0;i<6;i++){
            chrList.addAll(int1_7);
        }
        Collections.sort(chrList);
        String abd=String.join("", Collections.nCopies(7,"AABBDD"));
        chrToChrMap.put(0, "Un");
        for(int i=0;i<numOfChr.size();i++){
            chrToChrMap.put(numOfChr.get(i), String.valueOf(chrList.get(i))+abd.charAt(i));
        }
        chrToChrMap.put(43, "Mit");
        chrToChrMap.put(44, "Chl");
        return chrToChrMap;
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
     * @param inputVcfDir
     * @param outDir
     * @param chrConvertionRule
     */
    public static void mergeVCFtoLineage(String inputVcfDir, String outDir, ChrConvertionRule chrConvertionRule){
        File[] files=new File(inputVcfDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
        File[] f=Arrays.stream(files).filter(p).sorted().toArray(File[]::new);
        int[] chrIDArray=Arrays.stream(f).filter(p).map(File::getName).map(str->str.substring(3,6))
                .mapToInt(Integer::parseInt).toArray();
        int[][] lineage=new int[3][];
        lineage[0]=WheatLineage.getWheatLineageOf(WheatLineage.A);
        lineage[1]=WheatLineage.getWheatLineageOf(WheatLineage.B);
        lineage[2]=WheatLineage.getWheatLineageOf(WheatLineage.D);
        TIntArrayList[] indexArray=new TIntArrayList[3];
        for (int i = 0; i < 3; i++) {
            indexArray[i]=new TIntArrayList();
        }
        int index=Integer.MIN_VALUE;
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
            vcfArray[i].changeToRefChr(chrConvertionRule);
            vcfArray[i].write(outDir, outNames[i]);
        }
    }

    /**
     * merge chr001 chr002 chr003 chr004, ... to chr1A chr1B chr1D, ...
     * @param inputVcfDir
     * @param outDir
     * @param chrConvertionRule
     */
    public static void mergeVCFtoChr(String inputVcfDir, String outDir, ChrConvertionRule chrConvertionRule){
        File[] files=new File(inputVcfDir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> p= hidden.negate().and(f->f.getName().toLowerCase().startsWith("chr"));
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
        VCF vcf=null;
        String[] outNameArray= Arrays.stream(f).map(File::getName).map(str->str.substring(6)).toArray(String[]::new);
        String outName=null;
        for (int i = 0; i < indexList.size(); i=i+2) {
            vcf=new VCF(f[indexList.get(i)]);
            vcf.addVCF(new VCF(f[indexList.get(i+1)]));
            vcf.changeToRefChr(chrConvertionRule);
            outName=VCF.getchrToChrMap().get(chrIDArray[indexList.get(i)]);
            vcf.write(outDir, "chr"+outName+outNameArray[indexList.get(i)]);
        }
    }

    /**
     * Genotype for 0/0, 0/1, 1/1, or 0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles
     * @param vcfLine
     * @return maf or -1, if 3 or more alt alleles exist
     */
    public static double calculateMaf(String vcfLine){
        double[] alleleFrequency=new double[3]; // 0 1 2
        String[] genotypeArray={"0/0", "0/1", "0/2", "1/1", "1/2", "2/2"};
        long[] genotypeCount=new long[genotypeArray.length];
        for (int i = 0; i < genotypeCount.length; i++) {
            genotypeCount[i]=0;
        }
        List<String> vcfLineList=PStringUtils.fastSplit(vcfLine);
        Predicate<String> p=str->str.startsWith("./.");
        List<String> vcfl=vcfLineList.stream().skip(9).filter(p.negate()).map(str->str.substring(0, 3)).collect(Collectors.toList());
        Map<String, Long> map=vcfl.stream().collect(Collectors.groupingBy(s-> s, Collectors.counting()));
        List<String> keyList=new ArrayList<>(map.keySet());
        Collections.sort(keyList);
        int index=Integer.MIN_VALUE;
        for (int i = 0; i < keyList.size(); i++) {
            index=Arrays.binarySearch(genotypeArray, keyList.get(i));
            if (index < 0){
                System.out.println(keyList.get(i)+" genotype found");
                return -1d;
//                System.out.println("Only supports two or three alleles, program quit");
//                System.exit(1);
            }
            genotypeCount[index]=map.get(keyList.get(i));
        }
        double sum=Arrays.stream(genotypeCount).sum();
        alleleFrequency[0]=(genotypeCount[0]*2+genotypeCount[1]+genotypeCount[2])/(sum*2);
        alleleFrequency[1]=(genotypeCount[1]+genotypeCount[3]*2+genotypeCount[4])/(sum*2);
        alleleFrequency[2]=(genotypeCount[2]+genotypeCount[4]+genotypeCount[5]*2)/(sum*2);
        Arrays.sort(alleleFrequency);
        return NumberTool.format(alleleFrequency[1], 5);
    }

    /**
     * 位点为二等位时，计算Daf, 不包括三等位、四等位等
     * @param maf
     * @param majorAllele
     * @param minorAllele
     * @param ancestralAllele
     * @return daf or -1 if two or more alt alleles exist
     */
    public static double calculateDaf(double maf, String majorAllele, String minorAllele, String ancestralAllele){
        double daf=-1;
        if (majorAllele.equals(ancestralAllele)){
            daf=maf;
        }else if (minorAllele.equals(ancestralAllele)){
            daf=1-maf;
        }else {
            daf=-1;
        }
        return daf;
    }

    /**
     *
     * @param vcfDir
     * @param subsetFileDir
     * @param rate 0.001
     * @param numThreads 36
     */
    public static void getSubSetVcfFromFile(String vcfDir, String subsetFileDir, double rate, int numThreads){
        long totalStart=System.nanoTime();
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] f2=Arrays.stream(f1).map(File::getName).map(str->str.replaceAll(".vcf", ".subset.vcf"))
                .toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e->{
                BufferedReader br;
                BufferedWriter bw;
                long start=System.nanoTime();
                try {
                    if (f1[e].getName().endsWith("vcf")){
                        br=IOUtils.getTextReader(f1[e].getAbsolutePath());
                        bw=IOUtils.getTextWriter(new File(subsetFileDir, f2[e]).getAbsolutePath());
                    }else {
                        br=IOUtils.getTextGzipReader(f1[e].getAbsolutePath());
                        bw=IOUtils.getTextGzipWriter(new File(subsetFileDir, f2[e]).getAbsolutePath());
                    }
                    String line;
                    double r=-1d;
                    int count=0;
                    int total=0;
                    while ((line=br.readLine()).startsWith("##")){
                        bw.write(line);
                        bw.newLine();
                    }
                    bw.write(line);
                    bw.newLine();
                    while ((line=br.readLine())!=null){
                        total++;
                        r=Math.random();
                        if (r > rate) continue;
                        count++;
                        bw.write(line);
                        bw.newLine();
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("samping "+count+"("+total+") row from "
                            +f1[e].getName()+" into "+new File(subsetFileDir, f2[e]).getName()+" in "
                            +Benchmark.getTimeSpanMinutes(start)+" minutes");
                } catch (IOException ex) {
                    ex.printStackTrace();
                }
            });
        }
        System.out.println("samping "+vcfDir+" into "+subsetFileDir+" is completed in "+Benchmark.getTimeSpanHours(totalStart)+" hours");
    }

    /**
     *
     * @return 返回VCF文件所有taxa的总数
     */
    public int getNumberOfTaxa(){
        return header.size()-9;
    }

    public List<List<String>> getData(){return this.data;}

    /**
     * 根据VCF文件的CHR和POS对该对象排序
     */
    public void sort(){
        if (StringTool.isNumeric(data.get(0).get(0))){
            Comparator<List<String>> c=Comparator.comparing(e->Integer.valueOf(e.get(0)));
            data.sort(c.thenComparing(e->Integer.valueOf(e.get(1))));
        }else {
            Comparator<List<String>> c=Comparator.comparing(e->e.get(0));
            data.sort(c.thenComparing(e->Integer.valueOf(e.get(1))));
        }
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
            int nsValue=Integer.valueOf(PStringUtils.fastSplit(ns, "=").get(1));
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
            int nsValue=Integer.valueOf(PStringUtils.fastSplit(ns, "=").get(1));
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
     * @param chrConvertionRule
     */
    public void changeToRefChr(ChrConvertionRule chrConvertionRule){
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
            refChr= chrConvertionRule.getRefChrFromVCFChr(vcfChr);
            refPos= chrConvertionRule.getRefPosFromVCFChrPos(vcfChr, vcfPos);
            line.set(0, refChr);
            line.set(1, String.valueOf(refPos));
            res.get(i).addAll(line);
        }
        this.data=res;
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
            for(int i=0;i<header.size();i++){
                sb.append(header.get(i));
                sb.append("\t");
            }
            sb.deleteCharAt(sb.length()-1).append("\n");
            bw.write(sb.toString());
            sb=new StringBuilder();
            bw.write(sb.toString());
            for(int i=0;i<data.size();i++){
                for(int j=0;j<data.get(i).size();j++){
                    sb=sb.append(data.get(i).get(j)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1).append("\n");
                bw.write(sb.toString());
                sb=new StringBuilder();
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
    public void writeVcfToSplitedChr(String outputDir){
        this.sort();
        List<Integer> chrList=data.stream().flatMap(e->e.stream().limit(1)).mapToInt(Integer::valueOf).boxed()
                              .distinct().sorted().collect(Collectors.toCollection(ArrayList::new));
        Map<Integer, BufferedWriter> strToBufferedWriterMap=new HashMap<>();
        Integer key;
        BufferedWriter value;
        for(int i=0;i<chrList.size();i++){
            key=chrList.get(i);
            value=IOUtils.getTextWriter(outputDir+"/chr"+PStringUtils.getNDigitNumber(3, key)+".vcf");
            strToBufferedWriterMap.put(key, value);
        }
        StringBuilder sb=new StringBuilder();
        sb.append(meta);
        for(int i=0;i<header.size();i++){
            sb.append(header.get(i));
            sb.append("\t");
        }
        sb.deleteCharAt(sb.length()-1).append("\n");
        try{
            for(Map.Entry<Integer, BufferedWriter> entry:strToBufferedWriterMap.entrySet()){
                entry.getValue().write(sb.toString());
            }
            for(int i=0;i<data.size();i++){
                for(int j=0;j<data.get(i).size();j++){
                    key=Integer.valueOf(data.get(i).get(0));
                    value=strToBufferedWriterMap.get(key);
                    value.write(data.get(i).stream().collect(Collectors.joining("\t")));
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
}
