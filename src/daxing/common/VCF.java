package daxing.common;

import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class VCF {

    private String meta;
    private List<String> header;
    private List<List<String>> data;
    public static Map<Integer,String> numToChrMap=VCF.getNumToChrMap();

    public VCF(String inputFile){
        this.initilize(inputFile);
    }

    public VCF(File inputFile){
        this.initilize(inputFile.getAbsolutePath());
    }

    private void initilize(String inputFile){
        try(BufferedReader br=IOUtils.getTextReader(inputFile)){
            StringBuilder sb=new StringBuilder();
            String temp;
            while ((temp=br.readLine())!=null){
                if(temp.startsWith("##")){
                    sb.append(temp);
                    sb.append("\n");
                } else if(temp.startsWith("#")){
                    this.header=PStringUtils.fastSplit(temp);
                    break;
                }
            }
            this.meta =sb.toString();
            this.data=br.lines().map(PStringUtils::fastSplit).collect(Collectors.toList());
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private static Map<Integer, String> getNumToChrMap(){
        Map<Integer,String> numChrMap=new HashMap<>();
        List<Integer> numOfChr=IntStream.range(1,43).boxed().collect(Collectors.toList());
        List<Integer> int1_7= IntStream.range(1,8).boxed().collect(Collectors.toList());
        List<Integer> chrList=new ArrayList<>();
        for(int i=0;i<6;i++){
            chrList.addAll(int1_7);
        }
        Collections.sort(chrList);
        String abd=String.join("", Collections.nCopies(7,"AABBDD"));
        numChrMap.put(0, "ChrUn");
        for(int i=0;i<numOfChr.size();i++){
            numChrMap.put(numOfChr.get(i),"Chr"+chrList.get(i)+abd.charAt(i));
        }
        numChrMap.put(43, "Mit");
        numChrMap.put(44, "Chl");
        return numChrMap;
    }

    /**
     * 从VCF文件中随机抽取行，组成新文件
     * @param inputFile VCF输入文件
     * @param outFile 输出文件
     * @param numberOfRow 提取的行数
     */
    public static void extractRandomRowFromVCF(String inputFile, String outFile, Integer numberOfRow) {
        long start = System.nanoTime();
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            List<String> metaAndHeader= Files.newBufferedReader(Paths.get(inputFile))
                    .lines().limit(30)
                    .filter(index->index.startsWith("#"))
                    .collect(Collectors.toList());
            List<String> data=Files.newBufferedReader(Paths.get(inputFile))
                    .lines()
                    .filter(index -> (!index.startsWith("#")))
                    .collect(Collectors.toList());
            int[] randomIndex=RandomArray.getRandomNonrepetitionArray(numberOfRow,0,data.size());
            Arrays.sort(randomIndex);
            for(String str:metaAndHeader){
                bw.write(str);
                bw.newLine();
            }
            for(int i=0;i<randomIndex.length;i++){
                bw.write(data.get(randomIndex[i]));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
        System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanSeconds(start)) + " seconds");
    }

    /**
     * 从VCF文件中随机抽取行，在输入文件的目录下生成新的VCF文件（后缀是sample.vcf）
     * @param inputFile VCF输入文件
     * @param numberOfRow 提取的行数
     */
    public static void extractRandomRowFromVCF(String inputFile, Integer numberOfRow) {
        long start = System.nanoTime();
        String outFile=inputFile.replaceAll(".vcf$", "sample.vcf");
        try (BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            List<String> metaAndHeader=Files.newBufferedReader(Paths.get(inputFile))
                                            .lines().limit(30)
                                            .filter(e->e.startsWith("#"))
                                            .collect(Collectors.toList());
            List<String> data=Files.newBufferedReader(Paths.get(inputFile))
                                   .lines()
                                   .filter(index -> (!index.startsWith("#")))
                                   .collect(Collectors.toList());
            int[] randomIndex=RandomArray.getRandomNonrepetitionArray(numberOfRow,0,data.size());
            Arrays.sort(randomIndex);
            for(String str:metaAndHeader){
                bw.write(str);
                bw.newLine();
            }
            for(int i=0;i<randomIndex.length;i++){
                bw.write(data.get(randomIndex[i]));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
        System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanSeconds(start)) + " seconds");
    }

    /**
     * 根据VCF目录和染色体编号，返回对应染色体编号的各个VCF文件路径
     * @param vcfDir VCF目录
     * @param chrNumber 打断后的染色体编号
     * @return 染色体编号对应的所有VCF文件路径
     */
    public static List<String> getAllVcfPathInSpecfiedDir(String vcfDir, int[] chrNumber){
        int temp=(int)Arrays.stream(chrNumber).distinct().count();
        if(chrNumber.length>temp){
            System.out.println("please check your input array, it contains duplicate value");
            System.exit(1);
        }
        List<String> chrNumPathList=new ArrayList<>();
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Pattern p=Pattern.compile("\\d+");
        for(int i=0;i<files.length;i++){
            Matcher m=p.matcher(files[i].getName());
            if(m.find()){
                int num=Integer.valueOf(m.group());
                if(Arrays.stream(chrNumber).anyMatch(e->e==num)){
                    chrNumPathList.add(files[i].getAbsolutePath());
                }
            }
        }
        System.out.println("found "+chrNumPathList.size()+" chromosomes in "+vcfDir+", they are");
        chrNumPathList.stream().forEach(System.out::println);
        return chrNumPathList;
    }

    /**
     * 将输入目录下的所有VCF文件（如"chr001.vcf", "chr002.vcf"等）融合为一个"ChrAll.vcf"文件(不包含"chr000.vcf"、"chr043.vcf"和"chr044.vcf")
     * @param inputVcfDir VCF目录
     */
    public static void mergeNumVcf(String inputVcfDir){
        File[] filesOri=IOUtils.listRecursiveFiles(new File(inputVcfDir));
        File[] files=IOUtils.listFilesEndsWith(filesOri,"vcf");
        Arrays.sort(files);
        VCF vcf1=null;
        Integer flag=0;
        for(int i=0;i<files.length;i++){
            int a=StringTool.getNumFromString(files[i].getName());
            if((a==0)||(a>=43)) continue;
            vcf1=new VCF(files[i]);
            flag=i;
            break;
        }
        for(int i=0;i<files.length;i++){
            int a=StringTool.getNumFromString(files[i].getName());
            if((a==0)||(a>=43)) continue;
            if(i==flag.intValue()) continue;
            vcf1.addVCF(new VCF(files[i]));
        }
        vcf1.write(inputVcfDir,"ChrAll.vcf");
    }

    /**
     * 将输入目录下的所有VCF文件（如"chr001.vcf", "chr002.vcf"等）融合为一个"ChrAll.vcf"文件
     * @param inputVcfDir VCF目录
     * @param contains 是否包含("chr000.vcf"、"chr043.vcf"和"chr044.vcf")
     */
    public static void mergeNumVcf(String inputVcfDir, boolean contains){
        if(contains==true){
            File[] filesOri=IOUtils.listRecursiveFiles(new File(inputVcfDir));
            File[] files=IOUtils.listFilesEndsWith(filesOri,"vcf");
            Arrays.sort(files);
            VCF vcf1=new VCF(files[0]);
            for(int i=1;i<files.length;i++){
                vcf1.addVCF(new VCF(files[i]));
            }
            vcf1.write(inputVcfDir,"ChrAll.vcf");
        }else {
            VCF.mergeNumVcf(inputVcfDir);
        }
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
        Comparator<List<String>> c=Comparator.comparing(e->Integer.valueOf(e.get(0)));
        data.sort(c.thenComparing(e->Integer.valueOf(e.get(1))));
    }

    /**
     *
     * @return 返回所有位点的基因型缺失比例
     */
    public double[] getGenotypeMissingRate(){
        double[] rateOfAllTaxaWithMissingGenotype= new double[data.size()];
        for(int i=0,size=data.size();i<size;i++){
            double numberOfEachTaxaWithMissingGenotype=(double)data.get(i).stream().skip(9).filter(ele->ele.equals("./.")).count();
            rateOfAllTaxaWithMissingGenotype[i]=numberOfEachTaxaWithMissingGenotype/this.getNumberOfTaxa();
        }
        return rateOfAllTaxaWithMissingGenotype;
    }

    /**
     * 将指定的VCF与这个VCF进行融合（两个VCF的meta与header必须相同）
     * @param vcf 指定的VCF对象
     * @return 当VCF对象域data在调用后发生变化时返回true
     */
    public boolean addVCF(VCF vcf){
        return this.data.addAll(vcf.data);
    }

    public void write(String outFile){
        this.sort();
        try(BufferedWriter bw=IOUtils.getTextWriter(outFile)){
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
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void write(String outputDir, String fileName){
       this.write(outputDir+"/"+fileName);
    }

    public void writeVcfToSplitedChrNum(String outputDir){
        this.sort();
        List<Integer> chrNumList=data.stream().flatMap(e->e.stream().limit(1)).mapToInt(Integer::valueOf).boxed()
                              .distinct().collect(Collectors.toCollection(ArrayList::new));
        Collections.sort(chrNumList);
        Map<Integer, BufferedWriter> strToBufferedWriterMap=new HashMap<>();
        Integer key;
        BufferedWriter value;
        for(int i=0;i<chrNumList.size();i++){
            key=chrNumList.get(i);
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
