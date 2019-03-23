package daxing.common;

import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

public class VCF {

    private String meta;
    private List<String> header;
    private List<List<String>> data;

    public VCF(String inputFile){
        this.initilize(inputFile);
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
    public static List<String> getVcfFullPath(String vcfDir, int[] chrNumber){
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


    public static Map<String,Integer> mergeFractionalChr(String input){
        BufferedReader br;
        String temp;
        List<Integer> chr=new ArrayList<>();
        List<Integer> snpNum=new ArrayList<>();
        try {
            br=IOUtils.getTextReader(input);
            while ((temp=br.readLine())!=null){
                List<String> l=PStringUtils.fastSplit(temp);
                chr.add(Integer.valueOf(l.get(0)));

            }
        }catch (Exception e){
            e.printStackTrace();
        }

        return null;
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
}
