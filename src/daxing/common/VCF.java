package daxing.common;

import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Random;
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
     * @param numberOfRowsForExtract 提取的行数
     */
    public static void extractRandomRowFromVCF(String inputFile, String outFile, int numberOfRowsForExtract) {
        long start = System.nanoTime();
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            List<String> metaAndHeader= Files.newBufferedReader(Paths.get(inputFile))
                    .lines()
                    .filter(index->index.startsWith("#"))
                    .collect(Collectors.toList());
            List<String> data=Files.newBufferedReader(Paths.get(inputFile))
                    .lines()
                    .filter(index -> (!index.startsWith("#")))
                    .collect(Collectors.toList());
            Collections.shuffle(data, new Random());
            for(String str:metaAndHeader){
                bw.write(str);
                bw.newLine();
            }
            for(int i=0;i<numberOfRowsForExtract;i++){
                bw.write(data.get(i));
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
        System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanSeconds(start)) + " seconds");
    }

    /**
     * 返回VCF文件所有taxa的总数
     * @return 返回VCF文件所有taxa的总数
     */
    public int getNumberOfTaxa(){
        return header.size()-9;
    }

    public double[] getGenotypeMissingRate(){
        double[] rateOfAllTaxaWithMissingGenotype= new double[data.size()];
        for(int i=0,size=data.size();i<size;i++){
            double numberOfEachTaxaWithMissingGenotype=(double)data.get(i).stream().skip(9).filter(ele->ele.equals("./.")).count();
            rateOfAllTaxaWithMissingGenotype[i]=numberOfEachTaxaWithMissingGenotype/this.getNumberOfTaxa();
        }
        return rateOfAllTaxaWithMissingGenotype;
    }
}
