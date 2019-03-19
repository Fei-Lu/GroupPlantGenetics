package daxing.common;

import utils.Benchmark;
import utils.IOUtils;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;

public class VCF {

    /**
     * 从VCF文件中随机抽取行，组成新文件
     * @param inputFile 输入文件
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
        System.out.println("smallVCF is completed in " + String.format("%.4f", Benchmark.getTimeSpanSeconds(start)) + " seconds");
    }
}
