package daxing.common;

import utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class MyRandom {

    /**
     * 返回指定数值范围的随机数组（包含重复值）
     * @param arraySize 返回的数组大小
     * @param randomNumberOrigin 数组范围的最小值（包含）
     * @param randomNumberBound 数组范围的最大值（不包含）
     * @return 依据给定的数组范围返回一个随机数组（包含重复值), 如果给定的数值范围小于指定的数组大小, 则返回的数组大小为数值范围
     */
    public static int[] getRandomNumberArray(int arraySize, int randomNumberOrigin, int randomNumberBound){
        return new java.util.Random().ints(arraySize,randomNumberOrigin,randomNumberBound).toArray();
    }

    /**
     * 返回一个指定数值范围的随机数组（不包含重复值）
     * @param arraySize 返回的数组大小
     * @param randomNumberOrigin 数组范围的最小值（包含）
     * @param randomNumberBound 数组范围的最大值（不包含）
     * @return 依据给定的数组范围返回一个随机数组（不含重复值), 如果给定的数值范围小于指定的数组大小, 则返回的数组大小为数值范围
     */
    public static int[] getRandomNonrepetitionArray(int arraySize, int randomNumberOrigin, int randomNumberBound){
        List<Integer> list= IntStream.range(randomNumberOrigin,randomNumberBound).boxed().collect(Collectors.toList());
        Collections.shuffle(list);
        return list.stream().mapToInt(Integer::intValue).limit(arraySize).toArray();
    }

    /**
     * 从一个文件中随机抽取行，组成新文件
     * @param inputFile 输入文件
     * @param outFile 输出文件
     * @param numberOfRowsForExtract 提取的行数
     * @param head 输入文件是否包含表头
     */
    public static void extractRandomRowFromFile(String inputFile, String outFile, int numberOfRowsForExtract,
                                                boolean head){
        try(BufferedReader br= IOUtils.getTextReader(inputFile); BufferedWriter bw =IOUtils.getTextWriter(outFile)){
            long numberOfRows = Files.lines(Paths.get(inputFile)).count();
            List<String> l= new ArrayList<>(numberOfRowsForExtract);
            String line;
            String header=null;
            while ((line=br.readLine())!=null){
                if (head) {header=br.readLine();}
                line=br.readLine();
                l.add(line);
            }
            Collections.shuffle(l, new Random());
            if(head) {bw.write(header);}
            for(int i=0;i<numberOfRowsForExtract;i++){
                bw.write(l.get(i));
                bw.newLine();
            }
            bw.flush();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }

}
