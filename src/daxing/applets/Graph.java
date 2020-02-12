package daxing.applets;

import daxing.common.ArrayTool;
import daxing.common.NumberTool;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.Arrays;
import java.util.List;
import java.util.stream.DoubleStream;

public class Graph {

    /**
     * first column for graph
     * @param inputFile
     * @param binwidth
     * @param binsNum
     * @param columnIndex
     * @param outFile
     */
    public static void histogram(String inputFile, double binwidth, int binsNum, int columnIndex, String outFile){
        if (binwidth*binsNum!=1){
            System.out.println("please check your parameters");
            System.exit(1);
        }
        double[] bins= DoubleStream.iterate(binwidth, n->n+binwidth).map(d-> NumberTool.format(d, 2)).limit(binsNum).toArray();
        int[] count=new int[bins.length];
        for (int i = 0; i < count.length; i++) {
            count[i]=0;
        }
        int aa=0;
        String line=null;
        try (BufferedReader br = IOUtils.getTextReader(inputFile);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            List<String> lineList;
            br.readLine();
            bw.write("boundary"+"\t"+"count"+"\t"+"rate");
            bw.newLine();
            int index=Integer.MIN_VALUE;
            double maf=-1;
            while ((line=br.readLine())!=null){
                aa++;
                lineList= PStringUtils.fastSplit(line);
//                if (lineList.get(0).equals("NaN")) continue;
                maf=Double.parseDouble(lineList.get(columnIndex));
                index= Arrays.binarySearch(bins, maf);
                if (index<0){
                    index=-index-1;
                    count[index]++;
                }else {
                    count[index]++;
                }

            }
            double[] rate= ArrayTool.getElementPercent(count);
            StringBuilder sb;
            for (int i = 0; i < count.length; i++) {
                sb=new StringBuilder();
                sb.append(bins[i]).append("\t").append(count[i]).append("\t").append(rate[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
        }catch (Exception e){
            System.out.println(line);
            System.out.println(aa);
            e.printStackTrace();
        }
    }
}
