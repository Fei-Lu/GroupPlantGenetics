package daxing.common;


import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

/**
 * @author Daxing Xu
 */
public class PlotTools {

    public static void slidingWindow(String inputDir, String outDir, int columnIndex){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        String[] outNames=files.stream().map(File::getName)
                .map(str->str.replaceAll(".vcf", ".sliding.txt")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e-> PlotTools.slidingWindow(files.get(e).getAbsolutePath(),
                columnIndex, 10000000, 1000000, new File(outDir, outNames[e]).getAbsolutePath()));
    }

    /**
     * 将滑窗后生成的所有染色体合并起来
     * @param inputDir
     * @param outFile
     */
    public static void merge(String inputDir, String outFile){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        BufferedReader br;
        BufferedWriter bw=IOTool.getWriter(outFile);
        try {
            String line;
            List<String> temp;
            StringBuilder sb=new StringBuilder();
            bw.write("Chr\tStart\tEnd\tValue");
            bw.newLine();
            for (int i = 0; i < files.size(); i++) {
                br=IOTool.getReader(files.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                    sb.append(temp.get(2)).append("\t").append(temp.get(4));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * CHROM POS WEIR_AND_COCKERHAM_FST
     * 1D	113140	0.250977
     * 1D	119951	0.0114407
     * 1D	128693	-0.00703244
     * 1D	201064	0.287936
     * 1D	554081	0.930793
     * @param rowTable
     * @param columnIndex
     * @param windowSize
     * @param stepSize
     * @param outFile
     * @param sdMedian 是否输出标准差和中位数
     */
    public static void slidingWindow(RowTableTool<String> rowTable, int columnIndex, int windowSize, int stepSize,
                                     String outFile, boolean sdMedian){
        Set<String> chrNum=new HashSet<>(rowTable.getColumn(0));
        if (chrNum.size()!=1){
            System.out.println(chrNum.iterator().next()+" had multiple chromosomes, program quit");
            System.exit(1);
        }
        Comparator<List<String>> c=Comparator.comparing(l->Integer.parseInt(l.get(1)));
        rowTable.sortBy(c);
        List<String> posListStr=rowTable.getColumn(1);
        List<String> valueListStr=rowTable.getColumn(columnIndex);
        TDoubleArrayList keyList=new TDoubleArrayList();
        TDoubleArrayList valueList=new TDoubleArrayList();
        int cnt=0;
        for (int i = 0; i < valueListStr.size(); i++) {
            if (valueListStr.get(i).toUpperCase().equals("NAN")){
                cnt++;
                posListStr.remove(i);
                valueListStr.remove(i);
                i--;
            }
            keyList.add(Double.parseDouble(posListStr.get(i)));
            valueList.add(Double.parseDouble(valueListStr.get(i)));
        }
        System.out.print(chrNum.iterator().next()+" rowTable had "+cnt+" NaN rows");
        String valueName=rowTable.getColumnName(columnIndex).toUpperCase();
        String chrName=rowTable.getCell(0, 0);
        double[] key=keyList.toArray();
        double[] value=valueList.toArray();
        int maxDistace=(int)key[key.length-1];
        int num=(maxDistace-windowSize)/stepSize+2;
        double[] boundaryS= IntStream.iterate(0, n->n+stepSize).limit(num).mapToDouble(n->n-0.1).toArray();
        double[] boundaryL=IntStream.iterate(windowSize, n->n+stepSize).limit(num).mapToDouble(n->n-0.1).toArray();
        int[] indexS=new int[boundaryS.length];
        int[] indexL=new int[boundaryL.length];
        for (int i = 0; i < indexL.length; i++) {
            indexS[i]= Integer.MIN_VALUE;
            indexL[i]=Integer.MIN_VALUE;
        }
        int index=Integer.MAX_VALUE;
        for (int i = 0; i < boundaryS.length; i++) {
            index= Arrays.binarySearch(key, boundaryS[i]);
            indexS[i]=-index-1;
            index=Arrays.binarySearch(key, boundaryL[i]);
            indexL[i]=-index-1;
        }
        int countInBoundary=-1;
        StringBuilder sb;
        int count=0;
        DescriptiveStatistics stats;
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            if (sdMedian){
                bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tMEAN_"+valueName+"\tSD\tMEDIAN\n");
            }
            else {
                bw.write("CHROM\tBIN_START\tBIN_END\tN_VARIANTS\tMEAN_"+valueName+"\n");
            }
            for (int j = 0; j < indexS.length; j++) {
                countInBoundary = indexL[j] - indexS[j];
                if (countInBoundary < 1) {
                    count++;
                }
                sb = new StringBuilder();
                stats = new DescriptiveStatistics();
                for (int i = indexS[j]; i < indexL[j]; i++) {
                    stats.addValue(value[i]);
                }
                sb.append(chrName).append("\t");
                sb.append((int) (boundaryS[j] + 0.1)).append("\t").append((int) (boundaryL[j] + 0.1)).append("\t");
                sb.append(countInBoundary).append("\t");
                if (countInBoundary < 1) {
                    if (sdMedian){
                        sb.append("NaN").append("\t").append("NaN").append("\t").append("NaN").append("\t").append("NaN").append("\n");
                    }
                    else {
                        sb.append("NaN").append("\n");
                    }
                } else {
                    if (sdMedian){
                        sb.append(stats.getMean()).append("\t").append(stats.getStandardDeviation()).append("\t");
                        sb.append(stats.getPercentile(50)).append("\n");
                    }
                    else {
                        sb.append(NumberTool.format(stats.getMean(),5)).append("\n");
                    }
                }
                bw.write(sb.toString());
            }
            if (count == 0) {
                System.out.println(", "+count + " window count were 0");
            } else {
                System.out.println(", "+count + " windows count were 0, and its window MEAN, SD and MEDIAN will be setting NaN");
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void slidingWindow(String inputFile, int columnIndex, int windowSize, int stepSize,
                                     String outFile){
        RowTableTool<String> table=new RowTableTool<>(inputFile);
        slidingWindow(table, columnIndex, windowSize, stepSize, outFile, false);
    }

    public static void windowMean(String ingputFile, int binWidth_kb, int threshForDistance_Mb, String outFile){
        try (BufferedReader br = IOUtils.getTextReader(ingputFile);
             BufferedWriter bw=IOUtils.getTextWriter(outFile)) {
            String line;
            List<String> temp;
            int distance=Integer.MIN_VALUE;
            double r2=Double.MIN_VALUE;
            int thresh=binWidth_kb*1000;
            int limit=threshForDistance_Mb*1000000;
            int kb=binWidth_kb;
            DescriptiveStatistics r2Stats=new DescriptiveStatistics();
            br.readLine();
            bw.write("window_kb\tnumberInWindow\tmeanOfR2\n");
            StringBuilder sb;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                distance=Integer.parseInt(temp.get(0));
                if (distance > limit) break;
                r2=Double.parseDouble(temp.get(1));
                if (distance > thresh){
                    sb=new StringBuilder();
                    sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
                    bw.write(sb.toString());
                    bw.newLine();
                    thresh+=binWidth_kb*1000;
                    kb+=binWidth_kb;
                    r2Stats=new DescriptiveStatistics();
                }
                r2Stats.addValue(r2);
            }
            sb=new StringBuilder();
            sb.append(kb).append("\t").append(r2Stats.getN()).append("\t").append(r2Stats.getMean());
            bw.write(sb.toString());
            bw.newLine();
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public static void countHighLoadFrequency(String inputFile, int index, String outFile, int binNum, double binWidth){
        if (binWidth*binNum!=1){
            System.out.println("binNum*binWidth must be 1, program quit! ");
            System.exit(1);
        }
        double[] binStart= new double[binNum];
        double[] binEnd  = new double[binNum];
        int[] binCount=new int[binNum];
        double[] binRatio=null;
        for (int i = 0; i < binNum; i++) {
            binStart[i]=binWidth*i;
            binEnd[i]=binWidth*i+binWidth;
        }
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            br.readLine();
            bw.write("BinID\tBinStart\tBinEnd\tBinNum\tBinRatio");
            bw.newLine();
            String line;
            List<String> temp;
            int binStartIndex=-1;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                if (temp.get(index).equals("NA")) continue;
                int hit=Arrays.binarySearch(binStart, Double.parseDouble(temp.get(index)));
                binStartIndex= hit < 0 ? -hit-2 : hit;
                binCount[binStartIndex]++;
            }
            binRatio=ArrayTool.getElementPercent(binCount);
            StringBuilder sb;
            String binID;
            for (int i = 0; i < binCount.length; i++) {
                sb=new StringBuilder();
                binID="B"+PStringUtils.getNDigitNumber(5, i);
                sb.append(binID).append("\t").append(binStart[i]).append("\t").append(binEnd[i]).append("\t");
                sb.append(binCount[i]).append("\t").append(binRatio[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
