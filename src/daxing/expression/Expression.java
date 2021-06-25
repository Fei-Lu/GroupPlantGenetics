package daxing.expression;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;

import java.io.*;
import java.util.List;

public class Expression {

    /**
     * Pearsons Correlation Matrix
     * @param inputFile  Azhurnaya Development_tpm.tsv.gz
     * @param outFile matrix
     */
    public static void getPearsonsCorrelationMatrix(String inputFile, String outFile){
        long start=System.nanoTime();
        RowTableTool<String> table = new RowTableTool<>(inputFile);
        double[][] value=getMatrix(table);
        PearsonsCorrelation pearsonsCorrelation = new PearsonsCorrelation();
        long start1=System.nanoTime();
        RealMatrix realMatrix = pearsonsCorrelation.computeCorrelationMatrix(value);
        System.out.println("Pearsons correlation matrix Complicated in "+ Benchmark.getTimeSpanHours(start1)+ " " +
                "hours");
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb = new StringBuilder();
            List<String> header= table.getHeader();
            sb.append("GeneName").append("\t").append(String.join("\t", header));
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < realMatrix.getRowDimension(); i++) {
                sb.setLength(0);
                sb.append(header.get(i)).append("\t");
                for (int j = 0; j < realMatrix.getColumnDimension(); j++) {
                    sb.append(realMatrix.getEntry(i,j)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                if (i % 1000 ==0){
                    System.out.println(i + " genes had been written to "+ new File(outFile).getName());
                }
            }
            bw.flush();
            System.out.println("Total "+ header.size() + " genes had been written to " + new File(outFile).getName());
            System.out.println("Complicated in "+ Benchmark.getTimeSpanHours(start)+ " hours");
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private static double[][] getMatrix(RowTableTool<String> table){
        List<List<String>> matrixList=table.getCells();
        double[][] matrix = new double[matrixList.size()][matrixList.get(0).size()];
        double[] line;
        for (int i = 0; i < matrix.length; i++) {
            line = new double[matrix[i].length];
            for (int j = 0; j < matrix[i].length; j++) {
                line[j]= Double.parseDouble(matrixList.get(i).get(j));
            }
            matrix[i]=line;
        }
        return matrix;
    }

    /**
     * subset from Pearsons correlation matrix
     * @param inputFile
     * @param rate
     * @param outFile
     */
    public static void subset(String inputFile, double rate, String outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, headerList;
            headerList = PStringUtils.fastSplit(br.readLine());
            int threshold=1;
            double random;
            StringBuilder sb = new StringBuilder();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                for (int i = 1; i < temp.size(); i++) {
                    sb.setLength(0);
                    if (i > threshold) break;
                    random=Math.random();
                    if (random > rate) continue;
//                    sb.append(headerList.get(i)).append("\t").append(temp.get(0)).append("\t");
                    sb.append(temp.get(i));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                threshold++;
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * the connection degree of a gene (The number of genes exceeding the threshold)
     * @param pearsonsCorrelationFile
     * @param threshold
     * @param outFile
     */
    public static void calculateConnectionDegree(String pearsonsCorrelationFile, double threshold, String outFile){
        long start= System.nanoTime();
        try (BufferedReader br = IOTool.getReader(pearsonsCorrelationFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            double r;
            int cntPositive=0, cnt1Negative=0;
            StringBuilder sb = new StringBuilder();
            br.readLine();
            bw.write("Gene\tPositive\tNegative");
            bw.newLine();
            int cnt=0;
            while ((line=br.readLine())!=null){
                cnt++;
                temp=PStringUtils.fastSplit(line);
                sb.setLength(0);
                for (int i = 1; i < temp.size(); i++) {
                    r = Double.parseDouble(temp.get(i));
                    if (Double.isNaN(r)) continue;
                    if (Math.abs(r) < threshold) continue;
                    if (r < 0){
                        cnt1Negative++;
                    }else {
                        cntPositive++;
                    }
                }
                sb.append(temp.get(0)).append("\t");
                sb.append(cntPositive-1).append("\t").append(cnt1Negative);
                bw.write(sb.toString());
                bw.newLine();
                cntPositive=0;
                cnt1Negative=0;
                if (cnt % 1000 ==0){
                    System.out.println("Writing "+ cnt +" genes to "+ new File(outFile).getName());
                }
            }
            bw.flush();
            System.out.println("total "+ cnt +" genes had been written to "+ new File(outFile).getName());
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Connection degree based on threshold "+ threshold + " complicated in "+ Benchmark.getTimeSpanHours(start)+ " hours");
    }

    public static void calculateConnectionDegreeSum(String pearsonsCorrelationFile, String outFile){
        long start=System.nanoTime();
        try (BufferedReader br = IOTool.getReader(pearsonsCorrelationFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, headerList;
            double r;
            double totalNegativeR=0, totalPositiveR=0;
            int countNegativeR=0, countPositiveR=0;
            StringBuilder sb = new StringBuilder();
            headerList=PStringUtils.fastSplit(br.readLine());
            bw.write("Gene\tPositiveRSum\tPositiveRCount\tNegativeRSum\tNegativeRCount");
            bw.newLine();
            int cnt=0;
            while ((line=br.readLine())!=null){
                cnt++;
                temp=PStringUtils.fastSplit(line);
                sb.setLength(0);
                for (int i = 1; i < temp.size(); i++) {
                    r = Double.parseDouble(temp.get(i));
                    if (Double.isNaN(r)) continue;
                    if (r < 0){
                        totalNegativeR+=r;
                        countNegativeR++;
                    }else {
                        totalPositiveR+=r;
                        countPositiveR++;
                    }
                }
                sb.append(temp.get(0)).append("\t");
                sb.append(totalPositiveR-1).append("\t").append(countPositiveR-1).append("\t");
                sb.append(totalNegativeR).append("\t").append(countNegativeR);
                bw.write(sb.toString());
                bw.newLine();
                totalNegativeR=0;
                totalPositiveR=0;
                countNegativeR=0;
                countPositiveR=0;
                if (cnt % 1000 ==0){
                    System.out.println("Writing "+ cnt +" genes to "+ new File(outFile).getName());
                }
            }
            bw.flush();
            System.out.println("total "+ cnt +" genes had been written to "+ new File(outFile).getName());
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Connection degree based on sum" + " complicated in "+ Benchmark.getTimeSpanHours(start)+ " hours");
    }
}
