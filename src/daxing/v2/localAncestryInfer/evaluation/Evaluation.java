package daxing.v2.localAncestryInfer.evaluation;

import daxing.common.sh.CommandUtils;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class Evaluation {

    /**
     * detail information about Mean deviation and Pearson correlation see reference
     * Title: Detecting Structure of Haplotypes and Local Ancestry
     * Journal: Genetics
     * Year: 2014
     * Author: Yongtao Guan
     * DOI: 10.1534/genetics.113.160697
     */
    public enum Metric{
        MEAN_DEVIATION,
        PEARSON_CORRELATION
    }


    /**
     *
     * @param inferredValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @param actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @return Mean deviation or Pearson correlation, dim1 is different run, dim2 is different haplotype
     */
    public static double[][] evaluate(double[][][][] inferredValue, double[][][][] actualValue, Metric metric){
        double[][] result = new double[inferredValue.length][];
        switch (metric){
            case MEAN_DEVIATION:
                result = meanDeviation(inferredValue, actualValue);
                break;
            case PEARSON_CORRELATION:
                result = pearsonCorrelation(inferredValue, actualValue);
                break;
        }
        return result;

    }

    /**
     *
     * @param inferredValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @param actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @return Mean deviation, dim1 is different run, dim2 is different haplotype
     */
    public static double[][] meanDeviation(double[][][][] inferredValue, double[][][][] actualValue){
        double[][] totalDeviation = new double[inferredValue.length][inferredValue[0].length];
        for (int i = 0; i < inferredValue.length; i++) {
            for (int j = 0; j < inferredValue[i].length; j++) {
                for (int k = 0; k < inferredValue[i][j].length; k++) {
                    for (int l = 0; l < inferredValue[i][j][k].length; l++) {
                        totalDeviation[i][j]+=Math.abs(actualValue[i][j][k][l] - inferredValue[i][j][k][l]);
                    }
                }
            }
        }
        int[][] count = new int[inferredValue.length][inferredValue[0].length];
        for (int i = 0; i < inferredValue.length; i++) {
            for (int j = 0; j < inferredValue[i].length; j++) {
                count[i][j] = inferredValue[i][j].length * inferredValue[i][j][0].length;
            }
        }
        double[][] meanDeviation = new double[inferredValue.length][inferredValue[0].length];
        for (int i = 0; i < totalDeviation.length; i++) {
            for (int j = 0; j < totalDeviation[i].length; j++) {
                meanDeviation[i][j] = totalDeviation[i][j]/count[i][j];
            }
        }
        return meanDeviation;
    }

    /**
     *
     * @param inferredValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @param actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @return Pearson correlation, dim1 is different run, dim2 is different haplotype
     */
    public static double[][] pearsonCorrelation(double[][][][] inferredValue, double[][][][] actualValue){
        List<Callable<Double>> callableList= new ArrayList<>();
        for (int i = 0; i < inferredValue.length; i++) {
            for (int j = 0; j < inferredValue[i].length; j++) {
                for (int k = 0; k < inferredValue[i][j].length; k++) {
                    int finalI = i;
                    int finalJ = j;
                    int finalK = k;
                    callableList.add(() -> new PearsonsCorrelation().correlation(inferredValue[finalI][finalJ][finalK],
                            actualValue[finalI][finalJ][finalK]));
                }
            }
        }
        List<Double> res = CommandUtils.run_commands(callableList, 32);
        double[][] result = new double[inferredValue.length][inferredValue[0].length];
        for (int i = 0; i < res.size(); i++) {
            result[i/inferredValue[0].length][i%inferredValue[0].length] = res.get(i);
        }
        return result;
    }
}
