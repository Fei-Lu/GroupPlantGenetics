package xiaohan.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @ author: yxh
 * @ created: 2021-07-18 : 8:28 PM
 */
public class MathUtils {

    /**
     * @param array
     * get average or mean from an array
     * @return
     */
    public static double getmean(double[] array) {
        double average = 0;
        for (int i = 0; i < array.length; i++) {
            average += array[i];
        }
        average = average / array.length;
        return average;
    }

    /**
     * @param array
     * get standard variation or sd from an array
     * @return
     */
    public static double getsd(double[] array) {
        double sd = 0;
        double average = getmean(array);
        for (int i = 0; i < array.length; i++) {
            sd += (array[i] - average) * (array[i] - average);
        }
        sd = sd / (array.length - 1);
        sd = Math.sqrt(sd);
        return sd;
    }

    /**
     * @param array
     * get variance from an array
     * @return
     */
    public static double getvar(double[] array) {
        double var = 0;
        double average = getmean(array);
        for (int i = 0; i < array.length; i++) {
            var += (array[i] - average) * (array[i] - average);
        }
        var = var / (array.length - 1);
        return var;
    }

    /**
     * @param array
     * get zscores from an array
     * @return
     */
    public static double[] getzscore(double[] array) {
        double sd = getsd(array);
        double average = getmean(array);
        double[] zscore = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            zscore[i] = (array[i] - average) / sd;
        }
        return zscore;
    }

    /**
     * @param array
     * get median form an array
     * @return
     */

    public static double getmedian(double[] array) {
        double median = 0;
        Arrays.sort(array);
        if (array.length % 2 == 1) {
            median = array[array.length / 2];
        } else {
            median = (array[array.length / 2] + array[array.length / 2 - 1]) / 2;
        }
        return median;
    }

    /**
     * @param array
     * get range from an array
     * @return
     */
    public static double getrange(double[] array){
        double range = 0;
        Arrays.sort(array);
        range = array[array.length-1]-array[0];
        return range;
    }

    /**
     * @param array
     * @return
     */
    public static double getmax(double[] array){
        double max = 0;
        Arrays.sort(array);
        max = array[array.length-1]-array[0];
        return max;
    }

    /**
     * @param array
     * @return
     */
    public static double getmin(double[] array){
        double min = 0;
        Arrays.sort(array);
        min = array[0];
        return min;
    }

    /**
     * @param array
     * @return
     */
    public static double getmode(double[] array){
        double mode = 0;
        Arrays.sort(array);
        HashMap<Double,Integer> num = new HashMap<>();
        for (int i = 0; i < array.length;) {
            int count = 1;
            for (int j = i + 1; j < array.length; j++) {
                if (array[i] == array[j]) {
                    count++;
                } else {
                    break;
                }
            }
            num.put(array[i],count);
        }
        return mode;
    }

    public static double PearsonCor(double[] array1,double[] array2){
        double pearson = 0;
        if(array1.length!=array2.length){
            sortGeneByStartPosition();
        }else {
            double x = 0;
            double y = 0;
            double xy = 0;
            double x2 = 0;
            double y2 = 0;
            for (int i = 0; i < array1.length; i++) {
                x += array1[i];
                y += array2[i];
                xy += array1[i] * array2[i];
                x2 += array1[i] * array1[i];
                y2 += array2[i] * array2[i];
            }
            pearson = (xy - x*y) / (Math.sqrt(x2-x*x) * Math.sqrt(y2-y*y));
        }
        return pearson;
    }

    public static double SpearmanCor(double[] array1,double[] array2){
        double spearman = 0;
        return spearman;
    }

    public static double KendallCor(double[] array1,double[] array2){
        double kendall = 0;
        return kendall;
    }


    public static void sortGeneByStartPosition () {
        Boolean condition = false;
        System.out.println("X and Y have different lengths");
    }
}
