package xiaohan.utils;

/**
 * @ author: yxh
 * @ created: 2021-07-18 : 8:28 PM
 */
public class MathUtils {

    public static double getsd(double[] array){
        double sd = 0;
        double exp = 0;
        double average = 0;
        double res = 0;
        double zscore = 0;
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            average += exp;
        }
        average = average/array.length;
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            res = (exp - average)*(exp - average) ;
            sd += res;
        }
        sd = sd/(array.length-1);
        sd = Math.sqrt(sd);
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            zscore = (exp - average)/sd;
        }
        return sd;
    }

    public static double[] getzscore(double[] array){
        double sd = 0;
        double exp = 0;
        double average = 0;
        double res = 0;
        double[] zscore = new double[array.length];
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            average += exp;
        }
        average = average/array.length;
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            res = (exp - average)*(exp - average) ;
            sd += res;
        }
        sd = sd/(array.length-1);
        sd = Math.sqrt(sd);
        for (int i = 0; i < array.length; i++) {
            exp = array[i];
            zscore[i] = (exp - average)/sd;
        }
        return zscore;
    }
}
