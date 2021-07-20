package xiaohan.rareallele.getmedZ;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @ author: yxh
 * @ created: 2021-07-20 : 10:38 PM
 */
public class medZUtils {


    public static String medn(String[] array) {
        int len = 0;
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                len++;
            }
        }
        String length = String.valueOf(len);
        return length;
    }

    public static String medmedian(String[] array) {
        HashSet<String> DataSet = new HashSet<>();
        String median = null;
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                DataSet.add(array[i]);
            }
        }
        if (DataSet.isEmpty()) {
            median = "NA";
        } else {
            String[] data = DataSet.toArray(new String[0]);
            double[] datas = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                datas[i] = Double.parseDouble(data[i]);
            }
            Arrays.sort(datas);
            int index = data.length / 2;
            if (data.length % 2 == 1) {
                median = String.valueOf(datas[datas.length / 2]);
            } else {
                median = String.valueOf((datas[index] + datas[index - 1]) / 2);
            }
        }
        return median;
    }

    public static String[] pickoutlier(String[] array) {
        HashMap<String, Integer> IndexValues = new HashMap<>();
        for (int i = 0; i < array.length; i++) {
            IndexValues.put(array[i], i);
        }
        String max = null;
        String rank = null;
        HashSet<Double> DataSet = new HashSet<>();
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                DataSet.add(Double.parseDouble(array[i]));
            }
        }
        if (DataSet.isEmpty()) {
            max = "NA";
            rank = "NA";
        } else {
            Object[] data = DataSet.toArray();
            double[] abs = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                abs[i] = Math.abs((double) data[i]);
            }
            int indexmax = (int) getMaxIndex(abs)[1];
            max = String.valueOf((double) data[indexmax]);
            rank = String.valueOf(IndexValues.get(max));
        }
        String[] results = new String[]{max, rank};
        return results;
    }


    public static double[] getMaxIndex(double[] arr) {
        if (arr == null || arr.length == 0) {
            return null;//如果数组为空 或者是长度为0 就返回null
        }
        int maxIndex = 0;//假设第一个元素为最大值 那么下标设为0
        double[] arrnew = new double[2];//设置一个长度为2的数组用作记录，第一个元素存储最大值，第二个元素存储下标
        for (int i = 0; i < arr.length - 1; i++) {
            if (arr[maxIndex] < arr[i + 1]) {
                maxIndex = i + 1;
            }
        }
        arrnew[0] = arr[maxIndex];
        arrnew[1] = maxIndex;
        return arrnew;
    }

}
