package xiaohan.utils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;

public class SNPmappingInGene {
    //array in two values for start/end or three values for chr/start/end
    public static int[] bSearch(int[][] a, int key) {
        if (a == null || a.length <= 0 || key > a[a.length - 1][a[a.length - 1].length - 1]) {
            return new int[]{-1};
        }
        int low = 0;
        int high = a.length - 1;
        int mid = 0;
        ArrayList<Integer> res = new ArrayList<>();
        while (low <= high) {
            mid = low + (high - low) / 2;
            if (a[mid][0] > key) {
                high = mid - 1;
            } else if (a[mid][1] < key) {
                low = mid + 1;
            } else {
                break;
            }
        }

        if (a[mid][0] > key) {
            res.add(mid - 1);
        }
        res.add(mid);
        if (a[mid][1] < key) {
            res.add(mid + 1);
        }

        int[] ints = new int[res.size()];
        for (int i = 0; i < res.size(); i++) {
            ints[i] = res.get(i);
        }
        return ints;
//        return new int[]{-1};
    }

    public static int[] binarySearch(int[][] array, int value) {
        if (array == null || array.length <= 0 || value > array[array.length - 1][array[array.length - 1].length - 1]) {
            return new int[]{-1};
        }
        int colums = array[array.length - 1].length;
        int low = 0;
        int high = array.length - 1;
        int mid = 0;

        ArrayList<Integer> res = new ArrayList<>();

        while (low <= high) {
            mid = low + (high - low) / 2;
            if (array[mid][colums - 2] <= value && array[mid][colums - 1] >= value) {
                res.add(mid);
                low = mid + 1;
            } else if (array[mid][colums - 1] < value) {
                low = mid + 1;
            } else if (array[mid][colums - 2] > value) {
                high = mid - 1;
            } else break;
        }

        if(res.isEmpty()){
            res.add(-1);
        }

        int[] ints = new int[res.size()];
        for (int i = 0; i < res.size(); i++) {
            ints[i] = res.get(i);
        }
        return ints;
    }

    public static double[] Min(double[] arr){
        double min1 = arr[0];
        int index1 = 0;
        for (int i = 1; i < arr.length; i++) {
            if(arr[i] < min1){
                min1 = arr[i];
                index1 = i;
            }
        }
        double[] minMin = {min1,index1};
        return minMin;
    }

    public static double[] minMin(double[] arr){
        double min1 = arr[0], min2 = arr[1];
        int index1 = 0;
        int index2 = 1;
        if(min1 > min2){
            double temp = min1;
            min1 = min2;
            min2 = temp;
            index1 = 1;
            index2 = 0;
        }
        for (int i = 2; i < arr.length; i++) {
            if(arr[i] < min1){
                double temp = min1;
                min1 = arr[i];
                min2 = temp;
                index2 = index1;
                index1 = i;
            }else if(arr[i] < min2){
                min2 = arr[i];
                index2 = i;
            }
        }
        double[] minMin = {min1,min2,index1,index2};
        return minMin;
    }

    public static double[] minMinnumber(double[] arr,int number){
        double min1 = arr[0], min2 = arr[1];
        int index1 = 0;
        int index2 = 1;
        if(min1 > min2){
            double temp = min1;
            min1 = min2;
            min2 = temp;
            index1 = 1;
            index2 = 0;
        }
        for (int i = 2; i < arr.length; i++) {
            if(arr[i] < min1){
                double temp = min1;
                min1 = arr[i];
                min2 = temp;
                index2 = index1;
                index1 = i;
            }else if(arr[i] < min2){
                min2 = arr[i];
                index2 = i;
            }
        }
        double[] minMin = {min1,min2,index1,index2};
        return minMin;
    }
}
