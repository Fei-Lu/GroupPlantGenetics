/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.util.Random;

/**
 *
 * @author kanglipeng
 */
public class RandomArray {
    public static int[] randomArray(int min,int max,int n){  
    int len = max-min+1;  
      
    if(max < min || n > len){  
        return null;  
    }  
      
    //初始化给定范围的待选数组  
    int[] source = new int[len];  
       for (int i = min; i < min+len; i++){  
        source[i-min] = i;  
       }  
         
       int[] result = new int[n];  
       Random rd = new Random();  
       int index = 0;  
       for (int i = 0; i < result.length; i++) {  
        //待选数组0到(len-2)随机一个下标  
           index = Math.abs(rd.nextInt() % len--);  
           //将随机到的数放入结果集  
           result[i] = source[index];  
           //将待选数组中被随机到的数，用待选数组(len-1)下标对应的数替换  
           source[index] = source[len];  
       }  
       return result;  
}  

  public static int[] randomCommon(int min, int max, int n){  
    if (n > (max - min + 1) || max < min) {  
           return null;  
       }  
    int[] result = new int[n];  
    int count = 0;  
    while(count < n) {  
        int num = (int) (Math.random() * (max - min)) + min;  
        boolean flag = true;  
        for (int j = 0; j < n; j++) {  
            if(num == result[j]){  
                flag = false;  
                break;  
            }  
        }  
        if(flag){  
            result[count] = num;  
            count++;  
        }  
    }  
    return result;  
}    
  
    public static long[] randomLongCommon(long min, long max, int n){  
    if (n > (max - min + 1) || max < min) {  
           return null;  
       }  
    long[] result = new long[n];  
    int count = 0;  
    while(count < n) {  
        long num = (long) (Math.random() * (max - min)) + min;  
        boolean flag = true;  
        for (int j = 0; j < n; j++) {  
            if(num == result[j]){  
                flag = false;  
                break;  
            }  
        }  
        if(flag){  
            result[count] = num;  
            count++;  
        }  
    }  
    return result;  
}  
}

