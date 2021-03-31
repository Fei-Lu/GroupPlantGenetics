package xiaohan.eQTL;

import com.mitchellbosecke.pebble.extension.Test;
import xiaohan.eQTL.SNPmapping;

import java.util.Arrays;


public class javaPractice {
    public javaPractice(String...args){
//        this.attributionForCollections();
//        this.LeetCodePractice();
        this.testSNPmaping();
    }

    public void testSNPmaping(){
        int[][] array1 = {{1,3,4},{1,6,12},{1,6,14},{1,14,19}};
//        for (int i = 0; i < array1.length; i++) {
//            for (int j = 0; j < array1[0].length; j++) {
//                System.out.println(array1[i][j]);
//            }
//        }
//        int[] index = SNPmapping.bSearch(array1,11);
        int[] index1 = SNPmappingInGene.binarySearch(array1,5);
        if(index1[0]!=-1) {
            for (int i = 0; i < index1.length; i++) {
                System.out.println(index1[i]);
            }
        }

//        array1 = Arrays.copyOfRange(array1,2,3);
//        for (int i = 0; i < array1.length; i++) {
//            for (int j = 0; j < array1[0].length; j++) {
//                System.out.println(array1[i][j]);
//            }
//        }
//        int a = 3/2;
//        System.out.println(a);
    }

    public void LeetCodePractice(){
        new LeetCodePractice();
    }

    public void attributionForCollections(){
        new attributionForCollections();
    }

    public static void main (String[] args){
        new javaPractice(args);
    }
}
