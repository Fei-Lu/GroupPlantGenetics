package xiaohan.utils;

/**
 * Utils be used in search and sort
 *
 * @ author: yxh
 * @ created: 2021-11-30 : 9:47 AM
 */
public class Alogrithm {

    public Alogrithm(String[] args) {

    }

    public static int BinarySearch(Long[] arr,Long key,int low,int high){
        int a = -1;
        if(key < arr[low] || key > arr[high] || low > high){
            return -1;
        }

        int middle = (low + high) / 2;			//初始中间位置
        if(arr[middle] > key){
            //比关键字大则关键字在左区域
            return BinarySearch(arr, key, low, middle - 1);
        }else if(arr[middle] < key){
            //比关键字小则关键字在右区域
            return BinarySearch(arr, key, middle + 1, high);
        }else {
            return middle;
        }
    }

    public static void main(String[] args) {
        new Alogrithm(args);
    }
}