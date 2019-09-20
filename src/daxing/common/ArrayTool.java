package daxing.common;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ArrayTool {

    /**
     * 返回指定数值范围的随机数组（包含重复值）
     * @param arraySize 返回的数组大小
     * @param randomNumberOrigin  inclusive
     * @param randomNumberBound exclusive
     * @return 依据给定的数组范围返回一个随机数组（包含重复值), 如果给定的数值范围小于指定的数组大小, 则返回的数组大小为数值范围
     */
    public static int[] getRandomNumberArray(int arraySize, int randomNumberOrigin, int randomNumberBound){
        return new java.util.Random().ints(arraySize,randomNumberOrigin,randomNumberBound).toArray();
    }

    /**
     * 返回一个指定数值范围的随机数组（不包含重复值）
     * @param arraySize 返回的数组大小
     * @param randomNumberOrigin inclusive
     * @param randomNumberBound exclusive
     * @return 依据给定的数组范围返回一个随机数组（不含重复值), 如果给定的数值范围小于指定的数组大小, 则返回的数组大小为数值范围
     */
    public static int[] getRandomNonrepetitionArray(int arraySize, int randomNumberOrigin, int randomNumberBound){
        List<Integer> list= IntStream.range(randomNumberOrigin,randomNumberBound).boxed().collect(Collectors.toList());
        Collections.shuffle(list);
        return list.stream().mapToInt(Integer::intValue).limit(arraySize).toArray();
    }

    /**
     *
     * @param a array
     * @return the number of occurrences of each element in the array
     */
    public static Map<Integer, Long> calcalateRepetitionNum(int[] a){
        List<Integer> l= Arrays.stream(a).boxed().collect(Collectors.toList());
        return l.parallelStream().collect(Collectors.groupingByConcurrent(Function.identity(), Collectors.counting()));
    }

    /**
     *  将两个数组对应的index元素相加
     * @param a
     * @param b
     * @return
     */
    public static int[] add(int[] a, int[] b){
        if (a.length!=b.length){
            System.out.println(a+" and "+b+" length is not same");
            System.exit(1);
        }
        int[] c=new int[a.length];
        for (int i = 0; i < a.length; i++) {
            c[i]=a[i]+b[i];
        }
        return c;
    }

    /**
     * 将两个数组对应的index元素相加
     * @param a
     * @param b
     * @return
     */
    public static double[] add(double[] a, double[] b){
        if (a.length!=b.length){
            System.out.println(a+" and "+b+" length is not same");
            System.exit(1);
        }
        double[] c=new double[a.length];
        for (int i = 0; i < a.length; i++) {
            c[i]=a[i]+b[i];
        }
        return c;
    }

    /**
     * 返回数组每个元素的比例
     * @param a
     * @return
     */
    public static double[] getElementPercent(int[] a){
        double sum=Arrays.stream(a).sum();
        return Arrays.stream(a).mapToDouble(e->e/sum).toArray();
    }
}
