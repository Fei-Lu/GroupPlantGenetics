package daxing.common;

import java.util.Arrays;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * @author Daxing Xu
 */
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
        return new java.util.Random().ints(randomNumberOrigin, randomNumberBound).distinct().limit(arraySize).toArray();
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
     * the number of occurrences of each element in the array
     * @param a
     * @param <T>
     * @return
     */
    public static<T> Map<T, Long> caculateElementCount(T[] a){
        return Arrays.stream(a).collect(Collectors.groupingBy(t->t, Collectors.counting()));
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
        double sum=Arrays.stream(a).mapToDouble(Double::valueOf).sum();
        return Arrays.stream(a).mapToDouble(e->e/sum).toArray();
    }

    /**
     *
     * @param a
     * @return 数组最大元素的索引, 如果有多个最大元素, 则返回多个最大元素对应的最小索引
     */
    public static int getMaxIndex(int[] a){
        return IntStream.range(0, a.length).boxed().max(Comparator.comparing(e->a[e])).get();
    }

    /**
     *
     * @param a
     * @return 数组最大元素的索引, 如果有多个最大元素, 则返回多个最大元素对应的最小索引
     */
    public static int getMaxIndex(double[] a){
        return IntStream.range(0, a.length).boxed().max(Comparator.comparing(e->a[e])).get();
    }

    /**
     *
     * @param a
     * @return 数组最小元素的索引, 如果有多个最小元素, 则返回多个最小元素对应的最小索引
     */
    public static int getMinIndex(int[] a){
        return IntStream.range(0, a.length).boxed().min(Comparator.comparing(e->a[e])).get();
    }

    /**
     *
     * @param a
     * @return 数组最小元素的索引, 如果有多个最小元素, 则返回多个最小元素对应的最小索引
     */
    public static int getMinIndex(double[] a){
        return IntStream.range(0, a.length).boxed().min(Comparator.comparing(e->a[e])).get();
    }

}
