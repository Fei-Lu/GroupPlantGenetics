package daxing.common.utiles;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Daxing Xu
 */
public class CollectionTool {

    /**
     *
     * @param a
     * @param <T>
     * @return
     */
    public static <T> List<T> changeToList(T[] a){
        List<T> list=new ArrayList<>(a.length);
        Collections.addAll(list, a);
        return list;
    }

    /**
     *
     * @param list
     * @param <T>
     * @return the number of occurrences of each element in the list
     */
    public static<T> Map<T, Long> caculateElementCount(List<T> list){
        return list.stream().collect(Collectors.groupingBy(t->t, Collectors.counting()));
    }

    public static boolean hasIntersection(IntSet intSet1, IntSet intSet2){
        IntSet intSet = new IntOpenHashSet(intSet1);
        intSet.retainAll(intSet2);
        return intSet.size() > 0 ? true : false;
    }

}
