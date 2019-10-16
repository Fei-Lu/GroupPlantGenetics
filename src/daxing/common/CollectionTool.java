package daxing.common;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class CollectionTool<T> {

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

    /**
     *
     * @param aOrbOrd A B D
     * @return 1A, 2A, ...7A 或 1B, 2B, ...7B 或 1D, 2D, ...7D
     */
    public static List<String> wheatLineageOf(WheatLineage aOrbOrd){
        List<String> abd=new ArrayList<>();
        int[] chrA_Lineage= IntStream.range(1, 8).toArray();
        for (int i = 0; i < chrA_Lineage.length; i++) {
            abd.add(chrA_Lineage[i]+aOrbOrd.toString());
        }
        return abd;
    }
}
