package daxing.common;

import java.util.*;
import java.util.stream.IntStream;

public class CollectionTool {

    public static <T> List<T> changeToList(T[] a){
        List<T> list=new ArrayList<>(a.length);
        Collections.addAll(list, a);
        return list;
    }

    public static Map<Integer, Integer> changeToIntMap(Map<String, String> map){
        Map<Integer, Integer> m=new HashMap<>();
        for (Map.Entry<String, String> entry: map.entrySet()){
            String key=entry.getKey();
            String value=entry.getValue();
            if (StringTool.isNumeric(key) && StringTool.isNumeric(value)){
                m.put(Integer.valueOf(key), Integer.valueOf(value));
            }
        }
        return m;
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
