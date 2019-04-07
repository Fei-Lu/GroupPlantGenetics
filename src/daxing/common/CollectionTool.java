package daxing.common;

import java.util.*;

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
}
