package daxing.load.aoyue;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

public class Aoyue {
//    String [] key = {"M011","M110"};
//        Aoyue.sort(key);
//        Aoyue.models7(key);
    public static String models7(String[] input){
        String out=null;

        HashMap<String,String> hm = new HashMap<>();
        String[] value = {"#fd8582","#ffcbcc","#967bce","#cccdfe","darkgray","#4bcdc6","#cdfffc"};
        String[] key = {"M100","M011","M010","M101","M111","M001","M110"};
        for (int i = 0; i < value.length; i++) {
            hm.put(key[i],value[i]);
        }

        StringBuilder sb = new StringBuilder();
        sb.append("c(");
        for (int i = 0; i < input.length; i++) {
            String col = hm.get(input[i]);
            sb.append("'").append(col).append("',");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append(")");
        out = sb.toString();
        System.out.println(out);

        return out;
    }

    /**
     * 将因子变量进行自然顺序的排序，然后再调用颜色值
     * @param input
     */
    public static void sort(String[] input){
//        String[] order = {"M001","M010","M011", "M100", "M101", "M110", "M111"};
        List<String> l = new ArrayList<>();
        for (int i = 0; i < input.length; i++) {
            l.add(input[i]);
        }
        Collections.sort(l);
        for (int i = 0; i < l.size(); i++) {
            input[i] = l.get(i);
        }
    }
}
