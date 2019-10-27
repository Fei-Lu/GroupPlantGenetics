package daxing.common;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public enum WheatLineage {
    A, B, D;

    /**
     * 返回A B D对应的ChrID
     * @param aOrBOrD A B D
     * @return 1,2, ...
     */
    public static int[] getWheatLineageOf(WheatLineage aOrBOrD){
        int[] a= Stream.concat(IntStream.iterate(1, n->n+6).limit(7).boxed(),
                IntStream.iterate(2, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
        int[] b=Stream.concat(IntStream.iterate(3, n->n+6).limit(7).boxed(),
                IntStream.iterate(4, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
        int[] d=Stream.concat(IntStream.iterate(5, n->n+6).limit(7).boxed(),
                IntStream.iterate(6, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
        if (aOrBOrD.toString().equals("A")){
            return a;
        }else if (aOrBOrD.toString().equals("B")){
            return b;
        }else {
            return d;
        }
    }

    /**
     *
     * @param aOrBOrD A B D
     * @return 1A, 2A, ...7A 或 1B, 2B, ...7B 或 1D, 2D, ...7D
     */
    public static List<String> wheatLineageOf(WheatLineage aOrBOrD){
        List<String> abd=new ArrayList<>();
        int[] chrA_Lineage= IntStream.range(1, 8).toArray();
        for (int i = 0; i < chrA_Lineage.length; i++) {
            abd.add(chrA_Lineage[i]+aOrBOrD.toString());
        }
        return abd;
    }

    /**
     *
     * @return 1,2,3,4,7,8, ...
     */
    public static int[] ablineage(){
        int[] a=WheatLineage.getWheatLineageOf(WheatLineage.A);
        int[] b=WheatLineage.getWheatLineageOf(WheatLineage.B);
        return Stream.concat(Arrays.stream(a).boxed(), Arrays.stream(b).boxed()).mapToInt(Integer::intValue).sorted().toArray();
    }

    /**
     *
     * @return 5,6,11,12, ...
     */
    public static int[]  dLineage(){
        return WheatLineage.getWheatLineageOf(WheatLineage.D);
    }

    /**
     *
     * @return 1A, 1B, 2A, 2B, 3A, 3B, ...
     */
    public static List<String> abLineage(){
        List<String> a=WheatLineage.wheatLineageOf(WheatLineage.A);
        List<String> b=WheatLineage.wheatLineageOf(WheatLineage.B);
        a.addAll(b);
        return a.stream().sorted().collect(Collectors.toList());
    }
}
