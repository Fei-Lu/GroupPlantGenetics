package daxing.common;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public enum WheatLineage {
    A, B, D;

    /**
     * 返回A B D对应的ChrID
     * @param aOrBOrD A B D
     * @return
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
}
