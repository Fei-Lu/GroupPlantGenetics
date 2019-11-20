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
     *
     * @return 1,2,7,8, ... or 3,4,9,10, .. or 5,6,11,12, ...
     */
    public int[] getChrID(){
        int[] chrID=null;
        switch (this){
            case A:
                chrID=Stream.concat(IntStream.iterate(1, n->n+6).limit(7).boxed(),
                        IntStream.iterate(2, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
                break;
            case B:
                chrID=Stream.concat(IntStream.iterate(3, n->n+6).limit(7).boxed(),
                        IntStream.iterate(4, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
                break;
            case D:
                chrID=Stream.concat(IntStream.iterate(5, n->n+6).limit(7).boxed(),
                        IntStream.iterate(6, n->n+6).limit(7).boxed()).sorted().mapToInt(Integer::intValue).toArray();
                break;
            default:
                System.out.println("please input 'A' or 'B' or 'D'");

        }
        return chrID;
    }

    /**
     *
     * @return 1A, 2A, ...7A or 1B, 2B, ...7B or 1D, 2D, ...7D
     */
    public List<String> getChr(){
        List<String> chrs=new ArrayList<>();
        int[] chrID1_7= IntStream.range(1, 8).toArray();
        for (int i = 0; i < chrID1_7.length; i++) {
            chrs.add(chrID1_7[i]+this.name());
        }
        return chrs;
    }

    /**
     *
     * @return 1,2,3,4,7,8,9,10, ...
     */
    public static int[] ablineage(){
        int[] a=WheatLineage.valueOf("A").getChrID();
        int[] b=WheatLineage.valueOf("B").getChrID();
        return Stream.concat(Arrays.stream(a).boxed(), Arrays.stream(b).boxed()).mapToInt(Integer::intValue).sorted().toArray();
    }

    /**
     *
     * @return 5,6,11,12, ...
     */
    public static int[] dlineage(){
        return WheatLineage.valueOf("D").getChrID();
    }

    /**
     *
     * @return 1A, 1B, 2A, 2B, 3A, 3B, ...
     */
    public static List<String> abLineage(){
        List<String> a=WheatLineage.valueOf("A").getChr();
        List<String> b=WheatLineage.valueOf("B").getChr();
        a.addAll(b);
        return a.stream().sorted().collect(Collectors.toList());
    }

    /**
     *
     * @return 1D, 2D, 3D, ...
     */
    public static List<String> dLineage(){
        return WheatLineage.valueOf("D").getChr();
    }

    /**
     *
     * @return 1A, 1B, 1D, 2A, 2B, 2D, ...
     */
    public static List<String> abdLineage(){
        List<String> ab=WheatLineage.abLineage();
        List<String> d=WheatLineage.dLineage();
        ab.addAll(d);
        return ab.stream().sorted().collect(Collectors.toList());
    }
}
