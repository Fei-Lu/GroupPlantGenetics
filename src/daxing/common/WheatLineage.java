package daxing.common;

import pgl.infra.utils.PStringUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

/**
 * @author Daxing Xu
 */
public enum WheatLineage {

    A(0), B(1), D(2);

    int index;
    WheatLineage(int index) {
        this.index=index;
    }

    public int getIndex() {
        return index;
    }

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
     * @param chrID_or_Chr Chr or ChrID
     * @return
     */
    public Predicate<File> getPredicate(String chrID_or_Chr){
        List<String> type=Stream.of("ChrID", "Chr").collect(Collectors.toList());
        if (!type.contains(chrID_or_Chr)){
            System.out.println("error, check your chrID_or_Chr, it must be ChrID or Chr");
            System.exit(1);
        }
        Predicate<File> p=null;
        List<String> chrs;
        List<String> chrID;
        if (chrID_or_Chr.equals("Chr")){
            chrs=this.getChr();
            p=f->chrs.contains(f.getName().substring(3,5));
        }else {
            chrID=Arrays.stream(this.getChrID()).boxed().map(str-> PStringUtils.getNDigitNumber(3, str)).collect(Collectors.toList());
            p=f->chrID.contains(f.getName().substring(3,6));
        }
        return p;
    }

    /**
     *
     * @param chrID_or_Chr
     * @return
     */
    public static Predicate<File> getABPredicate(String chrID_or_Chr){
        List<String> type=Stream.of("ChrID", "Chr").collect(Collectors.toList());
        if (!type.contains(chrID_or_Chr)){
            System.out.println("error, check your chrID_or_Chr, it must be ChrID or Chr");
            System.exit(1);
        }
        Predicate<File> p=null;
        List<String> chrs=WheatLineage.abLineage();
        List<String> chrID=Arrays.stream(WheatLineage.ablineage()).boxed().map(s->String.valueOf(s)).collect(Collectors.toList());
        if (chrID_or_Chr.equals("Chr")){
            p=f->chrs.contains(f.getName().substring(3,5));
        }else {
            p=f->chrID.contains(f.getName().substring(3,6));
        }
        return p;
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
     * @return 1A, 1B, 1D, 2A, 2B, 2D, ...
     */
    public static List<String> abdLineage(){
        List<String> ab=WheatLineage.abLineage();
        List<String> d=WheatLineage.valueOf("D").getChr();
        ab.addAll(d);
        return ab.stream().sorted().collect(Collectors.toList());
    }
}
