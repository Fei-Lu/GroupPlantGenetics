package xiaohan.utils;

import java.util.List;
import java.util.Set;

public class CollectionsCalculator {


    public static Set getIntersection(Set A, Set B){
        A.retainAll(B);
        return A;
    }

    public static Set getComplement(Set A, Set B){
        A.removeAll(B);
        return A;
    }

    public static Set getUnionall(Set A, Set B){
        A.addAll(B);
        return A;
    }

    public static Set getUnionsorted(Set A, Set B){
        A.removeAll(B);
        A.addAll(B);
        return A;
    }

    public static Set[] getSubtraction(Set A, Set B){
        Set Acp = A;
        A.removeAll(B);
        B.removeAll(Acp);
        Set[] sets = new Set[2];
        sets[0] = A;
        sets[1] = B;
        return sets;
    }

    public static Set getSymmetricSubtraction(Set A, Set B){
        A.removeAll(B);
        return A;
    }

    public static List getIntersection(List A, List B){
        A.retainAll(B);
        return A;
    }

    public static List getComplement(List A, List B){
        A.removeAll(B);
        return A;
    }

    public static List getUnionall(List A, List B){
        A.addAll(B);
        return A;
    }

    public static List getUnionsorted(List A, List B){
        A.removeAll(B);
        A.addAll(B);
        return A;
    }

    public static List[] getSubtraction(List A, List B){
        List Acp = A;
        A.removeAll(B);
        B.removeAll(Acp);
        List[] Lists = new List[2];
        Lists[0] = A;
        Lists[1] = B;
        return Lists;
    }

    public static List getSymmetricSubtraction(List A, List B){
        A.removeAll(B);
        return A;
    }
}
