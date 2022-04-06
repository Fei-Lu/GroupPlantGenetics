package daxing.common.factors;

import pgl.infra.utils.wheat.RefV1Utils;

import java.util.Arrays;

public enum SubgenomeCombination {
    ABD, AB, D;

    public int[] getChrID(){
        switch (this){
            case ABD:
                return RefV1Utils.getChrIDs();
            case AB:
                return WheatLineage.ablineage();
            case D:
                return RefV1Utils.getChrIDsOfSubgenomeD();
        }
        return null;
    }

    public boolean containChrID(int chrID){
        int[] chrIDArray = getChrID();
        int index = Arrays.binarySearch(chrIDArray, chrID);
        if (index < 0) return false;
        return true;
    }
}
