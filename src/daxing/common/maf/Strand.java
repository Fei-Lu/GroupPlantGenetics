package daxing.common.maf;

import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

public enum Strand {

    Forward("+"), Reverse("-");

    String str;

    Strand(String str){
        this.str=str;
    }

    public String getStr() {
        return str;
    }

    private static final Map<String, Strand> strandMap= Arrays.stream(values()).collect(Collectors.toMap(Strand::getStr, e->e));

    public static Strand getInstanceFrom(String strand){
        return strandMap.get(strand);
    }
}
