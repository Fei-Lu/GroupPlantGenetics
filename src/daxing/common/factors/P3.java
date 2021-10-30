package daxing.common.factors;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum P3 {

    WILD_EMMER("WE",0), DOMESTICATED_EMMER("DE",1), FREE_THRESHING_TETRAPLOIDS("FTT",2), STRANGULATA("AT",3);

    private final String abbreviation;
    private final int index;

    P3(String p3Abbreviation, int index) {
        this.abbreviation=p3Abbreviation;
        this.index=index;
    }

    public int getIndex() {
        return index;
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    private static final Map<String, P3> abbreviationToEnumMap= Stream.of(values()).collect(Collectors.toMap(P3::getAbbreviation, e->e));

    public static Optional<P3> getInstanceFrom(String p3Abbreviation){
        return Optional.ofNullable(abbreviationToEnumMap.get(p3Abbreviation));
    }

}
