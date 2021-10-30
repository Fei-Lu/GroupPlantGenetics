package daxing.common.factors;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum P2 {

    CULTIVAR("CL", 0), LANDRACE("LR", 1);

    private final String abbreviation;
    private final int index;

    P2(String abbreviation, int index) {
        this.abbreviation=abbreviation;
        this.index=index;
    }

    public int getIndex() {
        return index;
    }

    public String getAbbreviation() {
        return abbreviation;
    }

    private static Map<String,P2> abbreviationToEnumMap= Stream.of(values()).collect(Collectors.toMap(P2::getAbbreviation, e->e));

    public static Optional<P2> getInstanceFrom(String abbreviation){
        return Optional.ofNullable(abbreviationToEnumMap.get(abbreviation));
    }
}
