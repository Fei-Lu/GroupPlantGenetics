package daxing.common.factors;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum LoadType {

    Syn("SYNONYMOUS",0), Non("NONSYNONYMOUS",1), Del("DELETERIOUS",2);

    private final String loadType;
    private final int index;

    LoadType(String loadType, int index) {
        this.loadType=loadType;
        this.index=index;
    }

    public String getLoadType() {
        return loadType;
    }

    public int getIndex() {
        return index;
    }

    private static final Map<String, LoadType> stringToEnumMap=
            Stream.of(values()).collect(Collectors.toMap(LoadType::getLoadType, e->e));

    public static Optional<LoadType> getInstanceFromString(String loadType){
        return Optional.ofNullable(stringToEnumMap.get(loadType));
    }
}
