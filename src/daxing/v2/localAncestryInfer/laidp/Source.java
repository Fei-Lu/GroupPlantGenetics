package daxing.v2.localAncestryInfer.laidp;

import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum Source {
    NATIVE_SOURCE_0(0, 0b1),
    INTROGRESSED_SOURCE_1(1, 0b10),
    INTROGRESSED_SOURCE_2(2, 0b100),
    INTROGRESSED_SOURCE_3(3, 0b1000),
    INTROGRESSED_SOURCE_4(4, 0b10000);

    private int index;
    private int feature;

    Source(int index, int feature) {
        this.index=index;
        this.feature=feature;
    }

    public int getFeature() {
        return feature;
    }

    public int getIndex() {
        return index;
    }

    public static int addSourceFeature(int feature0, int feature1){
        return feature0 | feature1;
    }

    public void setSourceFeature(Source source){
        int resFeature = this.feature | source.feature;
        this.feature=resFeature;
    }

    private static final Map<Integer, Source> indexToEnumMap =
            Stream.of(values()).collect(Collectors.toMap(Source::getIndex, e->e));

    public static Optional<Source> getInstanceFromIndex(int index){
        return Optional.ofNullable(indexToEnumMap.get(index));
    }


}
