package daxing.v2.localAncestryInfer;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum SourceType {


    WE(0, 0b1),
    DE(1, 0b10),
    FTT(2, 0b100),
    AT(3, 0b1000),
    NONE(4, 0b10000);

    private final int index;
    private final int feature;


    SourceType(int index, int feature) {
        this.index=index;
        this.feature=feature;
    }

    public int getIndex() {
        return index;
    }

    public int getSourceFeature() {
        return feature;
    }

    private static final Map<Integer, SourceType> indexToEnumMap =
            Stream.of(values()).collect(Collectors.toMap(SourceType::getIndex, e->e));

    public static Optional<SourceType> getInstanceFromIndex(int index){
        return Optional.ofNullable(indexToEnumMap.get(index));
    }

    private static final Map<Integer, SourceType> featureToEnumMap=
            Stream.of(values()).collect(Collectors.toMap(SourceType::getSourceFeature, e->e));

    public static Optional<SourceType> getInstanceFromFeature(int sourceFeature){
        return Optional.ofNullable(featureToEnumMap.get(sourceFeature));
    }

    public static int getSourceFeature(EnumSet<WindowSource.Source> sourceEnumSet){
        int feature=0;
        for (WindowSource.Source source: sourceEnumSet){
            feature = (SourceType.valueOf(source.name()).getSourceFeature())|feature;
        }
        return feature;
    }

    public static EnumSet<WindowSource.Source> getSourcesFrom(int sourceFeature){
        EnumSet<WindowSource.Source> sources= EnumSet.noneOf(WindowSource.Source.class);
        int feature;
        SourceType sourceType;
        for (int i = 1; i <= sourceFeature; i<<=1) {
            feature = sourceFeature & i;
            if (feature == 0) continue;
            sourceType = SourceType.getInstanceFromFeature(feature).get();
            sources.add(WindowSource.Source.valueOf(sourceType.name()));
        }
        if (sources.size()==0){
            System.out.println("check source feature!!!!");
        }
        return sources;
    }

    public static IntList getSingleSourceFeatureList(){
        IntList intList = new IntArrayList();
        for (int i = 0; i < SourceType.values().length; i++) {
            intList.add(SourceType.values()[i].getSourceFeature());
        }
        Collections.sort(intList);
        return intList;
    }
}
