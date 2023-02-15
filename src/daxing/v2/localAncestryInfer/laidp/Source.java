package daxing.v2.localAncestryInfer.laidp;

import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import java.util.Collections;
import java.util.EnumSet;
import java.util.Map;
import java.util.Optional;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public enum Source {
    NATIVE(0, 0b1),
    INTROGRESSED_1(1, 0b10),
    INTROGRESSED_2(2, 0b100),
    INTROGRESSED_3(3, 0b1000),
    INTROGRESSED_4(4, 0b10000),
    ADMIXED(5, 0b100000);

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

    private static final Map<Integer, Source> featureToEnumMap=
            Stream.of(values()).collect(Collectors.toMap(Source::getFeature, e->e));

    public static Optional<Source> getInstanceFromFeature(int sourceFeature){
        return Optional.ofNullable(featureToEnumMap.get(sourceFeature));
    }

    public static int getSourceFeature(EnumSet<Source> sourceEnumSet){
        int feature=0;
        for (Source source: sourceEnumSet){
            feature = (source.feature)|feature;
        }
        return feature;
    }

    public static EnumSet<Source> getSourcesFrom(int sourceFeature){
        EnumSet<Source> sources= EnumSet.noneOf(Source.class);
        int feature;
        Source source;
        for (int i = 1; i <= sourceFeature; i<<=1) {
            feature = sourceFeature & i;
            if (feature == 0) continue;
            source = Source.getInstanceFromFeature(feature).get();
            sources.add(Source.valueOf(source.name()));
        }
        if (sources.size()==0){
            System.out.println("check source feature!!!!");
        }
        return sources;
    }

    public static IntList getSingleSourceFeatureList(){
        IntList intList = new IntArrayList();
        for (int i = 0; i < Source.values().length; i++) {
            intList.add(Source.values()[i].getFeature());
        }
        Collections.sort(intList);
        return intList;
    }


}
