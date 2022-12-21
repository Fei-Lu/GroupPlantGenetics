package daxing.v2.localAncestryInfer;

import daxing.common.chrrange.ChrRange;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * Record the windows and source population
 */
public class WindowSource implements Comparable<WindowSource>{

    ChrRange chrRange;
    EnumSet<Source> sources;

    @Override
    public int compareTo(WindowSource o) {
        return this.chrRange.compareTo(o.chrRange);
    }

    public enum Source{
        WE(0), DE(1), FTT(2), AT(3), NONE(4);
        int index;
        Source(int index) {
            this.index=index;
        }

        public int getIndex() {
            return index;
        }

        private static Map<Integer, Source> index2SourceMap=
                Stream.of(values()).collect(Collectors.toMap(Source::getIndex, e->e));

        public static Optional<Source> getInstanceFromSubNum(int sourceIndex){
            return Optional.ofNullable(index2SourceMap.get(sourceIndex));
        }

        public static EnumSet<Source> getABSource(){
            return EnumSet.of(WE,DE,FTT,NONE);
        }

        public static EnumSet<Source> getDSource(){
            return EnumSet.of(AT,NONE);
        }
    }

    public WindowSource(ChrRange chrRange, EnumSet<Source> sources){
        this.chrRange=chrRange;
        this.sources =sources;
    }

    public ChrRange getChrRange() {
        return chrRange;
    }

    public EnumSet<Source> getSources() {
        return sources;
    }

    public Set<Source> getSourceSet(){
        Set<Source> set = new HashSet<>();
        EnumSet<Source> enumSet = this.getSources();
        for (Source source: enumSet){
            set.add(source);
        }
        return set;
    }

    public int getSourceNum(){
        return this.sources.size();
    }

    public boolean containAnySourceIn(Source source){
        return sources.contains(source);
    }

    /**
     * if this source contain any source in the specified collection, return true, else return false
     * @param sources
     * @return
     */
    public boolean containAnySourceIn(EnumSet<Source> sources){
        EnumSet<Source> thisSource=this.sources;
        sources.retainAll(thisSource);
        if (sources.size() > 0){
            return true;
        }
        return false;
    }

}
