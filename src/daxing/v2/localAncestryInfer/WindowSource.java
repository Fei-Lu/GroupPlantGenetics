package daxing.v2.localAncestryInfer;

import daxing.common.chrrange.ChrRange;
import java.util.EnumSet;

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

    public enum Source {
        WE(0), DE(1), FTT(2), AT(3), NONE(4);
        int index;
        Source(int index) {
            this.index=index;
        }

        public int getIndex() {
            return index;
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
