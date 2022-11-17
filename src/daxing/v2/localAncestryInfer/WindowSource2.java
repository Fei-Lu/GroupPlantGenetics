package daxing.v2.localAncestryInfer;

import daxing.common.chrrange.ChrRange;

public class WindowSource2 implements Comparable<WindowSource2> {

    ChrRange chrRange;
    int sourceFuture;

    @Override
    public int compareTo(WindowSource2 o) {
        return this.chrRange.compareTo(o.chrRange);
    }

    public WindowSource2(ChrRange chrRange, SourceType sourceType){
        this.chrRange=chrRange;
        this.sourceFuture = sourceType.getSourceFeature();
    }

    public ChrRange getChrRange() {
        return chrRange;
    }

    public int getSourceFuture() {
        return sourceFuture;
    }

    public void setSourceFuture(SourceType sourceType){
        this.sourceFuture = sourceFuture | sourceType.getSourceFeature();
    }

    public boolean isSourceTypeOf(SourceType sourceType){
        if ((this.sourceFuture & sourceType.getSourceFeature()) == 0) return false;
        return true;
    }

    /**
     * if this source contain any source in the specified sourceFuture, return true, else return false
     * @param sourceFuture
     * @return
     */
    public boolean containAnySourceIn(int sourceFuture){
        if ((this.sourceFuture & sourceFuture) == 0) return false;
        return true;
    }


}
