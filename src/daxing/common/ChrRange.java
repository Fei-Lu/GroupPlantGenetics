package daxing.common;

import pgl.infra.range.Range;
import pgl.infra.utils.wheat.RefV1Utils;

import java.util.Objects;

/**
 * @author Daxing Xu
 */
public class ChrRange implements Comparable<ChrRange>{

    String chr;
    int start;
    /**
     * Exclusive
     */
    int end;

    public ChrRange(String chr, int start, int end){
        this.chr=chr;
        this.start=start;
        this.end=end;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getLength(){
        return end-start;
    }

    public void setStart(int startOnRefChr){
        this.start=startOnRefChr;
    }

    public void setEnd(int endOnRefChr){
        this.end=endOnRefChr;
    }

    public int getChrID(){
        return RefV1Utils.getChrID(chr, start);
    }

    public int getVCFStart(){
        return RefV1Utils.getPosOnChrID(chr, start);
    }

    public int getVCFEnd(){
        return RefV1Utils.getPosOnChrID(chr, end);
    }

    public short getVCFStartChrID(){
        return (short) RefV1Utils.getChrID(chr, start);
    }

    public short getVCFEndChrID(){
        return (short) RefV1Utils.getChrID(chr, end-1);
    }

    public boolean contain(String chr, int pos){
        if (!this.getChr().equals(chr)) return false;
        if (pos < start) return false;
        return pos < end;
    }

    /**
     * Daxing's utils on ChrPos
     * @param chrPos chrPos
     * @return true if contain chrPos
     */
    public boolean contain(ChrPos chrPos){
        return contain(chrPos.chr, chrPos.pos);
    }

    @Override
    public int compareTo(ChrRange o) {
//        if (this.chr.compareToIgnoreCase(o.chr)==0){
//            if (this.start==o.start){
//                if (this.end==o.end) return 0;
//                if (this.end > o.end) return 1;
//                return -1;
//            }else if (this.start > o.start) return 1;
//            return -1;
//        }else if (this.chr.compareToIgnoreCase(o.chr) > 0) return 1;
//        return -1;
        if (this.chr.compareToIgnoreCase(o.chr)==0){
            return Integer.compare(this.start, o.start);
        }else if (this.chr.compareToIgnoreCase(o.chr) > 0) return 1;
        return -1;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        ChrRange chrRange = (ChrRange) o;
        return start == chrRange.start &&
                end == chrRange.end &&
                chr.equals(chrRange.chr);
    }

    @Override
    public int hashCode() {
        return Objects.hash(chr, start, end);
    }

    public String toString(){
        return this.getChr() + ":" + this.getStart() + "," + this.getEnd();
    }

    public boolean isOverlapped(ChrRange chrRange){
        if (!this.chr.equals(chrRange.chr)) return false;
        if (this.getStart() >= chrRange.getEnd()) return false;
        return this.getEnd() > chrRange.getStart();
    }

    public ChrRange getIntersection(ChrRange chrRange){
        if (this.isOverlapped(chrRange)) return null;
        if (this.getStart() > chrRange.getStart()) return new ChrRange(this.getChr(), this.getStart(), Math.min(this.getEnd(), chrRange.getEnd()));
        return new ChrRange(this.getChr(), chrRange.getStart(), Math.min(this.getEnd(), chrRange.getEnd()));
    }

    public static ChrRange changeToChrRange(Range range){
        int chrID=range.getRangeChromosome();
        int startOnChrID=range.getRangeStart();
        int endOnChrID=range.getRangeEnd();
        String chr=RefV1Utils.getChromosome(chrID, startOnChrID);
        int startOnRefChr=RefV1Utils.getPosOnChromosome(chrID, startOnChrID);
        int endOnRefChr=RefV1Utils.getPosOnChromosome(chrID, endOnChrID);
        return new ChrRange(chr, startOnRefChr, endOnRefChr);
    }
}
