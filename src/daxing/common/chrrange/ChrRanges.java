package daxing.common.chrrange;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import gnu.trove.list.array.TIntArrayList;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class ChrRanges {

    List<ChrRange> chrRanges;
    SortType sortType;
    public static final ChrRanges blankChrRanges=new ChrRanges();

    Swapper swapper = (indexA, indexB) -> {
        ChrRange temp = getChrRange(indexA);
        setChrRange(indexA, getChrRange(indexB));
        setChrRange(indexB, temp);
    };

    IntComparator compBySize = (indexA, indexB) ->{
        int sA = getChrRange(indexA).getLength();
        int sB = getChrRange(indexB).getLength();
        return sA-sB;
    };

    IntComparator compByPosition = (indexA, indexB) ->{
        String chrA = getChrRange(indexA).getChr();
        String chrB = getChrRange(indexB).getChr();
        if (chrA.equals(chrB)) {
            int startA = getChrRange(indexA).getStart();
            int startB = getChrRange(indexB).getStart();
            return startA-startB;
        }
        else {
            return chrA.compareTo(chrB);
        }
    };

    public enum SortType{
        SIZE, POSITION, UNSORTED
    }

    public ChrRanges(List<ChrRange> chrRangeList){
        this.chrRanges=chrRangeList;
        this.sortType=SortType.UNSORTED;
    }

    private ChrRanges(){
        this.chrRanges=new ArrayList<>();
        this.sortType=SortType.UNSORTED;
    }

    public List<ChrRange> getChrRanges() {
        return chrRanges;
    }

    public ChrRange getChrRange(int rangeIndex){
        return this.chrRanges.get(rangeIndex);
    }

    public void setSortType(SortType sortType){
        this.sortType=sortType;
    }

    /**
     *
     * @param rangeIndex
     * @param chrRange
     * @return the ChrRange previously at the specified position
     */
    public ChrRange setChrRange(int rangeIndex, ChrRange chrRange){
        ChrRange chrRange1=this.chrRanges.set(rangeIndex, chrRange);
        this.sortType=SortType.UNSORTED;
        return chrRange1;
    }

    /**
     *
     * @param rangeIndex
     * @return the ChrRange previously at the specified position
     */
    public ChrRange removeChrRange(int rangeIndex){
        ChrRange chrRange=this.chrRanges.remove(rangeIndex);
        return  chrRange;
    }

    public void addChrRange(int rangeIndex, ChrRange chrRange){
        this.chrRanges.add(rangeIndex, chrRange);
    }

    public int getRangeNumber(){
        return this.chrRanges.size();
    }

    public void sortBy(SortType sortType){
        assert sortType.equals(SortType.SIZE) || sortType.equals(SortType.POSITION) : "must sort by position or size," +
                " other do not support";
        switch (sortType){
            case SIZE:
                GenericSorting.quickSort(0, this.getRangeNumber(), compBySize, swapper);
                this.sortType = SortType.SIZE;
                break;
            case POSITION:
                GenericSorting.quickSort(0, this.getRangeNumber(), compByPosition, swapper);
                this.sortType = SortType.POSITION;
                break;
            case UNSORTED:
                break;
        }
    }

    public boolean contain(String refChr, int refPos){
        assert sortType==(SortType.POSITION) : "please sort by start position";
        ChrRange query = new ChrRange(refChr, refPos, refPos+1);
        int hit = Collections.binarySearch(this.getChrRanges(), query);
        hit = hit < -1 ? -hit-2 : hit;
        if (hit < 0) return false;
        TIntArrayList indexList=new TIntArrayList();
        for (int i = hit; i > -1 ; i--) {
            if (this.getChrRange(i).contain(refChr, refPos)) indexList.add(i);
        }
        if (indexList.size()==0) return false;
        return true;
    }

    /**
     *
     * @param refChr ref chr
     * @return negative if do not contain
     */
    public int getStartIndexOfChromosome(String refChr){
        assert this.sortType==SortType.POSITION : "please sort by position";
        ChrRange query= new ChrRange(refChr, Integer.MIN_VALUE, Integer.MIN_VALUE);
        int hit  = Collections.binarySearch(this.getChrRanges(), query);
        if (hit < 0) {
            int index = -hit-1;
            if (this.getChrRange(index).getChr().equals(refChr)) return index;
            return hit;
        }
        return hit;
    }

    /**
     * refChr的长度必须小于Integer.MAX_VALUE
     * @param refChr ref chr
     * @return negative if do not contain
     */
    public int getEndIndexOfChromosome(String refChr){
        assert this.sortType==SortType.POSITION : "please sort by position";
        ChrRange query = new ChrRange(refChr, Integer.MAX_VALUE, Integer.MAX_VALUE);
        int hit  = Collections.binarySearch(this.getChrRanges(), query);
        if (hit < -1) {
            int index = -hit-2;
            if (this.getChrRange(index).getChr().equals(refChr)) return index;
            return hit;
        }
        return hit;
    }

    /**
     *
     * @param refChr ref chr
     * @param refPos ref pos
     * @return new int[0] if do not contain
     */
    public int[] getIndicesContainsPosition (String refChr, int refPos) {
        assert sortType==(SortType.POSITION) : "please sort by start position";
        TIntArrayList indexList = new TIntArrayList();
        ChrRange query = new ChrRange(refChr, refPos, refPos+1);
        int hit = Collections.binarySearch(this.getChrRanges(), query);
        hit = hit < -1 ? -hit-2 : hit;
        if (hit < 0) return new int[0];
        for (int i = hit; i > -1; i--) {
            if (this.getChrRange(i).contain(refChr, refPos)) indexList.add(i);
        }
        indexList.sort();
        return indexList.toArray();
    }

    /**
     *
     * @param refChr ref chr
     * @param refPos ref position
     * @return null if do not contain
     */
    public ChrRanges getChrRangesContainsPosition(String refChr, int refPos){
        int[] indices=this.getIndicesContainsPosition(refChr, refPos);
        if (indices.length==0) return ChrRanges.blankChrRanges;
        List<ChrRange> chrRangeList=new ArrayList<>();
        for (int index : indices){
            chrRangeList.add(this.getChrRange(index));
        }
        return new ChrRanges(chrRangeList);
    }

    public ChrRanges getChrRangesByChromosome(String refChr){
        int startIndex = this.getStartIndexOfChromosome(refChr);
        if (startIndex < 0) return ChrRanges.blankChrRanges;
        int endIndex = this.getEndIndexOfChromosome(refChr);
        if (endIndex < 0) return ChrRanges.blankChrRanges;
        List<ChrRange> l = new ArrayList<>();
        for (int i = startIndex; i < endIndex+1; i++) {
            l.add(this.getChrRange(i));
        }
        return new ChrRanges(l);
    }

}
