package daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.util.Arrays;

/**
 * @author xudaxing
 */
public class SlidingWindowForLoadComplement {
    //incluside
    int[] starts = null;
    //exclusive
    int[] ends = null;
    TDoubleArrayList[] windowValueSynList=null;
    TDoubleArrayList[] windowValueNonList=null;
    TDoubleArrayList[] windowValueDelList=null;
    int windowSize = Integer.MIN_VALUE;
    int windowStep = Integer.MIN_VALUE;
    int chrLength = Integer.MIN_VALUE;

    public SlidingWindowForLoadComplement(int chrLength, int windowSize, int windowStep) {
        this.chrLength = chrLength;
        this.windowSize = windowSize;
        this.windowStep = windowStep;
        this.initialize();
    }

    /**
     * Return the starting positions of windows, inclusive
     * @return
     */
    public int[] getWindowStarts () {
        return this.starts;
    }

    /**
     * Return the ending positions of windows, exclusive
     * @return
     */
    public int[] getWindowEnds () {
        return this.ends;
    }

    public int getWindowNum(){
        return this.ends.length;
    }

    public TDoubleArrayList getSynWindowValue(int windowIndex){
        return this.windowValueSynList[windowIndex];
    }

    public TDoubleArrayList getNonWindowValue(int windowIndex){
        return this.windowValueNonList[windowIndex];
    }

    public TDoubleArrayList getDelWindowValue(int windowIndex){
        return this.windowValueDelList[windowIndex];
    }

    public int getGeneNum(int windowIndex){
        return this.getSynWindowValue(windowIndex).size();
    }


    /**
     * Return the first index of first window containing the position, inclusive
     * @param position
     * @return
     */
    public int getFirstWindowIndex (int position) {
        int index = Arrays.binarySearch(this.ends, position);
        if (index < 0) index = -index-1;
        else index++;
        return index;
    }

    /**
     * Return the last index of first window containing the position, exclusive
     * @param position
     * @return
     */
    public int getLastWindowIndex (int position) {
        int index = Arrays.binarySearch(this.starts, position);
        if (index < 0) index = -index-2;
        return index+1;
    }

    /**
     *
     * @param position
     * @param synNonDelValue
     */
    public void addValue(int position, double[] synNonDelValue){
        int firstIndex = -1;
        int lastIndex = -1;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex; i++) {
            this.windowValueSynList[i].add(synNonDelValue[0]);
            this.windowValueNonList[i].add(synNonDelValue[1]);
            this.windowValueDelList[i].add(synNonDelValue[2]);
        }
    }

    private void initialize () {
        TIntArrayList startList = new TIntArrayList();
        TIntArrayList endList = new TIntArrayList();
        int start = 1;
        int end = start+windowSize;
        while (end < chrLength) {
            startList.add(start);
            endList.add(end);
            start+=windowStep;
            end = start+windowSize;
        }
        startList.add(start);
        endList.add(end);
        this.starts = startList.toArray();
        this.ends = endList.toArray();
        this.windowValueSynList=new TDoubleArrayList[ends.length];
        this.windowValueNonList=new TDoubleArrayList[ends.length];
        this.windowValueDelList=new TDoubleArrayList[ends.length];
        for (int i = 0; i < this.windowValueSynList.length; i++) {
            this.windowValueSynList[i]=new TDoubleArrayList();
            this.windowValueNonList[i]=new TDoubleArrayList();
            this.windowValueDelList[i]=new TDoubleArrayList();
        }
    }
}