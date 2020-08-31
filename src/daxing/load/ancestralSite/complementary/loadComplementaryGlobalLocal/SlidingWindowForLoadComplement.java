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
    TDoubleArrayList[] windowValue1List =null;
    TDoubleArrayList[] windowValue2List =null;
    TDoubleArrayList[] windowValue3List =null;
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

    public TDoubleArrayList getWindow1Value(int windowIndex){
        return this.windowValue1List[windowIndex];
    }

    public TDoubleArrayList getWindow2Value(int windowIndex){
        return this.windowValue2List[windowIndex];
    }

    public TDoubleArrayList getWindow3Value(int windowIndex){
        return this.windowValue3List[windowIndex];
    }

    public int getCountInWindow(int windowIndex){
        return this.getWindow1Value(windowIndex).size();
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
     * @param value123
     */
    public void addValue(int position, double[] value123){
        int firstIndex = -1;
        int lastIndex = -1;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex; i++) {
            this.windowValue1List[i].add(value123[0]);
            this.windowValue2List[i].add(value123[1]);
            this.windowValue3List[i].add(value123[2]);
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
        this.windowValue1List =new TDoubleArrayList[ends.length];
        this.windowValue2List =new TDoubleArrayList[ends.length];
        this.windowValue3List =new TDoubleArrayList[ends.length];
        for (int i = 0; i < this.windowValue1List.length; i++) {
            this.windowValue1List[i]=new TDoubleArrayList();
            this.windowValue2List[i]=new TDoubleArrayList();
            this.windowValue3List[i]=new TDoubleArrayList();
        }
    }
}