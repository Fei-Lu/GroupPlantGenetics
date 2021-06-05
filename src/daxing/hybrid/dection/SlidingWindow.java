package daxing.hybrid.dection;

import gnu.trove.list.array.TIntArrayList;

import java.util.Arrays;

/**
 * modified from SimpleWindow
 */
public class SlidingWindow {
    //incluside genomic position 1-based
    int[] starts = null;
    //exclusive genomic position 1-based
    int[] ends = null;
    double[] windowValues = null;
    int windowSize;
    int windowStep;
    int chrLength;

    public SlidingWindow (int chrLength, int windowSize, int windowStep) {
        this.chrLength = chrLength;
        this.windowSize = windowSize;
        this.windowStep = windowStep;
        this.initialize();
    }

    public SlidingWindow (int chrLength, int windowSize) {
        this.chrLength = chrLength;
        this.windowSize = windowSize;
        this.windowStep = windowSize;
        this.initialize();
    }

    /**
     * Add count of specified positions to windows
     * @param positions
     */
    public void addPositionCount(int[] positions) {
        for (int position : positions) {
            this.addPositionCount(position);
        }
    }

    /**
     * Add count of a specified position to windows
     * @param position
     */
    public void addPositionCount(int position) {
        int firstIndex;
        int lastIndex;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex; i++) {
            this.windowValues[i]++;
        }
    }
    /**
     * Add values of specified positions in windows
     * @param positions
     * @param values
     */
    public void addIntValues (int[] positions, int[] values) {
        int firstIndex;
        int lastIndex;
        for (int i = 0; i < positions.length; i++) {
            firstIndex = this.getFirstWindowIndex(positions[i]);
            lastIndex = this.getLastWindowIndex(positions[1]);
            for (int j = firstIndex; j < lastIndex; j++) {
                this.windowValues[j]+=values[i];
            }
        }
    }

    /**
     * Add value of a specified position in windows
     * @param position
     * @param value
     */
    public void addIntValue (int position, int value) {
        int firstIndex;
        int lastIndex;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex; i++) {
            this.windowValues[i]+=value;
        }
    }

    /**
     * Add values of specified positions in windows
     * @param positions
     * @param values
     */
    public void addFloatValues (int[] positions, float[] values) {
        int firstIndex;
        int lastIndex;
        for (int i = 0; i < positions.length; i++) {
            firstIndex = this.getFirstWindowIndex(positions[i]);
            lastIndex = this.getLastWindowIndex(positions[1]);
            for (int j = firstIndex; j < lastIndex; j++) {
                this.windowValues[j]+=values[i];
            }
        }
    }

    /**
     * Add values of specified positions in windows
     * @param position
     * @param value
     */
    public void addFloatValue (int position, float value) {
        int firstIndex;
        int lastIndex;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex+1; i++) {
            this.windowValues[i]+=value;
        }
    }

    /**
     * Add values of specified positions in windows
     * @param positions
     * @param values
     */
    public void addDoubleValues (int[] positions, double[] values) {
        int firstIndex;
        int lastIndex;
        for (int i = 0; i < positions.length; i++) {
            firstIndex = this.getFirstWindowIndex(positions[i]);
            lastIndex = this.getLastWindowIndex(positions[1]);
            for (int j = firstIndex; j < lastIndex; j++) {
                this.windowValues[j]+=values[i];
            }
        }
    }

    /**
     * Add values of specified positions in windows
     * @param position
     * @param value
     */
    public void addDoubleValue (int position, double value) {
        int firstIndex;
        int lastIndex;
        firstIndex = this.getFirstWindowIndex(position);
        lastIndex = this.getLastWindowIndex(position);
        for (int i = firstIndex; i < lastIndex; i++) {
            this.windowValues[i]+=value;
        }
    }

    /**
     * Return the values in windows
     * @return
     */
    public int[] getWindowValuesInt () {
        int[] values = new int[this.windowValues.length];
        for (int i = 0; i < this.windowValues.length; i++) {
            values[i] = (int)this.windowValues[i];
        }
        return values;
    }

    /**
     * Return the values in windows
     * @return
     */
    public float[] getWindowValuesFloat () {
        float[] values = new float[this.windowValues.length];
        for (int i = 0; i < this.windowValues.length; i++) {
            values[i] = (float)this.windowValues[i];
        }
        return values;
    }

    /**
     * Clear window values to 0
     */
    public void clearWindowValues () {
        this.windowValues = new double[this.starts.length];
    }

    /**
     * Return the values in windows
     * @return
     */
    public double[] getWindowValuesDouble () {
        return this.windowValues;
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
        endList.add(chrLength+1);
        this.starts = startList.toArray();
        this.ends = endList.toArray();
        this.windowValues = new double[ends.length];
    }
}
