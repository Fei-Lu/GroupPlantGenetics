package daxing.common;

import gnu.trove.list.array.TIntArrayList;

public class RandomAccessFileTool {

    public static long getPointer(TIntArrayList everyLineLengthInFile, int rowIndex){
        long pointPosition=0;
        for (int i = 0; i < rowIndex; i++) {
            pointPosition+=everyLineLengthInFile.get(i)+1;
        }
        return pointPosition;
    }


}
