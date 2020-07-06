package daxing;

import daxing.selection.model.BlockOperation;

public class Start {

    public static void main(String[] args) {
        String triadsMafFile="/Users/xudaxing/Desktop/test/test.maf";
        String outFile="/Users/xudaxing/Desktop/test/res.txt";
        BlockOperation.initialize(triadsMafFile, outFile);
    }
}