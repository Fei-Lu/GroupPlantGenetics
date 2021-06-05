package daxing.applets;

import daxing.common.IOTool;
import daxing.common.SeqByte;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;

import java.io.BufferedReader;
import java.io.IOException;

public class VmapII_ElectronicEnzymeDigest {

    public static int digest(String chromosomeFile, String cutter1, String cutter2){
        TByteArrayList chrSeq=new TByteArrayList(1_000_000_000);
        try (BufferedReader bufferedReader = IOTool.getReader(chromosomeFile)) {
            String line;
            bufferedReader.readLine();
            SeqByte seqByte;
            while ((line= bufferedReader.readLine())!=null){
                seqByte=new SeqByte(line);
                chrSeq.add(seqByte.getSeqByte());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        SeqByte sequenceByte=new SeqByte(chrSeq.toArray());
        TIntArrayList indexCutter1 = sequenceByte.indexOfAll(cutter1.getBytes());
        TIntArrayList indexCutter2 = sequenceByte.indexOfAll(cutter2.getBytes());
        int count=0;
        for (int i = 0; i < indexCutter1.size(); i++) {
            for (int j = 0; j < indexCutter2.size(); j++) {
                if (indexCutter1.get(i) > indexCutter2.get(j)) continue;
                if (indexCutter2.get(j) - indexCutter1.get(i) > 500) break;
                count++;
            }
        }
        return count;
    }
}
