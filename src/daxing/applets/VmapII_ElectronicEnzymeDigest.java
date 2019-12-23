package daxing.applets;

import daxing.common.SeqByte;
import gnu.trove.list.array.TByteArrayList;
import utils.IOUtils;

import java.io.BufferedReader;
import java.io.IOException;

public class VmapII_ElectronicEnzymeDigest {

    public static void digest(String chromosomeFile, String cutter1, String cutter2){
        TByteArrayList chrSeq=new TByteArrayList(1_000_000_000);
        try (BufferedReader bufferedReader = IOUtils.getTextReader(chromosomeFile)) {
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

    }
}
