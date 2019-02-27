package daxing.barcode;

import analysis.pipeline.grt.TagUtils;
import com.koloboke.collect.map.hash.HashByteByteMap;
import format.dna.BaseEncoder;

/**
 *
 * @author xudaxing
 */
public class TagReads extends TagUtils{
    private static String polyA = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

    public static long[] getTagFromReads (String readR1, String readR2, HashByteByteMap ascIIByteMap, int tagLengthInLong) {
        int setReadLength = tagLengthInLong*BaseEncoder.longChunkSize;
        long[] tag = new long[tagLengthInLong*2];
        StringBuilder sb = new StringBuilder(readR1);
        if (sb.length()< setReadLength) {
            sb.append(polyA);
            readR1 = sb.substring(0, setReadLength);
        }
        sb = new StringBuilder(readR2);
        if (sb.length()<setReadLength) {
            sb.append(polyA);
            readR2 = sb.substring(0, setReadLength);
        }
        byte[] bArray = readR1.getBytes();
        for (int i = 0; i < bArray.length; i++) {
            bArray[i] = ascIIByteMap.get(bArray[i]);
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            tag[i] = BaseEncoder.getLongSeqFromSubByteArray(bArray, i*BaseEncoder.longChunkSize, (i+1)*BaseEncoder.longChunkSize);
        }
        bArray = readR2.getBytes();
        for (int i = 0; i < bArray.length; i++) {
            bArray[i] = ascIIByteMap.get(bArray[i]);
        }
        for (int i = 0; i < tagLengthInLong; i++) {
            tag[i+tagLengthInLong] = BaseEncoder.getLongSeqFromSubByteArray(bArray, i*BaseEncoder.longChunkSize, (i+1)*BaseEncoder.longChunkSize);
        }
        return tag;
    }

    public static String[] getReadsFromTag (long[] tag, short r1Length, short r2Length) {
        String[] reads = new String[2];
        int tagLengthInLong = tag.length/2;
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < tagLengthInLong; i++) {
            sb.append(BaseEncoder.getSequenceFromLong(tag[i]));
        }
        reads[0] = sb.toString().substring(0, r1Length);
        sb = new StringBuilder();
        for (int i = 0; i < tagLengthInLong; i++) {
            sb.append(BaseEncoder.getSequenceFromLong(tag[i+tagLengthInLong]));
        }
        reads[1] = sb.toString().substring(0, r2Length);
        return reads;
    }
}

