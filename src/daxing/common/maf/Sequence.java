package daxing.common.maf;

import com.koloboke.collect.map.hash.HashByteCharMap;
import com.koloboke.collect.map.hash.HashByteCharMaps;
import com.koloboke.collect.map.hash.HashCharByteMap;
import com.koloboke.collect.map.hash.HashCharByteMaps;

public class Sequence {

    /**
     * Bases in char
     */
    public static final char[] base = {'-', 'A', 'C', 'G', 'T'};

    /**
     * Bases in AscII
     */
    public static final byte[] baseAscII={45, 65, 67, 71, 84};

    public static final HashCharByteMap baseToAscIIMap =
            HashCharByteMaps.getDefaultFactory().newImmutableMap(base, baseAscII);

    public static final HashByteCharMap ascIIToBaseMap = HashByteCharMaps.getDefaultFactory().newImmutableMap(baseAscII, base);

    public static byte getAscIIFromChar(char base){
        return baseToAscIIMap.get(base);
    }

    public static char getBaseFromAscII(byte baseAscII){
        return ascIIToBaseMap.get(baseAscII);
    }
}
