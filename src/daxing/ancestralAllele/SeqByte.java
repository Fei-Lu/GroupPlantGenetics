package daxing.ancestralAllele;

import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import format.dna.DNAUtils;
import format.dna.SequenceInterface;
import gnu.trove.list.array.TCharArrayList;

import java.util.Arrays;

/**
 * this class use to store a DNA sequence of an alignment in MAF file format
 * it includes ATCGN-  and corresponding reverse complementary sequence (eg: -NCGAT)
 */

public class SeqByte implements SequenceInterface {

    byte[] seqByte = null;

    public SeqByte(String seq) {
        this.seqByte=seq.getBytes();
    }

    @Override
    public int getSequenceLength() {
        return seqByte.length;
    }

    @Override
    public double getProportionA() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 65) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionT() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionG() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 71) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getProportionC() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public double getGCContent() {
        int cnt = 0;
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 67 || this.seqByte[i] == 84) cnt++;
        }
        return (double)cnt/this.seqByte.length;
    }

    @Override
    public char getBase(int positionIndex) {
        return (char)this.seqByte[positionIndex];
    }

    @Override
    public String getSequence() {
        return new String(this.seqByte);
    }

    @Override
    public String getSequence(int startIndex, int endIndex) {
        return new String(this.seqByte, startIndex, endIndex-startIndex);
    }

    @Override
    public String getReverseComplementarySeq() {
        return this.getReverseComplementarySeq(0, this.getSequenceLength());
    }

    /**
     * ATCGN-   -NCGAT
     * @param startIndex inclusive
     * @param endIndex exclusive
     * @return reverse complementary sequence included 'N' and '-'
     */
    @Override
    public String getReverseComplementarySeq(int startIndex, int endIndex) {
        Byte[] baseByteWithN_Dash={45, 65, 67, 71, 78, 84};
        Byte[] compleBaseByteWithN_Dash={45, 84, 71, 67, 78, 65};
        HashByteByteMap baseCompleByteMap = HashByteByteMaps.getDefaultFactory().newImmutableMap(baseByteWithN_Dash, compleBaseByteWithN_Dash);
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = baseCompleByteMap.get(seqByte[endIndex-i-1]);
        }
        return new String(reverseByte);
    }

    @Override
    public boolean isThereNonACGTNBase() {
        for (int i = 0; i < this.getSequenceLength(); i++) {
            if (this.seqByte[i] == 78) return true;
        }
        return false;
    }

    @Override
    public boolean isThereN() {
        byte[] baseByteWithN = DNAUtils.getBaseWithNByteArray();
        for (int i = 0; i < this.getSequenceLength(); i++) {
            int index = Arrays.binarySearch(baseByteWithN, this.getBaseByte(i));
            if (index < 0) {
                System.out.println(this.getBase(i));
                return true;
            }
        }
        return false;
    }

    public byte getBaseByte (int positionIndex) {
        return this.seqByte[positionIndex];
    }

    /**
     *
     * @param startIndex
     * @param endIndex
     * @return
     */
    public String getReverseComplementarySeqWithoutDash(int startIndex, int endIndex){
        Byte[] baseByteWithN_Dash={65, 67, 71, 78, 84};
        Byte[] compleBaseByteWithN_Dash={84, 71, 67, 78, 65};
        HashByteByteMap baseCompleByteMap = HashByteByteMaps.getDefaultFactory().newImmutableMap(baseByteWithN_Dash, compleBaseByteWithN_Dash);;
        byte[] reverseByte = new byte[endIndex - startIndex];
        for (int i = 0; i < reverseByte.length; i++) {
            reverseByte[i] = baseCompleByteMap.get(seqByte[endIndex-i-1]);
        }
        return new String(reverseByte);
    }

    public String getReverseComplementarySeqWithoutDash() {
        return this.getReverseComplementarySeqWithoutDash(0, this.getSequenceLength());
    }

    public int getSequenceLengthWithoutDash() {
        return this.getSequenceWithoutDash().length();
    }

    public String getSequenceWithoutDash() {
        TCharArrayList seqByte=new TCharArrayList();
        for (int i = 0; i < this.getSequenceLength() ; i++) {
            if (this.getBase(i)=='-'){
                continue;
            }
            seqByte.add(this.getBase(i));
        }
        return new String(seqByte.toArray());
    }
}
