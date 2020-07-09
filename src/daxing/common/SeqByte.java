package daxing.common;

import com.google.common.primitives.Bytes;
import com.koloboke.collect.map.hash.HashByteByteMap;
import com.koloboke.collect.map.hash.HashByteByteMaps;
import com.koloboke.collect.map.hash.HashCharCharMap;
import com.koloboke.collect.map.hash.HashCharCharMaps;
import gnu.trove.list.array.TCharArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.DNAUtils;
import pgl.infra.dna.SequenceInterface;
import java.util.Arrays;
import java.util.BitSet;
import java.util.stream.IntStream;

/**
 * this class use to store a DNA sequence which includes "ATCGN-", it uses one byte to store a DNA base
 * @author Daxing Xu
 */

public class SeqByte implements SequenceInterface {

    byte[] seqByte = null;

    public SeqByte(String seq) {
        this.seqByte=seq.getBytes();
    }

    public SeqByte(byte[] seqByte){
        this.seqByte=seqByte;
    }

    public byte[] getSeqByte() {
        return seqByte;
    }

    /**
     * 返回所有子数组的索引
     * @param substr 子数组
     * @return
     */
    public TIntArrayList indexOfAll(byte[] substr){
        TIntArrayList indexes=new TIntArrayList();
        byte[] seqByte=this.seqByte;
        int index=0;
        byte[] des;
        int srcPos=0;
        int accumulator=0;
        int wordSize=0;
        while (index!=-1){
            des=new byte[seqByte.length-srcPos];
            System.arraycopy(seqByte, srcPos, des, 0, des.length);
            index=Bytes.indexOf(des, substr);
            if (index==-1) continue;
            accumulator=accumulator+index+wordSize;
            indexes.add(accumulator);
            wordSize=1;
            srcPos=accumulator+wordSize;
        }
        return indexes;
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
        byte[] baseByteWithN = DNAUtils.getBaseWithNAscIIArray();
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

    public char getReverseComplementaryBase(int index){
        char[] baseWithN_Dash={'A','T', 'C', 'G','-', 'N', 'a', 't', 'c', 'g', 'n'};
        char[] compleBaseWithN_Dash={'T', 'A', 'G', 'C', '-', 'N', 't', 'a', 'g', 'c', 'n'};
        HashCharCharMap baseCompleByteMap = HashCharCharMaps.getDefaultFactory().newImmutableMap(baseWithN_Dash, compleBaseWithN_Dash);
        return baseCompleByteMap.get(this.getBase(index));
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

    public int getDashCount(int startIndex, int endIndex){
        byte[] seqbyte=this.getSequence(startIndex, endIndex).getBytes();
        int countOfDash=0;
        for (int i = 0; i < seqbyte.length; i++) {
            if (seqbyte[i]==45){
                countOfDash++;
            }
        }
        return countOfDash;
    }

    /**
     *
     * @param position 0-based index of str without dash
     * @return index in str with dash
     */
    public int getIndexInDashStr(int position){
        String seqWithDash=this.getSequence();
        TIntArrayList indexOfDash= StringTool.getIndexOfSubStr(seqWithDash, "-");
        TIntArrayList allIndex=new TIntArrayList(IntStream.range(0, seqWithDash.length()).toArray());
        allIndex.removeAll(indexOfDash);
        return allIndex.get(position);
    }

    /**
     *
     * @param indexInDashSeq
     * @return index in Seq Without Dash
     */
    public int getIndexInSeqWithoutDash(int indexInDashSeq){
        String seqWithDash=this.getSequence();
        TIntArrayList indexOfDash= StringTool.getIndexOfSubStr(seqWithDash, "-");
        TIntArrayList allIndex=new TIntArrayList(IntStream.range(0, seqWithDash.length()).toArray());
        allIndex.removeAll(indexOfDash);
        return allIndex.binarySearch(indexInDashSeq);
    }

    /**
     *
     * @return
     */
    public BitSet makeDashBeZero(){
        byte[] dashACGTN={45, 65, 67, 71, 78, 84};
        byte[] dashBeZero={0, 1, 1, 1, 0, 1};
        HashByteByteMap dashACGTNToZeroMap = HashByteByteMaps.getDefaultFactory().newImmutableMap(dashACGTN, dashBeZero);
        BitSet bitSet = new BitSet();
        byte value;
        boolean res;
        for (int i = 0; i < seqByte.length; i++) {
            value=dashACGTNToZeroMap.get(seqByte[i]);
            res= value == 1 ? true : false;
            bitSet.set(i, res);
        }
        return bitSet;
    }

    @Override
    public SequenceInterface getSequenceInterface(int startIndex, int endIndex) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

}
