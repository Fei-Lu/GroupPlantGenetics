package daxing.selection.model;

import daxing.common.ChrPos;
import daxing.common.SeqByte;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;

/**
 * an alignment record without indel (block) in a maf file
 */
public class Block {

    String[] chrs;
    int[] seqStart;
    int[] seqLen;
    boolean[] ifMinus;
    SeqByte[] seqByte;

    public Block(List<String> blocks){
        chrs=new String[blocks.size()];
        seqStart=new int[blocks.size()];
        seqLen=new int[blocks.size()];
        ifMinus=new boolean[blocks.size()];
        seqByte=new SeqByte[blocks.size()];
        String[] temp;
        String subgenome;
        int chr;
        for (int i = 0; i < blocks.size(); i++) {
            temp= StringUtils.split(blocks.get(i), " ");
            subgenome=temp[1].substring(5,6);
            chr=Integer.parseInt(temp[1].substring(10,11));
            chrs[i]=chr+subgenome;
            seqStart[i]=Integer.parseInt(temp[2]);
            seqLen[i]=Integer.parseInt(temp[3]);
            ifMinus[i]=temp[4].equals("-") ? true : false;
            seqByte[i]=new SeqByte(temp[6].toUpperCase());
        }
    }

    public Block(String[] chrs, int[] seqStart, int[] seqLen, boolean[] ifMinus, SeqByte[] seqByte){
        this.chrs=chrs;
        this.seqStart=seqStart;
        this.seqLen=seqLen;
        this.ifMinus=ifMinus;
        this.seqByte=seqByte;
    }

    public Block(String chr, int pos){
        String[] chrs=new String[2];
        chrs[0]=chr;
        chrs[1]="1E";
        int[] seqStart=new int[2];
        seqStart[0]=pos-1;
        seqStart[1]=-1;
        int[] seqLen=new int[2];
        seqLen[0]=1;
        seqLen[1]=-1;
        boolean[] ifMinus=new boolean[2];
        ifMinus[0]=false;
        ifMinus[1]=false;
        SeqByte[] seqBytes=new SeqByte[2];
        this.chrs=chrs;
        this.seqStart=seqStart;
        this.seqLen=seqLen;
        this.ifMinus=ifMinus;
        this.seqByte=seqBytes;
    }

    public String[] getChrs() {
        return chrs;
    }

    public String[] getSubgenome(){
        String[] chrs=this.getChrs();
        return Arrays.stream(chrs).map(s->s.substring(1,2)).toArray(String[]::new);
    }

    public int[] getSeqStart() {
        return seqStart;
    }

    public int[] getSeqLen() {
        return seqLen;
    }

    public boolean[] getIfMinus() {
        return ifMinus;
    }

    public SeqByte[] getSeqByte() {
        return seqByte;
    }

    public String getChr(int seqByteIndex){
        return getChrs()[seqByteIndex];
    }

    public int getSeqStart(int index){
        return getSeqStart()[index];
    }

    public int getSeqLen(int index){
        return getSeqLen()[index];
    }

    public boolean getIfMinus(int index){
        return getIfMinus()[index];
    }

    public SeqByte getSeqByte(int index){
        return getSeqByte()[index];
    }

    public String getChr(String subgenome){
        int index=this.getSubIndex(subgenome);
        return this.getChr(index);
    }

    public int getSeqStart(String subgenome){
        int index=this.getSubIndex(subgenome);
        return this.getSeqStart(index);
    }

    public int getSeqLen(String subgenome){
        int index=this.getSubIndex(subgenome);
        return this.getSeqLen(index);
    }

    public boolean getIfMinus(String subgenome){
        int index=this.getSubIndex(subgenome);
        return this.getIfMinus(index);
    }

    public SeqByte getSeqByte(String subgenome){
        int index=this.getSubIndex(subgenome);
        return this.getSeqByte(index);
    }

    public int getBlockSize(){
        return getChrs().length;
    }

    public int getSubIndex(String subgenome){
        String[] subs=this.getSubgenome();
        return Arrays.binarySearch(subs, subgenome);
    }

    public boolean haveAlignment(String subgenome){
        int index=getSubIndex(subgenome);
        if (index < 0) return false;
        return true;
    }

    /**
     *
     * @param posA one-based position in subgenome A
     * @return PosB one-based position in subgenome B
     */
    public ChrPos getChrPosB(int posA){
        int index=posA-1-this.getSeqStart("A");
        int posBZeroBased=this.getSeqStart("B")+index;
        String chrB=this.getChr("B");
        int sizeB= RefV1Utils.getChromosomeLength(chrB);
        boolean ifBMinus=this.getIfMinus("B");
        int posBOneBased = ifBMinus ? sizeB - posBZeroBased : posBZeroBased + 1;
        return new ChrPos(chrB, posBOneBased);
    }

    /**
     *
     * @param posA one-based position in subgenome A
     * @return PosD one-based position in subgenome D
     */
    public ChrPos getChrPosD(int posA){
        int index=posA-1-this.getSeqStart("A");
        int posDZeroBased=this.getSeqStart("D")+index;
        String chrD=this.getChr("D");
        int sizeD= RefV1Utils.getChromosomeLength(chrD);
        boolean ifDMinus=this.getIfMinus("D");
        int posDOneBased = ifDMinus ? sizeD - posDZeroBased : posDZeroBased + 1;
        return new ChrPos(chrD, posDOneBased);
    }

    public List<Block> subsetBlockWithoutIndel(){
        List<Block> blocks=new ArrayList<>();
        List<int[]> rangeList=this.getFullyAlignmenRegion();
        int start, end, startIndexInSeqWithoutDash;
        String[] chrs;
        boolean[] ifMinus;
        int[] seqStart;
        int[] seqLen;
        SeqByte[] seqByte;
        for (int i = 0; i < rangeList.size(); i++) {
            start=rangeList.get(i)[0];
            end=rangeList.get(i)[1];
            chrs=this.getChrs();
            ifMinus=this.getIfMinus();
            seqStart=new int[this.getBlockSize()];
            seqLen=new int[this.getBlockSize()];
            seqByte=new SeqByte[this.getBlockSize()];
            for (int j = 0; j < this.getBlockSize(); j++) {
                startIndexInSeqWithoutDash=this.getSeqByte(j).getIndexInSeqWithoutDash(start);
                seqStart[j]=startIndexInSeqWithoutDash+this.getSeqStart(j);
                seqLen[j]=end-start;
                seqByte[j]=new SeqByte(this.getSeqByte(j).getSequence(start, end));
            }
            blocks.add(new Block(chrs, seqStart, seqLen, ifMinus, seqByte));
        }
        return blocks;
    }

    /**
     * fully alignment region: no dash among all sequences
     * @return
     */
    public List<int[]> getFullyAlignmenRegion(){
        BitSet bitSet0, bitSet;
        bitSet0=this.getSeqByte(0).makeDashBeZero();
        for (int i = 1; i < this.getBlockSize(); i++) {
            bitSet = this.getSeqByte(i).makeDashBeZero();
            bitSet0.and(bitSet);
        }
        List<int[]> rangeList=new ArrayList<>();
        int[] range;
        for (int i = bitSet0.nextSetBit(0); i > -1; i=bitSet0.nextSetBit(i+1)){
            range=new int[2];
            range[0]=i;
            i=bitSet0.nextClearBit(i+1);
            range[1]=i;
            rangeList.add(range);
        }
        return rangeList;
    }

    public String toString(){
        StringBuilder sb=new StringBuilder();
        String ifMinusstr = null;
        for (int i = 0; i < this.getBlockSize(); i++) {
            sb.append(this.chrs[i]).append(" ");
            sb.append(getNDigitNumber(9, this.seqStart[i])).append(" ").append(seqLen[i]).append(" ");
            ifMinusstr=ifMinus[i] ? "-" : "+";
            sb.append(ifMinusstr).append(" ").append(seqByte[i].getSequence()).append("\n");
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    /**
     * Return a string of number filled with 0 on the left
     * @param n
     * @param num
     * @return
     */
    public static String getNDigitNumber (int n, int num) {
        String s = String.valueOf(num);
        int cnt = n-s.length();
        for (int i = 0; i < cnt; i++) {
            s = " "+s;
        }
        return s;
    }
}
