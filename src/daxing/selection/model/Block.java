package daxing.selection.model;

import daxing.common.SeqByte;
import org.apache.commons.lang3.StringUtils;

import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;

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

    public String[] getChrs() {
        return chrs;
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

    public int getSeqStart(int seqByteIndex){
        return getSeqStart()[seqByteIndex];
    }

    public int getSeqLen(int seqByteIndex){
        return getSeqLen()[seqByteIndex];
    }

    public boolean getIfMinus(int seqByteIndex){
        return getIfMinus()[seqByteIndex];
    }

    public SeqByte getSeqByte(int seqByteIndex){
        return getSeqByte()[seqByteIndex];
    }

    public int getBlockSize(){
        return getChrs().length;
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
     * fully alignmen region: no dash among all sequences
     * @return
     */
    public List<int[]> getFullyAlignmenRegion(){
        BitSet bitSet0, bitSet;
        bitSet0=this.getSeqByte(0).makeDashBeOne();
        for (int i = 1; i < this.getBlockSize(); i++) {
            bitSet = this.getSeqByte(i).makeDashBeOne();
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
