package daxing.informal;

import daxing.common.SeqByte;

import java.util.Arrays;
import java.util.Comparator;

public class MAFrecord  {
    private long id;
    private int score;
    private String[] taxons;
    private String[] chr;
    private int startPos[];
    private int seqLen[];
    private boolean[] ifMinus;
    private int[] chrLen;
    private SeqByte[] seq;

    public MAFrecord(long id, int score, String[] taxons, String[] chr, int[] startPos, int[] seqLen, boolean[] ifMinus,
              int[] chrLen, SeqByte[] seq){
        this.id=id;
        this.score=score;
        this.taxons=taxons;
        this.chr=chr;
        this.startPos=startPos;
        this.seqLen=seqLen;
        this.ifMinus=ifMinus;
        this.chrLen=chrLen;
        this.seq=seq;
    }

    public long getId(){
        return this.id;
    }

    public int getScore() {
        return this.score;
    }

    public String[] getTaxon() {
        return this.taxons;
    }

    public String getTaxon(int indexOfTaxon){
        return this.getTaxon()[indexOfTaxon];
    }

    public String[] getChr(){
        return this.chr;
    }

    public String getChr(int indexOfTaxon){
        return this.getChr()[indexOfTaxon];
    }

    public int[] getStartPos(){
        return this.startPos;
    }

    public int getStartPos(int indexOfTaxon){
        return this.getStartPos()[indexOfTaxon];
    }

    public int[] getSeqLen() {
        return seqLen;
    }

    public int getSeqLen(int indexOfTaxon){
        return this.getSeqLen()[indexOfTaxon];
    }

    public boolean[] getIfMinus() {
        return ifMinus;
    }

    public boolean getIfMinus(int indexOfTaxon){
        return this.getIfMinus()[indexOfTaxon];
    }

    public int[] getChrLen() {
        return chrLen;
    }

    public int getChrLen(int indexOfTaxon){
        return this.getChrLen()[indexOfTaxon];
    }

    public SeqByte[] getSeq() {
        return seq;
    }

    public SeqByte getSeq(int indexOfTaxon){
        return this.getSeq()[indexOfTaxon];
    }

    /**
     *  if it is minus, return [large, small], else return [small, large]
     * @return
     */
    public int[] startEnd(int indexOfTaxon){
        int[] startEnd=new int[2];
        int startPos=this.getStartPos(indexOfTaxon);
        int seqLen=this.getSeqLen(indexOfTaxon);
        int chrLen=this.getChrLen(indexOfTaxon);
        if(this.getIfMinus(indexOfTaxon)){
            startEnd[0]= chrLen-startPos;
            startEnd[1]= startEnd[0]-seqLen+1;
            return startEnd;
        }
        startEnd[0]=startPos+1;
        startEnd[1]=startEnd[0]+seqLen-1;
        return startEnd;
    }

    /**
     * 判断指定的染色体位置是否超过比对范围
     * @param indexOfTaxon
     * @param chr 指定的染色体
     * @param allelePosition 指定的染色体位置
     * @return
     */
    public boolean isOutOfRange(int indexOfTaxon, String chr, int allelePosition){
        if (this.getChr(indexOfTaxon).equals(chr)){
            int[] startEnd=this.startEnd(indexOfTaxon);
            Arrays.sort(startEnd);
            if (startEnd[0]<=allelePosition && startEnd[1]>=allelePosition){
                return false;
            }
            return true;
        }
        return true;
    }

    /**
     *
     * @param indexOfTaxon1
     * @param indexOfTaxon2
     * @param chrOfTaxon1
     * @param allelePositionOfTaxon1
     * @param refAltAlleleOfTaxon1
     * @return allele of Taxon2
     */
    public String getAllele(int indexOfTaxon1, int indexOfTaxon2, String chrOfTaxon1, int allelePositionOfTaxon1, String[] refAltAlleleOfTaxon1){
        if (this.isOutOfRange(indexOfTaxon1, chrOfTaxon1, allelePositionOfTaxon1)) return null;
        Comparator<String> comparator=Comparator.comparing(e->e.length());
        Arrays.sort(refAltAlleleOfTaxon1, comparator.reversed());
        SeqByte seqByte=this.getSeq(indexOfTaxon2);
        int[] startEnd=this.startEnd(indexOfTaxon1);
        int index=allelePositionOfTaxon1-startEnd[0];
        int range=refAltAlleleOfTaxon1[0].length();
        String seq=seqByte.getSequence(index, index+range-1);
        for (int i = 0; i < refAltAlleleOfTaxon1.length; i++) {
            if (seq.startsWith(refAltAlleleOfTaxon1[i])){
                return refAltAlleleOfTaxon1[i];
            }
        }
        return null;
    }


}
