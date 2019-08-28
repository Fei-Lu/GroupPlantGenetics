package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;

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
     * @return 1-based coordinates
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
     * 根据给定Taxon的index和ChrPos，返回对应MAFrecord比对中另一个Taxon的allele
     * @param indexOfTaxon
     * @param chrPosOfTaxon 1-based coordinates in vcf
     * @param chrConvertionRule
     * @return allele of another taxon
     */
    public String getAllele(int indexOfTaxon, ChrPos chrPosOfTaxon, ChrConvertionRule chrConvertionRule){
        int allelePositionOfTaxon=chrConvertionRule.getRefPositionFromVCF(chrPosOfTaxon);
        int[] startEnd1=this.startEnd(indexOfTaxon);
        TIntArrayList indexForTaxons=new TIntArrayList();
        indexForTaxons.add(0);
        indexForTaxons.add(1);
        indexForTaxons.remove(indexOfTaxon);
        SeqByte seqByte2=this.getSeq(indexForTaxons.get(0));
        int index=allelePositionOfTaxon-startEnd1[0];
        char base2=seqByte2.getBase(Math.abs(index));
        if (base2=='-'){
            return null;
        }
        return String.valueOf(base2);
    }

}
