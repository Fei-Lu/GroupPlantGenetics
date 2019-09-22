package daxing.ancestralAllele;

import daxing.common.ChrConvertionRule;
import daxing.common.StringTool;
import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.stream.IntStream;

public class MAFrecord  {
    private long id;
    private long score;
    private String[] taxons;
    private String[] chr;
    private int startPos[];
    private int seqLen[];
    private boolean[] ifMinus;
    private int[] chrLen;
    private SeqByte[] seq;

    public MAFrecord(long id, long score, String[] taxons, String[] chr, int[] startPos, int[] seqLen, boolean[] ifMinus,
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

    public MAFrecord(String[] chr, int[] startPos, int[] seqLen, boolean[] ifMinus,
                     int[] chrLen, SeqByte[] seq){
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

    public long getScore() {
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

    public int getStartPos1_based(int indexOfTaxon){
        int[] startEnd=this.startEnd(indexOfTaxon);
        return Arrays.stream(startEnd).min().getAsInt();
    }

    /**
     * 根据给定Taxon的index和ChrPos，返回对应MAFrecord比对中另一个Taxon的allele
     * @param indexOfTaxon 0 or 1
     * @param chrPosOfTaxon 1-based coordinates in vcf
     * @param chrConvertionRule
     * @return allele of another taxon
     */
    public String getOutgroupAllele(int indexOfTaxon, ChrPos chrPosOfTaxon, ChrConvertionRule chrConvertionRule){
        int allelePosition1_Based=chrConvertionRule.getRefPosFromVCFChrPos(chrPosOfTaxon);
        int[] startEnd=this.startEnd(indexOfTaxon);
        TIntArrayList indexForTaxons=new TIntArrayList();
        indexForTaxons.add(0);
        indexForTaxons.add(1);
        indexForTaxons.remove(indexOfTaxon);
        SeqByte seqByte=this.getSeq(indexOfTaxon);
        SeqByte seqByteOfOutGroup=this.getSeq(indexForTaxons.get(0));
        int index=Math.abs(allelePosition1_Based-startEnd[0]);
        int indexInDashStr=seqByte.getIndexInDashStr(index);
        char base=seqByteOfOutGroup.getBase(indexInDashStr);
        return String.valueOf(base);
    }

    /**
     * 三个参数指同一个taxon的index, ChrPos, ChrConvertionRule
     * @param indexOfTaxon wheat
     * @param chrPosOfTaxon wheat
     * @param chrConvertionRuleOfTaxon wheat
     * @return allele of this taxon
     */
    public String getRefAllele(int indexOfTaxon, ChrPos chrPosOfTaxon, ChrConvertionRule chrConvertionRuleOfTaxon){
        int allelePosition1_Based=chrConvertionRuleOfTaxon.getRefPosFromVCFChrPos(chrPosOfTaxon);
        int[] startEnd=this.startEnd(indexOfTaxon);
        SeqByte seqByte=this.getSeq(indexOfTaxon);
        int index=Math.abs(allelePosition1_Based-startEnd[0]);
        int indexInDashStr=seqByte.getIndexInDashStr(index);
        char base=seqByte.getBase(indexInDashStr);
        return String.valueOf(base);
    }

    /**
     *
     * @param indexOfTaxon wheat
     * @return
     */
    public Map<ChrPos, String[]> getRefCoordinateOfOutGroup(int indexOfTaxon, ChrConvertionRule chrConvertionRule){
        TIntArrayList indexForTaxons=new TIntArrayList();
        indexForTaxons.add(0);
        indexForTaxons.add(1);
        indexForTaxons.remove(indexOfTaxon);
        SeqByte seqByte1=this.getSeq(indexOfTaxon);
        SeqByte seqByte2=this.getSeq(indexForTaxons.get(0));
        String seqWithDash1=seqByte1.getSequence();
        String seqWithDash2=seqByte2.getSequence();
        int[] indexOfDash1= StringTool.getIndexOfSubStr(seqWithDash1, "-");
        int[] indexOfDash2= StringTool.getIndexOfSubStr(seqWithDash2, "-");
        TIntHashSet allIndex1=new TIntHashSet(IntStream.range(0, seqWithDash1.length()).toArray());
        TIntHashSet allIndex2=new TIntHashSet(IntStream.range(0, seqWithDash2.length()).toArray());
        allIndex1.removeAll(indexOfDash1);
        allIndex2.removeAll(indexOfDash2);
        allIndex1.retainAll(allIndex2);
        TIntArrayList alleleIndex=new TIntArrayList(allIndex1);
        alleleIndex.sort();
        Map<ChrPos, String[]> outgroupAllele=new HashMap<>();
        ChrPos chrPos;
        short chr;
        int pos;
        String[] refOutgroupAllele;
        int[] startEnd=this.startEnd(indexOfTaxon);
        int[] seqCoordinate1_based;
        if (this.getIfMinus()[indexOfTaxon]){
            seqCoordinate1_based=IntStream.iterate(startEnd[0], n->n-1).limit(startEnd[0]-startEnd[startEnd.length-1]+1).toArray();
        }else {
            seqCoordinate1_based=IntStream.iterate(startEnd[0], n->n+1).limit(startEnd[startEnd.length-1]-startEnd[0]+1).toArray();
        }
        for (int i = 0; i < alleleIndex.size(); i++) {
            pos=seqCoordinate1_based[seqByte1.getIndexInSeqWithoutDash(alleleIndex.get(i))];
            chr=(short) chrConvertionRule.getVCFChrFromRefChrPos(this.getChr(indexOfTaxon), pos);
            chrPos=new ChrPos(chr, chrConvertionRule.getVCFPosFromRefChrPos(this.getChr(indexOfTaxon), pos));
            refOutgroupAllele=new String[2];
            if (this.getIfMinus(indexOfTaxon)) {
                refOutgroupAllele[0]=String.valueOf(seqByte1.getReverseComplementaryBase(alleleIndex.get(i)));
                refOutgroupAllele[1]=String.valueOf(seqByte2.getReverseComplementaryBase(alleleIndex.get(i)));
            }else {
                refOutgroupAllele[0]=String.valueOf(seqByte1.getBase(alleleIndex.get(i)));
                refOutgroupAllele[1]=String.valueOf(seqByte2.getBase(alleleIndex.get(i)));
            }
            if (refOutgroupAllele[0].equals("N")) continue;
            if (refOutgroupAllele[0].equals("n")) continue;
            if (refOutgroupAllele[1].equals("N")) continue;
            if (refOutgroupAllele[1].equals("n")) continue;
            outgroupAllele.put(chrPos, refOutgroupAllele);
        }
        return outgroupAllele;
    }

}
