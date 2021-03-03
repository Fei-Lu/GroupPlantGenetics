package daxing.load.complementary;

import daxing.common.ChrRange;
import daxing.common.PGF;
import daxing.common.WheatLineage;
import java.util.Arrays;
import java.util.List;

public class TriadsBlock {

    String triadsID;
    String[] geneName;
    int[] triadsCDSLen;
    List<String>[] blockGeneName;
    ChrRange[] chrRanges;

    public TriadsBlock(String triadsID, String[] geneName, int[] triadsCDSLen, List<String>[] blockGeneName){
        this.triadsID=triadsID;
        this.geneName=geneName;
        this.triadsCDSLen=triadsCDSLen;
        this.blockGeneName=blockGeneName;
        this.chrRanges=null;
    }

    public TriadsBlock(String triadsID, String[] geneName, int[] triadsCDSLen, List<String>[] blockGeneName,
                ChrRange[] chrRange){
        this.triadsID=triadsID;
        this.geneName=geneName;
        this.triadsCDSLen=triadsCDSLen;
        this.blockGeneName=blockGeneName;
        this.chrRanges=chrRange;
    }

    public TriadsBlock(ChrRange chrRange){
        this.triadsID=null;
        this.geneName=null;
        this.triadsCDSLen=null;
        this.blockGeneName=null;
        ChrRange[] chrRanges=new ChrRange[WheatLineage.values().length];
        String[] abd={"A","B","D"};
        int subIndex=Arrays.binarySearch(abd, chrRange.getChr().substring(1,2));
        for (int i = 0; i < chrRanges.length; i++) {
            if (subIndex==i){
                chrRanges[i]=chrRange;
            }else {
                chrRanges[i]=null;
            }
        }
        this.chrRanges=chrRanges;
    }

    public String getTriadsID() {
        return triadsID;
    }

    public String[] getGeneName() {
        return geneName;
    }

    public int[] getTriadsCDSLen() {
        return triadsCDSLen;
    }

    public List<String>[] getBlockGeneName() {
        return blockGeneName;
    }

    public ChrRange[] getChrRanges() {
        return chrRanges;
    }

    public boolean containsGene(String geneName){
        WheatLineage wheatLineage=WheatLineage.valueOf(geneName.substring(8,9));
        if (this.getBlockGeneName()[wheatLineage.getIndex()].contains(geneName)) return true;
        return false;
    }

    public String toString(){
        StringBuilder sb=new StringBuilder();
        String[] geneName=this.getGeneName();
        int[] cdsLen=this.getTriadsCDSLen();
        List<String>[] blockGeneName=this.getBlockGeneName();
        sb.append(this.getTriadsID()).append("\t");
        for (int i = 0; i < geneName.length; i++) {
            sb.append(geneName[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t");
        for (int i = 0; i < cdsLen.length; i++) {
            sb.append(cdsLen[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t");
        for (int i = 0; i < blockGeneName.length; i++) {
            sb.append(blockGeneName[i].size()).append(",");
        }
        sb.deleteCharAt(sb.length()-1);
        sb.append("\t");
        for (int i = 0; i < blockGeneName.length; i++) {
            for (int j = 0; j < blockGeneName[i].size(); j++) {
                sb.append(blockGeneName[i].get(j)).append(",");
            }
            sb.deleteCharAt(sb.length()-1);
            sb.append("\t");
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    public static ChrRange[] getChrRange(TriadsBlock triadsBlock, PGF pgf){
        pgf.sortGeneByName();
        ChrRange[] chrRanges=new ChrRange[WheatLineage.values().length];
        Arrays.fill(chrRanges, null);
        List<String>[] blockGeneName=triadsBlock.getBlockGeneName();
        int geneIndex;
        ChrRange tempChrRange, chrRange;
        String chr = null;
        for (int i = 0; i < blockGeneName.length; i++) {
            int startMini=Integer.MAX_VALUE, endMax=Integer.MIN_VALUE;
            for (int j = 0; j < blockGeneName[i].size(); j++) {
                geneIndex=pgf.getGeneIndex(blockGeneName[i].get(j));
                tempChrRange=ChrRange.changeToChrRange(pgf.getGene(geneIndex).getGeneRange());
                startMini=tempChrRange.getStart() < startMini ? tempChrRange.getStart() : startMini;
                endMax=tempChrRange.getEnd() > endMax ? tempChrRange.getEnd() : endMax;
                chr=tempChrRange.getChr();
            }
            chrRange=new ChrRange(chr, startMini, endMax);
            chrRanges[i]=chrRange;
        }
        return chrRanges;
    }

}
