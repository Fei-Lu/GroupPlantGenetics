package daxing.load.complementary;

import daxing.common.WheatLineage;

import java.util.List;

public class TriadsBlock {

    String triadsID;
    String[] geneName;
    int[] triadsCDSLen;
    List<String>[] blockGeneName;

    TriadsBlock(String triadsID, String[] geneName, int[] triadsCDSLen, List<String>[] blockGeneName){
        this.triadsID=triadsID;
        this.geneName=geneName;
        this.triadsCDSLen=triadsCDSLen;
        this.blockGeneName=blockGeneName;
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

    public void makeBlock(List<String>[] blockGeneName){
        this.blockGeneName=blockGeneName;
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
}
