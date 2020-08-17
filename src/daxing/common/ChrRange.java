package daxing.common;

import pgl.infra.utils.wheat.RefV1Utils;

/**
 * @author Daxing Xu
 */
public class ChrRange {

    String chr;
    int start;
    int end;

    public ChrRange(String chr, int start, int end){
        this.chr=chr;
        this.start=start;
        this.end=end;
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getLength(){
        return end-start;
    }

    public int getChrID(){
        return RefV1Utils.getChrID(chr, start);
    }

    public int getVCFStart(){
        return RefV1Utils.getPosOnChrID(chr, start);
    }

    public int getVCFEnd(){
        return RefV1Utils.getPosOnChrID(chr, end);
    }

    public short getVCFStartChrID(){
        return (short) RefV1Utils.getChrID(chr, start);
    }

    public short getVCFEndChrID(){
        return (short) RefV1Utils.getChrID(chr, end);
    }
}
