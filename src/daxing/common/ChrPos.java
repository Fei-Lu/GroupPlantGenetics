package daxing.common;

public class ChrPos {

    String chr;
    int pos;

    public ChrPos(String chr, int pos){
        this.chr=chr;
        this.pos=pos;
    }

    public String getChr() {
        return chr;
    }

    public int getPos() {
        return pos;
    }

    @Override
    public String toString() {
        return "ChrPos{" +
                "chr='" + chr + '\'' +
                ", pos=" + pos +
                '}';
    }
}
