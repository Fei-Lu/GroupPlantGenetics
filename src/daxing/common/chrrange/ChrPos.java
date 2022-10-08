package daxing.common.chrrange;

public class ChrPos implements Comparable<ChrPos>{

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

    @Override
    public int compareTo(ChrPos o) {
        if (this.chr.equals(o.chr)){
            if (this.pos < o.pos){
                return -1;
            }else if (this.pos==o.pos){
                return 0;
            }else {
                return 1;
            }
        }else {
            return this.chr.compareTo(o.chr);
        }
    }
}
