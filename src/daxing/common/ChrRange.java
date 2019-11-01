package daxing.common;

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
}
