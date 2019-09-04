package daxing.filterSNP;

import format.position.ChrPos;

public class Dot extends ChrPos {

    private double depth;
    private double sd;

    public Dot(short chr, int pos, double depth, double sd){
        super(chr, pos);
        this.depth=depth;
        this.sd = sd;
    }

    public double getDepth(){
        return depth;
    }

    public double getSd() {
        return sd;
    }
}
