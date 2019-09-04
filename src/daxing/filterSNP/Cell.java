package daxing.filterSNP;

import java.util.ArrayList;
import java.util.List;

public class Cell {

    private double depthBoundary;
    private double depthWindow;
    private double sdBoundary;
    private double sdWindow;
    private List<Dot> dotList;

    public Cell(double depthBoundary, double depthWindow, double sdBoundary, double sdWindow){
        this.depthBoundary=depthBoundary;
        this.depthWindow=depthWindow;
        this.sdBoundary=sdBoundary;
        this.sdWindow=sdWindow;
        this.dotList=new ArrayList<>();
    }

    public void add(Dot dot){
        this.dotList.add(dot);
    }

    public double getDepthBoundary() {
        return depthBoundary;
    }

    public double getDepthWindow() {
        return depthWindow;
    }

    public double getSdBoundary() {
        return sdBoundary;
    }

    public double getSdWindow() {
        return sdWindow;
    }

    public List<Dot> getDotList() {
        return dotList;
    }

    public int size(){
        return this.getDotList().size();
    }
}
