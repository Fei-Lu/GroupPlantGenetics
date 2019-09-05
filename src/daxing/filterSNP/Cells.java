package daxing.filterSNP;

import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.stream.DoubleStream;

public class Cells {

    private List<List<Cell>> cellList;
    private double[] boundaryOfDepth;
    private double depthWindow;
    private double[] boundaryOfSD;
    private double sdWindow;
    private int countOfCell;

    public Cells(double maxOfdepth, double maxOfSD, int numberOfCell){
        this.countOfCell=numberOfCell*numberOfCell;
        this.initializeCells(maxOfdepth, maxOfSD, numberOfCell);
    }

    private void initializeCells(double maxOfdepth, double maxOfSD, long numberOfCell){
        List<List<Cell>> cellss=new ArrayList<>();
        List<Cell> cells;
        double depthWindow=maxOfdepth/numberOfCell;
        double[] boundaryOfDepth=DoubleStream.iterate(0, n->n+depthWindow).limit(numberOfCell).toArray();
        double sdWindow=maxOfSD/numberOfCell;
        double[] boundaryOfSD=DoubleStream.iterate(0, n->n+sdWindow).limit(numberOfCell).toArray();
        for (int i = 0; i < boundaryOfDepth.length; i++) {
            cells=new ArrayList<>();
            for (int j = 0; j < boundaryOfSD.length; j++) {
                cells.add(new Cell(boundaryOfDepth[i], depthWindow, boundaryOfSD[j], sdWindow));
            }
            cellss.add(cells);
        }
        this.cellList=cellss;
        this.boundaryOfDepth=boundaryOfDepth;
        this.depthWindow=depthWindow;
        this.boundaryOfSD=boundaryOfSD;
        this.sdWindow=sdWindow;
    }

    public double[] getBoundaryOfDepth(){
        return boundaryOfDepth;
    }

    public double getDepthWindow() {
        return depthWindow;
    }

    public double[] getBoundaryOfSD() {
        return boundaryOfSD;
    }

    public double getSdWindow() {
        return sdWindow;
    }

    public List<List<Cell>> getCellList() {
        return cellList;
    }

    public Cell getCell(int indexOfDepth, int indexOfSD){
        return this.getCellList().get(indexOfDepth).get(indexOfSD);
    }

    public Cell getCell(int[] indexOfDepthSD){
        return this.getCell(indexOfDepthSD[0], indexOfDepthSD[1]);
    }

    public int getCountOfCell() {
        return countOfCell;
    }

    public int[] binarySearch(Dot dot){
        double[] boundaryOfDepth=this.getBoundaryOfDepth();
        double depthWindow=this.getDepthWindow();
        double[] boundaryOfSD=this.boundaryOfSD;
        double sdWindow=this.sdWindow;
        int indexInDepth= Arrays.binarySearch(boundaryOfDepth, dot.getDepth());
        int indexInSD=Arrays.binarySearch(boundaryOfSD, dot.getSd());
        if (indexInDepth<0){
            indexInDepth=-indexInDepth-2;
        }
        if (indexInSD<0){
            indexInSD=-indexInSD-2;
        }
        int[] indexInDepthSD=new int[2];
        indexInDepthSD[0]=indexInDepth;
        indexInDepthSD[1]=indexInSD;
        return indexInDepthSD;
    }

    public double getDotNumOfAllCell(){
        Cell cell;
        double size=0D;
        for (int i = 0; i < this.getCellList().size(); i++) {
            for (int j = 0; j < this.getCellList().get(0).size(); j++) {
                cell=this.getCell(i, j);
                size=size+cell.size();
            }
        }
        return size;
    }

    private List<Cell> getAllCell(){
        List<Cell> cells=new ArrayList<>();
        for (int i = 0; i < this.getCellList().size(); i++) {
            for (int j = 0; j < this.getCellList().get(0).size(); j++) {
                cells.add(this.getCell(i,j));
            }
        }
        return cells;
    }

    public void write(String ouputDir){
        File cellSizeFile=new File(ouputDir, "cellSize.txt");
        File positionSubDir=new File(ouputDir, "position");
        positionSubDir.mkdir();
        BufferedWriter bwCellSizeFile= IOUtils.getTextWriter(cellSizeFile.getAbsolutePath());
        BufferedWriter[] bwForDotPosition=new BufferedWriter[this.countOfCell];
        for (int i = 0; i < bwForDotPosition.length; i++) {
            bwForDotPosition[i]=IOUtils.getTextWriter(new File(positionSubDir, "cell"+ PStringUtils.getNDigitNumber(5, i) +"Position.txt").getAbsolutePath());
        }
        try{
            bwCellSizeFile.write("Cell"+"\t"+"Size"+"\t"+"cumulativePercentage"+"\n");
            for (int i = 0; i < bwForDotPosition.length; i++) {
                bwForDotPosition[i].write("Chr"+"\t"+"POS"+"\t"+"Depth"+"\t"+"SD");
                bwForDotPosition[i].newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        try{
            List<Cell> cells=this.getAllCell();
            Comparator<Cell> cellComparator=Comparator.comparing(c->c.size());
            cellComparator=cellComparator.reversed();
            Collections.sort(cells, cellComparator);
            Cell cell;
            int count=0;
            double rate=0;
            short chr;
            int pos;
            double depth;
            double sd;
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < cells.size(); i++) {
                cell=cells.get(i);
                rate=rate+cell.size()/this.getDotNumOfAllCell();
                bwCellSizeFile.write(count+"\t"+cell.size()+"\t"+rate);
                bwCellSizeFile.newLine();
                for (int k = 0; k < cell.size(); k++) {
                    sb=new StringBuilder();
                    chr=cell.getDotList().get(k).getChromosome();
                    pos=cell.getDotList().get(k).getPosition();
                    depth=cell.getDotList().get(k).getDepth();
                    sd=cell.getDotList().get(k).getSd();
                    sb.append(chr).append("\t").append(pos).append("\t").append(depth)
                            .append("\t").append(sd);
                    bwForDotPosition[count].write(sb.toString());
                    bwForDotPosition[count].newLine();
                }
                count++;
            }
            bwCellSizeFile.flush();
            bwCellSizeFile.close();
            for (int i = 0; i < bwForDotPosition.length; i++) {
                bwForDotPosition[i].flush();
                bwForDotPosition[i].close();
            }
            System.out.println("Total sites in this file is "+ (int)this.getDotNumOfAllCell());
        }catch (Exception e){
            e.printStackTrace();
        }
        this.writeForGraph(new File(ouputDir, "cellGraph.txt").getAbsolutePath());
    }

    public void writeForGraph(String outputfile){
        try(BufferedWriter bw=IOUtils.getTextWriter(outputfile)){
            bw.write("Depth"+"\t"+"SD"+"\t"+"Density"+"\n");
            Cell cell;
            double cellDepth;
            double cellSD;
            double cellSizeDensity;
            StringBuilder sb;
            for (int i = 0; i < this.getCellList().size(); i++) {
                for (int j = 0; j < this.getCellList().get(0).size(); j++) {
                    cell=this.getCellList().get(i).get(j);
                    cellDepth=cell.getMeanOfDepth();
                    cellSD=cell.getMeanOfSD();
                    cellSizeDensity=(cell.size())/(this.getDotNumOfAllCell());
                    sb=new StringBuilder();
                    sb.append(cellDepth).append("\t").append(cellSD).append("\t").append(cellSizeDensity);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}
