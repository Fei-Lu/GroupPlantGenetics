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

    public void write(String outPutDir){
        File cellSizeFile=new File(outPutDir, "cellSize.txt");
        File positionSubDir=new File(outPutDir, "position");
        positionSubDir.mkdir();
        BufferedWriter bwCellSizeFile= IOUtils.getTextWriter(cellSizeFile.getAbsolutePath());
        BufferedWriter[] bwForDotPosition=new BufferedWriter[this.countOfCell];
        for (int i = 0; i < bwForDotPosition.length; i++) {
            bwForDotPosition[i]=IOUtils.getTextWriter(new File(positionSubDir, "cell"+ PStringUtils.getNDigitNumber(3, i) +"Position.txt").getAbsolutePath());
        }
        try{
            bwCellSizeFile.write("Cell"+"\t"+"Size"+"\n");
            for (int i = 0; i < bwForDotPosition.length; i++) {
                bwForDotPosition[i].write("POS");
                bwForDotPosition[i].newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        try{
            Cell cell;
            int count=0;
            for (int i = 0; i < this.cellList.size(); i++) {
                for (int j = 0; j < cellList.get(0).size(); j++) {
                    cell=this.getCell(i,j);
                    bwCellSizeFile.write(count+"\t"+cell.size());
                    bwCellSizeFile.newLine();
                    for (int k = 0; k < cell.size(); k++) {
                        bwForDotPosition[count].write(String.valueOf(cell.getDotList().get(k).getPosition()));
                        bwForDotPosition[count].newLine();
                    }
                    count++;
                }
            }
            bwCellSizeFile.flush();
            bwCellSizeFile.close();
            for (int i = 0; i < bwForDotPosition.length; i++) {
                bwForDotPosition[i].flush();
                bwForDotPosition[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        this.writeForGraph(new File(outPutDir, "cellGraph.txt").getAbsolutePath());

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
