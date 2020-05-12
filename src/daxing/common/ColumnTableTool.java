package daxing.common;

import pgl.infra.table.ColumnTable;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ColumnTableTool<T> extends ColumnTable<T> {

    public ColumnTableTool(String infileS) {
        super(infileS);
    }

    public ColumnTableTool(List<String> header, List<List<T>> cells){
        super(header,cells);
    }

    public static <T> ColumnTableTool<T> newInstance(List<String> header, List<T> column){
        List<List<T>> cells=new ArrayList<>();
        cells.add(column);
        return new ColumnTableTool<>(header, cells);
    }

    public RowTableTool<T> transformToRowTable(){
        List<List<T>> rowCells=new ArrayList<>();
        List<T> row;
        for (int i = 0; i < this.getRowNumber(); i++) {
            row=new ArrayList<>();
            for (int j = 0; j < this.getColumnNumber(); j++) {
                row.add(this.getCell(i, j));
            }
            rowCells.add(row);
        }
        return new RowTableTool<>(this.header, rowCells);
    }

    public List<List<T>> getCells(){
        return this.cells;
    }

    public double[] getColumnAsDoubleArray(String columnName) {
        int columnIndex=this.hiMap.get(columnName);
        List<T> column=this.getColumn(columnIndex);
        double[] values=new double[column.size()];
        T ob;
        for (int i = 0; i < column.size(); i++) {
            ob=column.get(i);
            values[i]=Double.parseDouble((String)ob);
        }
        return values;
    }

    /**
     * 选择列组成ColumnTable
     */
    public ColumnTableTool<T> selectColumn(int[] columnIndex){
        Arrays.sort(columnIndex);
        List<String> header=new ArrayList<>();
        List<List<T>> data=new ArrayList<>();
        for (int i = 0; i < columnIndex.length; i++) {
            data.add(this.getColumn(columnIndex[i]));
            header.add(this.header.get(columnIndex[i]));
        }
        return new ColumnTableTool<>(header, data);
    }
}
