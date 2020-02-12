package daxing.common;

import pgl.format.table.ColumnTable;

import java.util.List;

public class ColumnTableTool<T> extends ColumnTable<T> {

    public ColumnTableTool(String infileS) {
        super(infileS);
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
}
