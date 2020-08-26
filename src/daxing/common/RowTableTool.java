package daxing.common;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import pgl.infra.table.RowTable;
import org.apache.commons.lang.StringUtils;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Consumer;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * @author Daxing Xu
 * @param <T>
 */
public class RowTableTool<T> extends RowTable<T> {

    public RowTableTool(String infileS){
        super(infileS);
    }

    public RowTableTool(String infileS, String delimiter){
        super(infileS, delimiter);
    }

    public RowTableTool(List<String> header, List<List<T>> cells){
        super(header, cells);
    }

    public List<List<T>> getCells(){
        return this.cells;
    }

    /**
     * Removes all of the rows of this rowTable that satisfy the given predicate
     * @param p 针对表格进行过滤的函数式接口
     * @return true, if any row were removed
     */
    public boolean removeIf(Predicate<List<T>> p){
        List<List<T>> cells=this.cells;
        boolean f=cells.removeIf(p);
        this.cells=cells;
        return f;
    }

    /**
     * Performs the given action for each row
     * @param consumer The action to be performed for each element
     */
    public void forEach(Consumer<List<T>> consumer){
        List<List<T>> cells=this.cells;
        cells.forEach(consumer);
    }

    public RowTableTool<T> filter(Predicate<List<T>> p){
        List<List<T>> cells=this.cells;
        List<List<T>> newCells=cells.stream().filter(p).collect(Collectors.toList());
        List<String> header=this.header;
        return new RowTableTool<>(header, newCells);
    }

    /**
     *
     * @param comparator 针对表格进行排序的比较器对象
     */
    public void sortBy(Comparator<List<T>> comparator){
        this.cells.sort(comparator);
    }

    public Map<T, T> getHashMap(int columnIndexA, int columnIndexB){
        List<T> columnA=this.getColumn(columnIndexA);
        List<T> columnB=this.getColumn(columnIndexB);
        Map<T, T> map=new HashMap<>();
        for (int i = 0; i < columnA.size(); i++) {
            map.put(columnA.get(i), columnB.get(i));
        }
        return map;
    }

    public boolean contain(String str, int columnIndex){
        List<String> column= (List<String>) this.getColumn(columnIndex);
        return column.contains(str);
    }

    public void write(String outfile){
        this.writeTextTable(outfile, IOFileFormat.Text);
    }

    public void write(File outFile){
        write(outFile.getAbsolutePath());
    }

    public void write(String outFile, IOFileFormat ioFileFormat){
        writeTextTable(outFile, ioFileFormat);
    }

    public void write(File outFile, IOFileFormat ioFileFormat){
        write(outFile.getAbsolutePath(), ioFileFormat);
    }

    /**
     * 表格合并
     * @param rowTable
     * @return
     */
    public boolean add(RowTableTool<T> rowTable){
        List<String> header_1=this.header;
        List<String> header_2=rowTable.getHeader();
        if (! header_1.equals(header_2)){
            System.out.println("error, two rowTables were not equal");
            System.exit(1);
        }
        List<List<T>> data=this.cells;
        boolean res=data.addAll(rowTable.cells);
        this.cells=data;
        return res;
    }

    /**
     *
     * @param columnKeyA
     * @param columnKeyB
     * @param columnValue
     * @return
     */
    public Table<T, T, T> getTable(int columnKeyA, int columnKeyB, int columnValue){
        Table<T, T, T> table=HashBasedTable.create();
        List<List<T>> data=this.cells;
        T key1, key2, value;
        for (int i = 0; i < data.size(); i++) {
            key1=data.get(i).get(columnKeyA);
            key2=data.get(i).get(columnKeyB);
            value=data.get(i).get(columnValue);
            table.put(key1, key2, value);
        }
        return table;
    }

    /**
     * 默认前两列为key
     * @param columnIndex value
     * @return
     */
    public Table<T,T,T> getTable(int columnIndex){
        return getTable(0, 1, columnIndex);
    }

    /**
     *
     * @param columnName
     * @return
     */
    public double[] getColumnAsDoubleArray(String columnName){
        int columnIndex=this.hiMap.get(columnName);
        return this.getColumnAsDoubleArray(columnIndex);
    }

    /**
     *
     * @param columnName
     * @return
     */
    public int[] getColumnAsIntArray(String columnName){
        int columnIndex=this.hiMap.get(columnName);
        return this.getColumnAsIntArray(columnIndex);
    }

    /**
     * extract the specific column for big table
     * @param inFile
     * @param columnIndex
     * @return
     */
    public static Set<String> getColumnSet(String inFile, int columnIndex){
        Set<String> s=new HashSet<>();
        try (BufferedReader br = IOTool.getReader(inFile)) {
            String line;
            String[] temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= StringUtils.split(line);
                s.add(temp[columnIndex]);
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return s;
    }

    public static List<String> getColumnList(String inFile, int columnIndex, String sep){
        List<String> list=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(inFile)) {
            String line;
            String[] temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= StringUtils.split(line, sep);
                list.add(temp[columnIndex]);
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return list;
    }

    public static List<String> getColumnList(String inFile, int columnIndex){
        return getColumnList(inFile, columnIndex, "\t");
    }

    public static Map<String,String> getMap(String tableFile, int columnIndexA, int columnIndexB, String sep){
        List<String> listA=getColumnList(tableFile, columnIndexA, sep);
        List<String> listB=getColumnList(tableFile, columnIndexB, sep);
        Map<String,String> map=new HashMap<>();
        for (int i = 0; i < listA.size(); i++) {
            map.put(listA.get(i), listB.get(i));
        }
        return map;
    }

    public static Map<String,String> getMap(String tableFile, int columnIndexA, int columnIndexB){
        return getMap(tableFile, columnIndexA, columnIndexB, "\t");
    }

    /**
     * merge multiple table into one table with a new column
     * note: each row name of the new column is equal to table file name without file extension
     * @param tablesInputDir
     * @param newColumnName new column name
     * @param outFile
     */
    public static void meregMultipleTable(String tablesInputDir, String newColumnName, String outFile){
        File[] files=new File(tablesInputDir).listFiles();
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate().and(File::isFile)).sorted().toArray(File[]::new);
        String[] groupNames= Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("[.].+$", ""))
                                    .toArray(String[]::new);
        BufferedReader br;
        BufferedWriter bw;
        if (outFile.endsWith("gz")){
            bw=IOTool.getTextGzipWriter(outFile);
        }else {
            bw=IOTool.getWriter(outFile);
        }
        String line;
        StringBuilder sb;
        String header;
        boolean first=true;
        try {
            for (int i = 0; i < f1.length; i++) {
                br=IOTool.getReader(f1[i].getAbsolutePath());
                header=br.readLine();
                if (first){
                    bw.write(header+"\t"+newColumnName+"\n");
                    first=false;
                }
                while ((line=br.readLine())!=null){
                    sb=new StringBuilder();
                    sb.append(line).append("\t").append(groupNames[i]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * tableFile
     * Chr	Pos	Ancestral	Ref	Urartu	Barley
     * 1	304	C	C	C	C
     * 1	305	A	A	A	A
     * 1	306	C	C	C	C
     * 1	307	G	G	G	G
     * @param tableFile
     * @param columnKey1
     */
    public static Table<String, String, String> getTable(String tableFile, int columnKey1, int columnKey2, int columnValue){
        BufferedReader br;
        if (tableFile.endsWith("gz")){
            br=IOUtils.getTextGzipReader(tableFile);
        }else {
            br=IOUtils.getTextReader(tableFile);
        }
        Table<String, String, String> table= HashBasedTable.create();
        try {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                table.put(temp.get(columnKey1), temp.get(columnKey2), temp.get(columnValue));
            }
            br.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return table;
    }

    /**
     *
     * @param tableFile
     * @param columnIndex
     * @return
     */
    public static Table<String, String, String> getTable(String tableFile, int columnIndex){
        return getTable(tableFile, 0, 1, columnIndex);
    }

    /**
     * merge multiple tables
     * @param rowTableFiles
     * @param rowTableOutFile
     */
    public static void mergeRowTables(List<File> rowTableFiles, String rowTableOutFile){
        BufferedReader br;
        BufferedWriter bw;
        String header, line;
        try {
            br=IOTool.getReader(rowTableFiles.get(0));
            bw=IOUtils.getTextWriter(rowTableOutFile);
            header=br.readLine();
            bw.write(header);
            bw.newLine();
            while ((line=br.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            br.close();
            for (int i = 1; i < rowTableFiles.size(); i++) {
                br=IOTool.getReader(rowTableFiles.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    bw.write(line);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
