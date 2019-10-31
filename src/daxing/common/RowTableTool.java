package daxing.common;

import format.table.RowTable;
import org.apache.commons.lang.StringUtils;
import utils.IOUtils;
import utils.Tuple;

import java.io.BufferedReader;
import java.util.*;
import java.util.stream.IntStream;

public class RowTableTool<T> extends RowTable<T> {

    public RowTableTool(String infileS){
        super(infileS);
    }

    public Map<T, T> getMap(int columnIndex){
        Set<T> s=new HashSet<>(this.getColumn(0));
        if (s.size()<this.getColumn(0).size()){
            System.out.println("0 column has duplicated keys");
        }
        List<T> key=this.getColumn(0);
        List<T> value=this.getColumn(columnIndex);
        Map<T,T> map=new HashMap<>();
        IntStream.range(0,key.size()).forEach(index-> map.put(key.get(index), value.get(index)));
        return map;
    }

    public Tuple<T[], T[]> getTuple(int columnIndexA, int columnIndexB){
        Set<T> s=new HashSet<>(this.getColumn(columnIndexA));
        if (s.size()<this.getColumn(columnIndexA).size()){
            System.out.println(columnIndexA+" column has duplicated keys");
        }
        List<T> key=this.getColumn(columnIndexA);
        List<T> value=this.getColumn(columnIndexB);
        T[] k= (T[]) key.toArray();
        T[] v= (T[]) value.toArray();
        return new Tuple(k, v);
    }

    /**
     *
     * @param comparator 针对表格进行排序的比较器对象
     */
    public void sortBy(Comparator<List<T>> comparator){
        this.cells.sort(comparator);
    }

    /**
     * extract the specific column for big table
     * @param inFile
     * @param columnIndex
     * @return
     */
    public static Set<String> getColumnSet(String inFile, int columnIndex){
        Set<String> s=new HashSet<>();
        try (BufferedReader br = IOUtils.getTextReader(inFile)) {
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
}
