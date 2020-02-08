package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.stream.IntStream;

/**
 * have column name, but header can be false;
 * dxy, r2
 * 默认的数据类型为double, 可以以int类型输出
 * @author Daxing Xu
 */
public class Matrix {

    private double[][] data;
    private List<String> id;  //id 不能重复

    public Matrix(String inputFile, boolean header){
        this.readFile(inputFile, header);
    }

    public Matrix(List<String> id, double[][] data){
        this.id=id;
        this.data=data;
    }

    private void readFile(String inputFile, boolean header){
        try (BufferedReader br = IOUtils.getTextReader(inputFile)) {
            List<double[]> dataList=new ArrayList<>();
            String line;
            List<String> lineList;
            if (header){
                List<String> names=PStringUtils.fastSplit(br.readLine());
                names.remove(0);
                this.id=names;
                while ((line=br.readLine())!=null){
                    lineList= PStringUtils.fastSplit(line);
                    dataList.add(lineList.stream().skip(1).mapToDouble(Double::parseDouble).toArray());
                }
                data=new double[dataList.size()][];
                for (int i = 0; i < dataList.size(); i++) {
                    data[i]=dataList.get(i);
                }
            }else {
                this.id=new ArrayList<>();
                while ((line=br.readLine())!=null){
                    lineList= PStringUtils.fastSplit(line);
                    id.add(lineList.get(0));
                    dataList.add(lineList.stream().skip(1).mapToDouble(Double::parseDouble).toArray());
                }
                data=new double[dataList.size()][];
                for (int i = 0; i < dataList.size(); i++) {
                    data[i]=dataList.get(i);
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    public double[][] getData() {
        return data;
    }

    public List<String> getId() {
        return id;
    }

    /**
     *
     * @param taxon
     * @return
     */
    public double[] getRow(String taxon){
        int index=this.indexOf(taxon);
        return this.getData()[index];
    }

    public int indexOf(String taxon){
        return this.getId().indexOf(taxon);
    }

    /**
     *
     * @param id1
     * @param id2
     * @return
     */
    public double getValue(String id1, String id2){
        int index1=this.getId().indexOf(id1);
        int index2=this.getId().indexOf(id2);
        return this.getData()[index1][index2];
    }

    public void removeTaxons(String[] taxons){
        TIntArrayList taxonsIndex=new TIntArrayList(IntStream.range(0, this.getId().size()).toArray());
        int[] taxonsIndexWillBeRemoved=new int[taxons.length];
        List<String> id=this.getId();
        for (int i = 0; i < taxons.length; i++) {
            taxonsIndexWillBeRemoved[i]=id.indexOf(taxons[i]);
            if (taxonsIndexWillBeRemoved[i]<0){
                System.out.println("please check your parameter taxons "+taxons[i]+" did not contain in this matric");
                System.exit(1);
            }
        }
        taxonsIndex.removeAll(taxonsIndexWillBeRemoved);
        double[][] data=this.getData();
        List<Double> temp;
        List<double[]> dataRetainList=new ArrayList<>();
        for (int i = 0; i < taxonsIndex.size(); i++) {
            temp=new ArrayList<>();
            for (int j = 0; j < taxonsIndex.size(); j++) {
                temp.add(data[taxonsIndex.get(i)][taxonsIndex.get(j)]);
            }
            dataRetainList.add(temp.stream().mapToDouble(Double::doubleValue).toArray());
        }
        double[][] dataRetain=new double[dataRetainList.size()][];
        for (int i = 0; i < dataRetainList.size(); i++) {
            dataRetain[i]=dataRetainList.get(i);
        }
        this.data=dataRetain;
        List<String> idRetain=new ArrayList<>();
        for (int i = 0; i < taxonsIndex.size(); i++) {
            idRetain.add(this.getId().get(taxonsIndex.get(i)));
        }
        this.id=idRetain;
    }

    public Matrix getSubMatrix(String[] subID){
        Set<String> s=new HashSet<>(CollectionTool.changeToList(subID));
        if (s.size()<subID.length){
            System.out.println("subID parameter can not contain duplicated ID");
            System.exit(1);
        }
        List<Integer> indexList=new ArrayList<>();
        int index=-1;
        for (int i = 0; i < subID.length; i++) {
            index=this.getId().indexOf(subID[i]);
            if (index<0){
                System.out.println("please check your subID parameter, "+subID[i]+" not contain in matrix");
                System.exit(1);
            }
            indexList.add(index);
        }
        double[][] data=this.getData();
        List<String> id=this.getId();
        List<String> subid= new ArrayList<>();
        List<double[]> dataList=new ArrayList<>();
        List<Double> doubleList;
        for (int i = 0; i < indexList.size(); i++) {
            doubleList=new ArrayList<>();
            subid.add(id.get(indexList.get(i)));
            for (int j = 0; j < indexList.size(); j++) {
                doubleList.add(data[indexList.get(i)][indexList.get(j)]);
            }
            dataList.add(doubleList.stream().mapToDouble(Double::doubleValue).toArray());
        }
        data=new double[dataList.size()][];
        for (int i = 0; i < dataList.size(); i++) {
            data[i]=dataList.get(i);
        }
        return new Matrix(subid, data);
    }

    public void sortBy(String[] sortedID){
        int size= (int) Arrays.stream(sortedID).distinct().count();
        if (sortedID.length!=size || sortedID.length!= this.getId().size()){
            System.out.println("sortedID array is not right, program will quit");
            System.exit(1);
        }
        List<Integer> indexList=new ArrayList<>();
        for (int i = 0; i < sortedID.length; i++) {
            indexList.add(this.id.indexOf(sortedID[i]));
        }
        double[][] data=this.getData();
        List<String> id=this.getId();
        List<String> sortedIDList= new ArrayList<>();
        List<double[]> dataList=new ArrayList<>();
        List<Double> doubleList;
        for (int i = 0; i < indexList.size(); i++) {
            doubleList=new ArrayList<>();
            sortedIDList.add(id.get(indexList.get(i)));
            for (int j = 0; j < indexList.size(); j++) {
                doubleList.add(data[indexList.get(i)][indexList.get(j)]);
            }
            dataList.add(doubleList.stream().mapToDouble(Double::doubleValue).toArray());
        }
        data=new double[dataList.size()][];
        for (int i = 0; i < dataList.size(); i++) {
            data[i]=dataList.get(i);
        }
        this.data=data;
        this.id=sortedIDList;
    }

    /**
     * 按照sortedID的顺序输出
     * @param sortedID 不能有重复值，且不能含有对象ID之外的ID
     * @param outFile
     */
    public void WriteBySortedID(String[] sortedID, String outFile){
        int size= (int) Arrays.stream(sortedID).distinct().count();
        if (sortedID.length!=size || sortedID.length!= this.getId().size()){
            System.out.println("sortedID array is not right, program will quit");
            System.exit(1);
        }
        List<Integer> indexList=new ArrayList<>();
        for (int i = 0; i < sortedID.length; i++) {
            indexList.add(this.id.indexOf(sortedID[i]));
        }
        int[] index=indexList.stream().mapToInt(Integer::intValue).toArray();
        try {
            BufferedWriter bw=IOUtils.getTextWriter(outFile);
            StringBuilder sb=new StringBuilder();
            sb.append("ID\t");
            for (int i = 0; i <index.length ; i++) {
                sb.append(this.getId().get(index[i])).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < indexList.size(); i++) {
                sb=new StringBuilder();
                sb.append(this.getId().get(index[i]));
                for (int j = 0; j < indexList.size(); j++) {
                    sb.append("\t").append(data[index[i]][index[j]]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void write(String outFile){
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("ID\t").append(String.join("\t", this.getId()));
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < data.length; i++) {
                sb=new StringBuilder();
                sb.append(this.getId().get(i));
                for (int j = 0; j < data[i].length; j++) {
                    sb.append("\t").append(data[i][j]);
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     * 以整数类型输出
     * @param outFile
     */
    public void writeAsInteger(String outFile){
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("ID\t").append(String.join("\t", this.getId()));
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < data.length; i++) {
                sb=new StringBuilder();
                sb.append(this.getId().get(i));
                for (int j = 0; j < data[i].length; j++) {
                    sb.append("\t").append(Math.round(data[i][j]));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public static void sortMatrixBySortedGroup(String matirxFile, String groupDir, String[] groupNames, String outFile){
        Map<String, Integer> map=new HashMap<>();
        for (int i = 0; i < groupNames.length; i++) {
            map.put(groupNames[i], i);
        }
        File[] files=IOUtils.listRecursiveFiles(new File(groupDir));
        Comparator<File> comparator=Comparator.comparing(f->map.get(f.getName().replaceAll("\\.txt$", "")));
        Arrays.sort(files, comparator);
        String[] ids= Arrays.stream(files).map(File::getAbsolutePath).map(IOUtils::getTextReader).flatMap(BufferedReader::lines)
                .toArray(String[]::new);
        Matrix m=new Matrix(matirxFile, false);
        m.WriteBySortedID(ids, outFile);
    }

    /**
     * 将某矩阵的某一行写出
     * @param dxyMatrix
     * @param taxon
     * @param outFile
     */
    public static void writeFordxyGraph(String dxyMatrix, String taxon, String outFile){
        Matrix matrix=new Matrix(dxyMatrix, false);
        double[] cs=matrix.getRow(taxon);
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write("taxon\tdxy");
            bw.newLine();
            for (int i = 0; i < cs.length-1; i++) {
                bw.write(matrix.getId().get(i)+"\t"+ cs[i]);
                bw.newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }

    }

}
