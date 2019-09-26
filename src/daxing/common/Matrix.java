package daxing.common;

import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * have column name, but header can be false
 */
public class Matrix {

    private double[][] data;
    private List<String> id;

    public Matrix(String inputFile, boolean header){
        this.readFile(inputFile, header);
    }

    private void readFile(String inputFile, boolean header){
        try (BufferedReader br = IOUtils.getTextReader(inputFile)) {
            this.id=new ArrayList<>();
            List<double[]> dataList=new ArrayList<>();
            String line;
            List<String> lineList;
            while ((line=br.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                id.add(lineList.get(0));
                dataList.add(lineList.stream().skip(1).mapToDouble(Double::parseDouble).toArray());
            }
            data=new double[dataList.size()][];
            for (int i = 0; i < dataList.size(); i++) {
                data[i]=dataList.get(i);
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
            sb.append("\t").append(String.join("\t", this.getId()));
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

}
