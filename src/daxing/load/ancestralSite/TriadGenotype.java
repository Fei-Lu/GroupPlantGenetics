package daxing.load.ancestralSite;

import daxing.common.ColumnTableTool;
import daxing.common.IOTool;
import daxing.common.NumberTool;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class TriadGenotype {

    public static String[] loadGenotypeArraySorted ={"HHH","HHL","HLH","HLL","LHH","LHL","LLH","LLL"};
    public static String[] modelSorted ={"M001","M010","M011","M100","M101","M110","M111"};
    public static double[] expectedModelRatio={ 1d/12, 1d/12, 2d/12, 1d/12, 2d/12, 2d/12, 3d/12};
    public static double loadThresh;

    public ColumnTableTool<String> columnTableTool;

    public TriadGenotype(String triadGenotypeFile){
        this.columnTableTool=new ColumnTableTool<>(triadGenotypeFile);
    }

    public TriadGenotype(ColumnTableTool<String> columnTable, double loadThresh){
        this.columnTableTool=columnTable;
        TriadGenotype.loadThresh=loadThresh;
    }

    public void write(String outFile, IOFileFormat ioFileFormat){
        this.columnTableTool.writeTextTable(outFile, ioFileFormat);
    }

    public static TriadGenotype mergeHexaploid(String inputDir, double loadThresh){
        List<File> files= IOTool.getVisibleDir(inputDir);
        ColumnTableTool<String> columnTable=null;
        List<String> loadAColumn, loadBColumn, loadDColumn, loadModelColumn, triadID;
        List<List<String>> cells=new ArrayList<>(410);
        List<String> loadABD;
        StringBuilder sb;
        StringBuilder headerSB=new StringBuilder();
        headerSB.append("Triad");
        String taxonName, loadA, loadB, loadD, loadGenotype;
        String loadAGenotype, loadBGenotype, loadDGenotype;
        for (int i = 0; i < files.size(); i++) {
            columnTable=new ColumnTableTool<>(files.get(i).getAbsolutePath());
            taxonName=PStringUtils.fastSplit(files.get(i).getName(),".triad").get(0);
            headerSB.append("\t").append(taxonName);
            loadAColumn=columnTable.getColumn(5);
            loadBColumn=columnTable.getColumn(6);
            loadDColumn=columnTable.getColumn(7);
            loadModelColumn=columnTable.getColumn(8);
            loadABD=new ArrayList<>(20000);
            for (int j = 0; j < loadAColumn.size(); j++) {
                sb=new StringBuilder();
                loadA=loadAColumn.get(j);
                loadB=loadBColumn.get(j);
                loadD=loadDColumn.get(j);
                if (loadA.equals("NA") || loadB.equals("NA") || loadD.equals("NA")){
                    loadGenotype="NA";
                }else {
                    loadAGenotype=Double.parseDouble(loadA) > loadThresh ? "H" : "L";
                    loadBGenotype=Double.parseDouble(loadB) > loadThresh ? "H" : "L";
                    loadDGenotype=Double.parseDouble(loadD) > loadThresh ? "H" : "L";
                    loadGenotype=String.join("", loadAGenotype,loadBGenotype,loadDGenotype);
                }
                sb.append(loadA).append(",").append(loadB).append(",").append(loadD).append(":").append(loadGenotype);
                sb.append(":").append(loadModelColumn.get(j));
                loadABD.add(sb.toString());
            }
            cells.add(loadABD);
        }
        triadID=columnTable.getColumn(0);
        cells.add(0, triadID);
        List<String> headerList= PStringUtils.fastSplit(headerSB.toString());
        ColumnTableTool<String> outTable=new ColumnTableTool<>(headerList, cells);
        List<String> highLoadFrequency=new ArrayList<>();
        String highLoadFre;
        for (int i = 0; i < outTable.getRowNumber(); i++) {
            highLoadFre= calculateHighLoadFrequency(outTable.getRow(i));
            highLoadFrequency.add(highLoadFre);
        }
        outTable.insertColumn("HighLoadFrequency",1, highLoadFrequency);
//        File outFile=new File(outDir, "hexaploidLoadThresh"+loadThresh+".txt.gz");
        return new TriadGenotype(outTable, loadThresh);
//        outTable.writeTextTable(outFile.getAbsolutePath(), IOFileFormat.TextGzip);
    }

    private static String calculateHighLoadFrequency(List<String> row){
        String loadGenotype;
        int[] countLoadGenotype={3,2,2,1,2,1,1,0};
        double countH=0, count=0;
        for (int i = 1; i < row.size(); i++) {
            loadGenotype=PStringUtils.fastSplit(row.get(i), ":").get(1);
            if (loadGenotype.equals("NA")) continue;
            count++;
            int index= Arrays.binarySearch(loadGenotypeArraySorted, loadGenotype);
            countH=countH+countLoadGenotype[index];
        }
        if (count==0) return "NA";
        return String.valueOf(NumberTool.format(countH/(count*3), 5));
    }

    public void writeModel(String outFile){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
//            List<String> header=columnTableTool.getHeader();
            bw.write(String.join("\t", modelSorted));
            bw.newLine();
            StringBuilder sb;
            List<String> temp;
            int[] count;
            for (int i = 0; i < columnTableTool.getRowNumber(); i++) {
                sb=new StringBuilder();
                sb.append(columnTableTool.getCell(i, 0)).append("\t");
                sb.append(columnTableTool.getCell(i, 1)).append("\t");
                count=new int[7];
                for (int j = 2; j < columnTableTool.getRow(i).size(); j++) {
                    temp=PStringUtils.fastSplit(columnTableTool.getCell(i,j),":");
                    int index=Arrays.binarySearch(modelSorted, temp.get(2));
                    if (index < 0) continue;
                    count[index]=count[index]+1;
                }
                for (int j = 0; j < count.length; j++) {
                    sb.append(count[j]).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public void writeSelectionCoefficientModel(String outFile){
        String header="Triad\tHighFitnessModel\tM001\tM010\tM011\tM100\tM101\tM110\tM111";
        ColumnTableTool<String> columnTable=this.columnTableTool;
        String[] na={"NA","NA","NA","NA","NA","NA","NA","NA"};
        StringBuilder sb;
        String triadID, highFitnessModel;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < columnTable.getRowNumber(); i++) {
                triadID=columnTable.getRow(i).get(0);
                sb=new StringBuilder();
                sb.append(triadID).append("\t");
                if (ifAllNA(columnTable.getRow(i))){
                    sb.append(String.join("\t", na));
                    bw.write(sb.toString());
                    bw.newLine();
                }else {
                    highFitnessModel=calculateHighFitnessModel(columnTable.getRow(i));
                    sb.append(highFitnessModel).append("\t");
                    double[] oeRatio=calculateOERatio(columnTable.getRow(i));
                    int maxIndex=IntStream.range(0, oeRatio.length).boxed().max(Comparator.comparing(e->oeRatio[e])).get();
                    for (int j = 0; j < oeRatio.length; j++) {
                        sb.append(1-(oeRatio[j]/oeRatio[maxIndex])).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private boolean ifAllNA(List<String> row){
        List<String> temp;
        Set<String> set=new HashSet<>();
        for (int i = 2; i < row.size(); i++) {
            temp=PStringUtils.fastSplit(row.get(i), ":");
            set.add(temp.get(2));
        }
        Set<String> naSet=new HashSet<>();
        naSet.add("NA");
        return set.equals(naSet);
    }

    private String calculateHighFitnessModel(List<String> row){
        double[] oeRatio=calculateOERatio(row);
        int maxIndex=IntStream.range(0, oeRatio.length).boxed().max(Comparator.comparing(e->oeRatio[e])).get();
        return modelSorted[maxIndex];
    }

    private double[] calculateOERatio(List<String> row){
        List<String> temp;
        int index;
        int total=0;
        double[] count=new double[7];
        for (int i = 2; i < row.size(); i++) {
            temp=PStringUtils.fastSplit(row.get(i), ":");
            index=Arrays.binarySearch(modelSorted, temp.get(2));
            if (index < 0) continue;
            total++;
            count[index]=count[index]+1;
        }
        double[] expectedCount=new double[7];
        for (int i = 0; i < 7; i++) {
            expectedCount[i]=total*expectedModelRatio[i];
        }
        double[] oeRatio=new double[7];
        for (int i = 0; i < 7; i++) {
            oeRatio[i]=count[i]/expectedCount[i];
        }

        return oeRatio;
    }

    public void writeSelectionCoefficientLH(String outFile){
        String header="Triad\tHighFitnessGenotypeLH\tHHH\tHHL\tHLH\tHLL\tLHH\tLHL\tLLH\tLLL";
        ColumnTableTool<String> columnTable=this.columnTableTool;
        String[] na={"NA","NA","NA","NA","NA","NA","NA","NA"};
        StringBuilder sb;
        String triadID, highLoadFrequency, highFitnessGenotypeLH;
        double[] selectionCoefficientLH;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write(header);
            bw.newLine();
            for (int i = 0; i < columnTable.getRowNumber(); i++) {
                sb=new StringBuilder();
                triadID=columnTable.getRow(i).get(0);
                highFitnessGenotypeLH=calculateHighFitnessGenotypeLH(columnTable.getRow(i));
                highLoadFrequency=columnTable.getRow(i).get(1);
                sb.append(triadID).append("\t").append(highFitnessGenotypeLH).append("\t");
                if (highLoadFrequency.equals("NA") || Double.parseDouble(highLoadFrequency)==0){
                    sb.append(String.join("\t", na));
                    bw.write(sb.toString());
                    bw.newLine();
                }else {
                    double fre= Double.parseDouble(highLoadFrequency);
                    selectionCoefficientLH=calculateSelectionCoefficientLH(columnTable.getRow(i), fre);
                    for (int j = 0; j < selectionCoefficientLH.length; j++) {
                        sb.append(selectionCoefficientLH[j]).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }



    private String calculateHighFitnessGenotypeLH(List<String> row){
        String highLoadFrequency=row.get(1);
        if (highLoadFrequency.equals("NA") || Double.parseDouble(highLoadFrequency)==0) return "NA";
        double[] oeRatio=calculateOERatioLH(row, Double.parseDouble(highLoadFrequency));
        int indexMaxOERatio= IntStream.range(0, oeRatio.length).boxed().max(Comparator.comparing(e->oeRatio[e])).get();
        return loadGenotypeArraySorted[indexMaxOERatio];
    }

    private double[] calculateSelectionCoefficientLH(List<String> row, double highLoadFrequency){
        double[] oeRatio=calculateOERatioLH(row, highLoadFrequency);
        int indexMaxOERatio= IntStream.range(0, oeRatio.length).boxed().max(Comparator.comparing(e->oeRatio[e])).get();
        double[] s=new double[8];
        for (int i = 0; i < 8; i++) {
            s[i]=1-(oeRatio[i]/oeRatio[indexMaxOERatio]);
        }
        return s;
    }

    private double[] calculateOERatioLH(List<String> row, double highLoadFrequency){
        double[] countObserved=new double[8];
        double[] countExpected=null;
        List<String> temp;
        int index, total=0;
        for (int i = 2; i < row.size(); i++) {
            temp=PStringUtils.fastSplit(row.get(i), ":");
            index=Arrays.binarySearch(loadGenotypeArraySorted, temp.get(1));
            if (index < 0) continue;
            countObserved[index]=countObserved[index]+1;
            total++;
        }
        countExpected=countExpectedLoadGenotype(highLoadFrequency, total);
        double[] oeRatio=new double[8];
        for (int i = 0; i < 8; i++) {
            oeRatio[i]=countObserved[i]/countExpected[i];
        }
        return oeRatio;
    }

    private double[] countExpectedLoadGenotype(double highLoadFrequency, int taxonNum){
        double[] expectedLoadGenotypedRatio=new double[8];
        double hhh=highLoadFrequency*highLoadFrequency*highLoadFrequency;
        double lll=(1-highLoadFrequency)*(1-highLoadFrequency)*(1-highLoadFrequency);
        double lhh=(1-highLoadFrequency)*highLoadFrequency*highLoadFrequency;
        double hll=highLoadFrequency*(1-highLoadFrequency)*(1-highLoadFrequency);
        expectedLoadGenotypedRatio[0]=hhh;
        expectedLoadGenotypedRatio[1]=lhh;
        expectedLoadGenotypedRatio[2]=lhh;
        expectedLoadGenotypedRatio[3]=hll;
        expectedLoadGenotypedRatio[4]=lhh;
        expectedLoadGenotypedRatio[5]=hll;
        expectedLoadGenotypedRatio[6]=hll;
        expectedLoadGenotypedRatio[7]=lll;
        double[] expectedCount=new double[8];
        for (int i = 0; i < expectedCount.length; i++) {
            expectedCount[i]=taxonNum*expectedLoadGenotypedRatio[i];
        }
        return expectedCount;
    }
}
