package daxing.load;

import daxing.common.ColumnTableTool;
import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class TriadGenotype {

    public static String[] loadGenotypeArray={"HHH","HHL","HLH","HLL","LHH","LHL","LLH","LLL"};

    public static void mergeHexaploid(String inputDir, double loadThresh, String outDir, String vmapIIGroupFile){
        List<File> files= IOTool.getVisibleDir(inputDir);
        RowTableTool<String> groupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> groupMap=groupTable.getHashMap(0,15);
        Predicate<File> cultivar= f->groupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals(
                "OthrerHexaploid");
        List<File> notCultivarFiles=files.stream().collect(Collectors.toList());
        ColumnTableTool<String> columnTable=null;
        List<String> loadAColumn, loadBColumn, loadDColumn, loadModelColumn, triadID;
        List<List<String>> cells=new ArrayList<>(410);
        List<String> loadABD;
        StringBuilder sb;
        StringBuilder headerSB=new StringBuilder();
        headerSB.append("Triad");
        String taxonName, loadA, loadB, loadD, loadGenotype;
        String loadAGenotype, loadBGenotype, loadDGenotype;
        for (int i = 0; i < notCultivarFiles.size(); i++) {
            columnTable=new ColumnTableTool<>(notCultivarFiles.get(i).getAbsolutePath());
            taxonName=PStringUtils.fastSplit(notCultivarFiles.get(i).getName(),".").get(0);
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
        List<String> minorLoadFrequency=new ArrayList<>();
        for (int i = 0; i < outTable.getRowNumber(); i++) {
            minorLoadFrequency.add(calculateMinorLoadFrequency(outTable.getRow(i)));
        }
        outTable.insertColumn("HighLoadFrequency",1, minorLoadFrequency);
        File outFile=new File(outDir, "hexaploidLoadThresh"+loadThresh+".txt.gz");
        outTable.writeTextTable(outFile.getAbsolutePath(), IOFileFormat.TextGzip);
    }

    private static String calculateMinorLoadFrequency(List<String> row){
        String loadGenotype;
        int[] countLoadGenotype={3,2,2,1,2,1,1,0};
        double countH=0, count=0;
        for (int i = 1; i < row.size(); i++) {
            loadGenotype=PStringUtils.fastSplit(row.get(i), ":").get(1);
            if (loadGenotype.equals("NA")) continue;
            count++;
            int index= Arrays.binarySearch(loadGenotypeArray, loadGenotype);
            countH=countH+countLoadGenotype[index];
        }
        if (count==0) return "NA";
        return String.valueOf(NumberTool.format(countH/(count*3), 5));
    }
}
