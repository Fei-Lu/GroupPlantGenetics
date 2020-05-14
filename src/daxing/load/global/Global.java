package daxing.load.global;

import daxing.common.ColumnTableTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.load.Standardization;
import pgl.infra.table.ColumnTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class Global {

    public static void global(String inputDir, String vmapIIGroupFile, String outDir, String nonsynOrDel){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=vmapIIGroupTable.getHashMap(0,15);
        Predicate<File> cultivarPredict= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Cultivar");
        Predicate<File> landraceEUPredict= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Landrace_Europe");
        List<File> cultivarFiles=files.stream().filter(cultivarPredict).collect(Collectors.toList());
        List<File> landraceEUFiles=files.stream().filter(landraceEUPredict).collect(Collectors.toList());
        List<ColumnTable<String>> cultivarTables=new ArrayList<>();
        List<ColumnTable<String>> landraceEUTables=new ArrayList<>();
        RowTableTool<String> rowTableTool;
        ColumnTable<String> columnTable;
        List<String> triadID, loadA, loadB, loadD;
        List<List<String>> cells;
        String header="TriadID\tloadA\tloadB\tloadD";
        int[] indexNonsyn={5,6,7};
        int[] indexDel={9,10,11};
        int[] indexSyn={1,2,3};
        int[] index=null;
        if (nonsynOrDel.equals("nonsyn")){
            index=indexNonsyn;
        }else if (nonsynOrDel.equals("del")){
            index=indexDel;
        }else if (nonsynOrDel.equals("syn")){
            index=indexSyn;
        }
        for (int i = 0; i < cultivarFiles.size(); i++) {
            rowTableTool=new RowTableTool<>(cultivarFiles.get(i).getAbsolutePath());
            triadID=rowTableTool.getColumn(0);
            loadA=rowTableTool.getColumn(indexNonsyn[0]);
            loadB=rowTableTool.getColumn(indexNonsyn[1]);
            loadD=rowTableTool.getColumn(index[2]);
            cells=new ArrayList<>();
            cells.add(triadID);
            cells.add(loadA);
            cells.add(loadB);
            cells.add(loadD);
            columnTable=new ColumnTable<>(PStringUtils.fastSplit(header), cells);
            cultivarTables.add(columnTable);
        }
        for (int i = 0; i < landraceEUFiles.size(); i++) {
            rowTableTool=new RowTableTool<>(landraceEUFiles.get(i).getAbsolutePath());
            triadID=rowTableTool.getColumn(0);
            loadA=rowTableTool.getColumn(index[0]);
            loadB=rowTableTool.getColumn(index[1]);
            loadD=rowTableTool.getColumn(index[2]);
            cells=new ArrayList<>();
            cells.add(triadID);
            cells.add(loadA);
            cells.add(loadB);
            cells.add(loadD);
            columnTable=new ColumnTable<>(PStringUtils.fastSplit(header), cells);
            landraceEUTables.add(columnTable);
        }
        RowTableTool<String> cultivar=average(cultivarTables);
        RowTableTool<String> landraceEU=average(landraceEUTables);
        cultivar.write(new File(outDir, "cultivar.txt.gz"), IOFileFormat.TextGzip);
        landraceEU.write(new File(outDir, "landraceEU.txt.gz"),IOFileFormat.TextGzip);
        List<File> cultivarLandraceFiles=new ArrayList<>();
        cultivarLandraceFiles.add(new File(outDir, "cultivar.txt.gz"));
        cultivarLandraceFiles.add(new File(outDir, "landraceEU.txt.gz"));
        ColumnTableTool<String> cultivarColumnTable=new ColumnTableTool<>(cultivarLandraceFiles.get(0).getAbsolutePath());
        ColumnTableTool<String> landraceEUColumnTable= new ColumnTableTool<>(cultivarLandraceFiles.get(1).getAbsolutePath());
        List<ColumnTable<String>> globalColumn=new ArrayList<>();
        globalColumn.add(cultivarColumnTable);
        globalColumn.add(landraceEUColumnTable);
        RowTableTool<String> global=average(globalColumn);
        List<List<String>> cultivarColumns=cultivarColumnTable.getCells();
        List<List<String>> landraceEUColumns=landraceEUColumnTable.getCells();
        String cultivarHeader="CultivarLoadA\tCultivarLoadB\tCultivarLoadD\tCultivarRegion";
        String landraceEUHeader="Landrace_EULoadA\tLandrace_EULoadB\tLandrace_EULoadD\tLandrace_EURegion";
        for (int i = 1; i < cultivarColumns.size(); i++) {
            global.addColumn(PStringUtils.fastSplit(cultivarHeader).get(i-1), cultivarColumns.get(i));
        }
        for (int i = 1; i < landraceEUColumns.size() ; i++) {
            global.addColumn(PStringUtils.fastSplit(landraceEUHeader).get(i-1), landraceEUColumns.get(i));
        }
        addDistance(global);
        global.write(new File(outDir, nonsynOrDel+"Global.txt.gz"), IOFileFormat.TextGzip);
        new File(outDir, "cultivar.txt.gz").delete();
        new File(outDir, "landraceEU.txt.gz").delete();
    }

    private static RowTableTool<String> average(List<ColumnTable<String>> tableList){
        List<String> a,b,d;
        String averageA, averageB, averageD, region;
        List<List<String>> cells=new ArrayList<>();
        List<String> cell;
        double[] totalABD;
        int size= tableList.get(0).getRowNumber();
        for (int i = 0; i < tableList.get(0).getRowNumber(); i++) {
            a=new ArrayList<>();
            b=new ArrayList<>();
            d=new ArrayList<>();
            int aSize=0,bSize=0,dSize=0;
            double totalA=0,totalB=0,totalD=0;
            for (int j = 0; j < tableList.size(); j++) {
                a.add(tableList.get(j).getCell(i, 1));
                b.add(tableList.get(j).getCell(i, 2));
                d.add(tableList.get(j).getCell(i, 3));
            }
            for (int j = 0; j < a.size(); j++) {
                if (a.get(j).equals("NA")) continue;
                aSize++;
                totalA+=Double.parseDouble(a.get(j));
            }
            for (int j = 0; j < b.size(); j++) {
                if (b.get(j).equals("NA")) continue;
                bSize++;
                totalB+=Double.parseDouble(b.get(j));
            }
            for (int j = 0; j < d.size(); j++) {
                if (d.get(j).equals("NA")) continue;
                dSize++;
                totalD+=Double.parseDouble(d.get(j));
            }
            averageA = aSize==0 ? "NA" : String.valueOf(NumberTool.format(totalA/aSize, 5));
            averageB = bSize==0 ? "NA" : String.valueOf(NumberTool.format(totalB/bSize, 5));
            averageD = dSize==0 ? "NA" : String.valueOf(NumberTool.format(totalD/dSize, 5));
            totalABD=new double[3];
            totalABD[0]=totalA/aSize;
            totalABD[1]=totalB/bSize;
            totalABD[2]=totalD/dSize;
            if (averageA.equals("NA") || averageB.equals("NA") || averageD.equals("NA")){
                region="NA";
            }else {
                region=Standardization.getNearestPointIndex(totalABD).getRegion();
            }
            cell=new ArrayList<>();
            cell.add(tableList.get(0).getCell(i, 0));
            cell.add(averageA);
            cell.add(averageB);
            cell.add(averageD);
            cell.add(region);
            cells.add(cell);
        }
        String header="TriadID\tA\tB\tD\tRegion";
        return new RowTableTool<>(PStringUtils.fastSplit(header), cells);
    }

    private static void addDistance(RowTableTool<String> global){
        double[] abdGlobal, abdC, abdLR;
        List<String> distanceCultivar=new ArrayList<>(), distanceLandraceEU=new ArrayList<>();
        double cDistance, lrDistance;
        String distanceC, distanceLR;
        for (int i = 0; i < global.getRowNumber(); i++) {
            if (global.getCell(i, 4).equals("NA") || global.getCell(i, 8).equals("NA") || global.getCell(i, 12).equals("NA")){
                distanceCultivar.add("NA");
                distanceLandraceEU.add("NA");
                continue;
            }
            if (global.getCell(i, 4).equals("M000") || global.getCell(i, 8).equals("M000") || global.getCell(i, 12).equals("M000")){
                distanceCultivar.add("NA");
                distanceLandraceEU.add("NA");
                continue;
            }
            abdGlobal=new double[3];
            abdC=new double[3];
            abdLR=new double[3];
            abdGlobal[0]=Double.parseDouble(global.getCell(i, 1));
            abdGlobal[1]=Double.parseDouble(global.getCell(i, 2));
            abdGlobal[2]=Double.parseDouble(global.getCell(i, 3));
            abdC[0]=Double.parseDouble(global.getCell(i, 5));
            abdC[1]=Double.parseDouble(global.getCell(i, 6));
            abdC[2]=Double.parseDouble(global.getCell(i, 7));
            abdLR[0]=Double.parseDouble(global.getCell(i, 9));
            abdLR[1]=Double.parseDouble(global.getCell(i, 10));
            abdLR[2]=Double.parseDouble(global.getCell(i, 11));
            cDistance= Standardization.getDistance(abdGlobal, abdC);
            lrDistance=Standardization.getDistance(abdGlobal, abdLR);
            distanceC=String.valueOf(NumberTool.format(cDistance, 5));
            distanceLR=String.valueOf(NumberTool.format(lrDistance, 5));
            distanceCultivar.add(distanceC);
            distanceLandraceEU.add(distanceLR);
        }
        global.insertColumn("CultivarDistance", 9, distanceCultivar);
        global.addColumn("Landrace_EUDistance",  distanceLandraceEU);
    }
}
