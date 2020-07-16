package daxing.load.global;

import daxing.common.ColumnTableTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.load.ancestralSite.Standardization;
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

    public static void globalPopulation(String inputDir, String vmapIIGroupFile, String outDir, String nonsynOrDel){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=vmapIIGroupTable.getHashMap(0,15);
        Predicate<File> cultivarPredict= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Cultivar");
        Predicate<File> landraceEUPredict= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Landrace_Europe");
        List<File> cultivarFiles=files.stream().filter(cultivarPredict).collect(Collectors.toList());
        List<File> landraceEUFiles=files.stream().filter(landraceEUPredict).collect(Collectors.toList());
        List<ColumnTableTool<String>> cultivarTables=new ArrayList<>();
        List<ColumnTableTool<String>> landraceEUTables=new ArrayList<>();
        RowTableTool<String> rowTableTool;
        ColumnTableTool<String> columnTable;
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
            columnTable=new ColumnTableTool<>(PStringUtils.fastSplit(header), cells);
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
            columnTable=new ColumnTableTool<>(PStringUtils.fastSplit(header), cells);
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
        List<ColumnTableTool<String>> globalColumn=new ArrayList<>();
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
        RowTableTool<String> globalAlledDistance=addDistance(global);
        globalAlledDistance.write(new File(outDir, nonsynOrDel+"Global.txt.gz"), IOFileFormat.TextGzip);
        new File(outDir, "cultivar.txt.gz").delete();
        new File(outDir, "landraceEU.txt.gz").delete();
    }

    private static RowTableTool<String> average(List<ColumnTableTool<String>> tableList){
        List<String> a,b,d;
        String averageA, averageB, averageD, region;
        List<List<String>> cells=new ArrayList<>();
        List<String> cell;
        double[] totalABD;
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

    private static RowTableTool<String> addDistance(RowTableTool<String> global){
        double[] abdGlobal, abd;
        double distance;
        List<List<String>> cells=new ArrayList<>();
        List<String> newHeader=new ArrayList<>();
        List<String> header=global.getHeader();
        StringBuilder headerSB=new StringBuilder();
        headerSB.append(header.get(0)).append("\t").append(header.get(1)).append("\t");
        headerSB.append(header.get(2)).append("\t").append(header.get(3)).append("\t");
        headerSB.append(header.get(4)).append("\t");
        String taxonName;
        for (int i = 5; i < header.size(); i=i+4) {
            headerSB.append(header.get(i)).append("\t");
            headerSB.append(header.get(i+1)).append("\t");
            headerSB.append(header.get(i+2)).append("\t");
            headerSB.append(header.get(i+3)).append("\t");
            taxonName=header.get(i+3).substring(0, header.get(i+3).length()-1);
            headerSB.append(taxonName+"Distance").append("\t");
        }
        headerSB.deleteCharAt(headerSB.length()-1);
        newHeader=PStringUtils.fastSplit(headerSB.toString());
        List<String> line;
        StringBuilder sb;
        for (int i = 0; i < global.getRowNumber(); i++) {
            sb=new StringBuilder(1000);
            sb.append(global.getCell(i, 0)).append("\t");
            sb.append(global.getCell(i,1)).append("\t");
            sb.append(global.getCell(i, 2)).append("\t");
            sb.append(global.getCell(i, 3)).append("\t");
            sb.append(global.getCell(i, 4)).append("\t");
            if (global.getCell(i, 4).equals("NA") || global.getCell(i, 4).equals("M000")){
                for (int j = 5; j < global.getColumnNumber(); j=j+4) {
                    sb.append(global.getCell(i, j)).append("\t");
                    sb.append(global.getCell(i, j+1)).append("\t");
                    sb.append(global.getCell(i, j+2)).append("\t");
                    sb.append(global.getCell(i, j+3)).append("\t");
                    sb.append("NA").append("\t");
                }
            }else {
                abdGlobal=new double[3];
                abd=new double[3];
                abdGlobal[0]=Double.parseDouble(global.getCell(i, 1));
                abdGlobal[1]=Double.parseDouble(global.getCell(i, 2));
                abdGlobal[2]=Double.parseDouble(global.getCell(i, 3));
                for (int j = 5; j < global.getColumnNumber(); j=j+4) {
                    sb.append(global.getCell(i, j)).append("\t");
                    sb.append(global.getCell(i, j+1)).append("\t");
                    sb.append(global.getCell(i, j+2)).append("\t");
                    sb.append(global.getCell(i, j+3)).append("\t");
                    if (global.getCell(i, j+3).equals("NA") || global.getCell(i, j+3).equals("M000")){
                        sb.append("NA").append("\t");
                    }else {
                        abd[0]=Double.parseDouble(global.getCell(i, j));
                        abd[1]=Double.parseDouble(global.getCell(i, j+1));
                        abd[2]=Double.parseDouble(global.getCell(i, j+2));
                        distance=Standardization.getDistance(abdGlobal, abd);
                        sb.append(distance).append("\t");
                    }
                }
            }
            sb.deleteCharAt(sb.length()-1);
            line=PStringUtils.fastSplit(sb.toString());
            cells.add(line);
        }
        return new RowTableTool<>(newHeader, cells);
    }

    public static void globalIndividual(String inputDir, String outDir, String nonsynOrDel){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        List<ColumnTableTool<String>> hexaploidTables=new ArrayList<>();
        List<String> taxonNames=new ArrayList<>();
        RowTableTool<String> rowTableTool;
        ColumnTableTool<String> columnTable;
        List<String> triadID, loadA, loadB, loadD, region;
        double[] loadABD;
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
        String taxon;
        for (int i = 0; i < files.size(); i++) {
            taxon=PStringUtils.fastSplit(files.get(i).getName(), "triad.").get(0);
            taxonNames.add(taxon);
            rowTableTool=new RowTableTool<>(files.get(i).getAbsolutePath());
            triadID=rowTableTool.getColumn(0);
            loadA=rowTableTool.getColumn(index[0]);
            loadB=rowTableTool.getColumn(index[1]);
            loadD=rowTableTool.getColumn(index[2]);
            region=new ArrayList<>();
            for (int j = 0; j < loadA.size(); j++) {
                if (loadA.get(j).equals("NA") || loadB.get(j).equals("NA") || loadD.get(j).equals("NA")){
                    region.add("NA");
                }else {
                    loadABD=new double[3];
                    loadABD[0]=Double.parseDouble(loadA.get(j));
                    loadABD[1]=Double.parseDouble(loadB.get(j));
                    loadABD[2]=Double.parseDouble(loadD.get(j));
                    if (loadABD[0]==0 && loadABD[1]==0 && loadABD[2]==0){
                        region.add("M000");
                    }else {
                        region.add(Standardization.getNearestPointIndex(loadABD).getRegion());
                    }
                }
            }
            cells=new ArrayList<>();
            cells.add(triadID);
            cells.add(loadA);
            cells.add(loadB);
            cells.add(loadD);
            cells.add(region);
            columnTable=new ColumnTableTool<>(PStringUtils.fastSplit(header), cells);
            hexaploidTables.add(columnTable);
        }
        RowTableTool<String> global=average(hexaploidTables);
        List<List<String>> individualColumns;
        String individualHeader;
        List<String> individualHeaderList;
        for (int i = 0; i < hexaploidTables.size(); i++) {
            individualHeader=getIndividualHeader(taxonNames.get(i));
            individualColumns=hexaploidTables.get(i).getCells();
            individualHeaderList=PStringUtils.fastSplit(individualHeader);
            for (int j = 1; j < individualColumns.size(); j++) {
                global.addColumn(individualHeaderList.get(j), individualColumns.get(j));
            }
        }
        RowTableTool<String> globalAddedDistance=addDistance(global);
        globalAddedDistance.write(new File(outDir, nonsynOrDel+"Global.txt.gz"), IOFileFormat.TextGzip);
    }

    private static String getIndividualHeader(String taxonName){
        String header="\tLoadA\tLoadB\tLoadD\tRegion";
        List<String> temp=PStringUtils.fastSplit(header);
        StringBuilder sb=new StringBuilder();
        for (int i = 0; i < temp.size(); i++) {
            sb.append(taxonName).append(temp.get(i)).append("\t");
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
