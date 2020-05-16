package daxing.load.neutralSiteLoad;

import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.load.ancestralSite.Standardization;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class Utils {

    public static void mergeIndividual(String inputDir,String vmapIIGroupFile, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=vmapIIGroupTable.getHashMap(0,15);
        Predicate<File> cultivarPredict=f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Cultivar");
        Predicate<File> landraceEUPredict= f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(),".").get(0)).equals("Landrace_Europe");
        List<File> cultivarFiles=files.stream().filter(cultivarPredict).collect(Collectors.toList());
        List<File> landraceEUFiles=files.stream().filter(landraceEUPredict).collect(Collectors.toList());
        List<RowTableTool<String>> cultivarTables=new ArrayList<>();
        List<RowTableTool<String>> landraceEUTables=new ArrayList<>();
        for (int i = 0; i < cultivarFiles.size(); i++) {
            cultivarTables.add(exactLoad(cultivarFiles.get(i)));
        }
        for (int i = 0; i < landraceEUFiles.size(); i++) {
            landraceEUTables.add(exactLoad(landraceEUFiles.get(i)));
        }
        RowTableTool<String> cultivarAveraged=average(cultivarTables);
        RowTableTool<String> landraceEUAveraged=average(landraceEUTables);
        List<RowTableTool<String>> cultivarLandraceEUTables=new ArrayList<>();
        cultivarLandraceEUTables.add(cultivarAveraged);
        cultivarLandraceEUTables.add(landraceEUAveraged);
        String[] subDirs={"001_tissue","002_global"};
        for (int i = 0; i < subDirs.length; i++) {
            new File(outDir, subDirs[i]).mkdir();
        }
        List<String> cultivarRegion=getRegion(cultivarAveraged);
        List<String> landraceEURegion=getRegion(landraceEUAveraged);
        cultivarAveraged.addColumn("CultivarRegion", cultivarRegion);
        landraceEUAveraged.addColumn("LandraceEURegion", landraceEURegion);
        List<File> files1= IOTool.getVisibleDir(outDir);
        cultivarAveraged.write(new File(files1.get(0), "cultivar.nonsynLoad.txt.gz"), IOFileFormat.TextGzip);
        landraceEUAveraged.write(new File(files1.get(0), "landraceEU.nonsynLoad.txt.gz"), IOFileFormat.TextGzip);
        List<File> files2=IOTool.getVisibleDir(new File(outDir, subDirs[0]).getAbsolutePath());
        merge(files2, new File(outDir, subDirs[1]).getAbsolutePath());
    }

    private static void merge(List<File> files,String outDir){
        RowTableTool<String> cultivar=new RowTableTool<>(files.get(0).getAbsolutePath());
        RowTableTool<String> landrace=new RowTableTool<>(files.get(1).getAbsolutePath());
        List<RowTableTool<String>> cultivarLandraceTable=new ArrayList<>();
        cultivarLandraceTable.add(cultivar);
        cultivarLandraceTable.add(landrace);
        RowTableTool<String> global=average(cultivarLandraceTable);
        List<String> regionGlobal=getRegion(global);
        global.addColumn("GlobalRegion", regionGlobal);
        String[] cultivarNames={"CultivarLoadA","CultivarLoadB","CultivarLoadD","CultivarRegion"};
        String[] landraceEUNames={"Landrace_EULoadA","Landrace_EULoadB","Landrace_EULoadD","Landrace_EURegion"};
        addColumn(global, cultivar, cultivarNames);
        addColumn(global, landrace, landraceEUNames);
//        for (int i = 0; i < files.size(); i++) {
//            files.get(i).delete();
//        }
        global.write(new File(outDir,"global.txt.gz"), IOFileFormat.TextGzip);
    }

    private static RowTableTool<String> exactLoad(File individualTriadFile){
        RowTableTool<String> table=new RowTableTool<>(individualTriadFile.getAbsolutePath());
        List<String> triadID=table.getColumn(0);
        int[] cdsLenArray=table.getColumnAsIntArray("cdsLen");
        int[] numDerivedInNonsynArray=table.getColumnAsIntArray("numDerivedInNonsyn");
        int[] numDerivedInGeneLocalArray=table.getColumnAsIntArray("numDerivedInGeneLocal");
        int[] numDerivedInHGDeleteriousArray=table.getColumnAsIntArray("numDerivedInHGDeleterious");
        int cdsLen, numDerivedInNonsyn, numDerivedInGeneLocal, numDerivedInHGDeleterious;
        String[] loadNonsynArray=new String[cdsLenArray.length];
        double load;
        for (int i = 0; i < cdsLenArray.length; i++) {
            cdsLen=cdsLenArray[i];
            numDerivedInNonsyn=numDerivedInNonsynArray[i];
            numDerivedInGeneLocal=numDerivedInGeneLocalArray[i];
            numDerivedInHGDeleterious=numDerivedInHGDeleteriousArray[i];
            if (numDerivedInNonsyn==0 || (numDerivedInGeneLocal-numDerivedInNonsyn)==0){
                loadNonsynArray[i]="NA";
            }else {
                load=1000*100*((double)numDerivedInHGDeleterious)/(cdsLen*(numDerivedInGeneLocal-numDerivedInNonsyn));
                loadNonsynArray[i]= String.valueOf(NumberTool.format(load, 5));
            }
        }
        List<List<String>> cells=new ArrayList<>();
        List<String> triadIDABD;
        for (int i = 0; i < loadNonsynArray.length; i=i+3) {
            triadIDABD=new ArrayList<>();
            triadIDABD.add(triadID.get(i));
            triadIDABD.add(loadNonsynArray[i]);
            triadIDABD.add(loadNonsynArray[i+1]);
            triadIDABD.add(loadNonsynArray[i+2]);
            cells.add(triadIDABD);
        }
        String str="TriadID\tGlobalLoadA\tGlobalLoadB\tGlobalLoadD";
        List<String> header= PStringUtils.fastSplit(str);
        return new RowTableTool<>(header, cells);
    }

    private static RowTableTool<String> average(List<RowTableTool<String>> tableList){
        List<String> a,b,d;
        String averageA, averageB, averageD;
        List<List<String>> cells=new ArrayList<>();
        List<String> cell;
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
            cell=new ArrayList<>();
            cell.add(tableList.get(0).getCell(i, 0));
            cell.add(averageA);
            cell.add(averageB);
            cell.add(averageD);
            cells.add(cell);
        }
        String header="TriadID\tA\tB\tD";
        return new RowTableTool<>(PStringUtils.fastSplit(header), cells);
    }

    private static List<String> getRegion(RowTableTool<String> table){
        List<String> region=new ArrayList<>();
        double[] abd;
        String aStr, bStr, dStr;
        for (int i = 0; i < table.getRowNumber(); i++) {
            aStr=table.getCell(i, 1);
            bStr=table.getCell(i, 2);
            dStr=table.getCell(i, 3);
            if (aStr.equals("NA") || bStr.equals("NA") || dStr.equals("NA")){
                region.add("NA");
                continue;
            }
            abd=new double[3];
            abd[0]=Double.parseDouble(aStr);
            abd[1]=Double.parseDouble(bStr);
            abd[2]=Double.parseDouble(dStr);
            region.add(Standardization.getNearestPointIndex(abd).getRegion());
        }
        return region;
    }

    private static RowTableTool<String> addColumn(RowTableTool<String> table1, RowTableTool<String> table2,
                                                  String[] columnName){
        List<List<String>> columnList=new ArrayList<>();
        for (int i = 1; i < table2.getColumnNumber(); i++) {
            columnList.add(table2.getColumn(i));
        }
        for (int i = 0; i < columnList.size(); i++) {
            table1.addColumn(columnName[i],columnList.get(i));
        }
        return table1;
    }

    public static void calculateDistance(String inputFile, String outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextGzipWriter(outFile)) {
            String header=br.readLine();
            List<String> newHeader=PStringUtils.fastSplit(header);
            newHeader.add(9, "CultivarDist");
            newHeader.add("Landrace_EUDistance");
            bw.write(String.join("\t", newHeader));
            bw.newLine();
            String line;
            List<String> temp;
            double[] globalABD, landraceABD, cultivarABD;
            double cultivarDistance, landraceDistance;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                if (temp.contains("NA") || temp.contains("M000")){
                    temp.add(9, "NA");
                    temp.add("NA");
                    bw.write(String.join("\t",temp));
                    bw.newLine();
                    continue;
                }
                globalABD=new double[3];
                cultivarABD=new double[3];
                landraceABD=new double[3];
                globalABD[0]=Double.parseDouble(temp.get(1));
                globalABD[1]=Double.parseDouble(temp.get(2));
                globalABD[2]=Double.parseDouble(temp.get(3));
                cultivarABD[0]=Double.parseDouble(temp.get(5));
                cultivarABD[1]=Double.parseDouble(temp.get(6));
                cultivarABD[2]=Double.parseDouble(temp.get(7));
                landraceABD[0]=Double.parseDouble(temp.get(9));
                landraceABD[1]=Double.parseDouble(temp.get(10));
                landraceABD[2]=Double.parseDouble(temp.get(11));
                cultivarDistance=Standardization.getDistance(globalABD, cultivarABD);
                landraceDistance=Standardization.getDistance(globalABD, landraceABD);
                temp.add(9, String.valueOf(NumberTool.format(cultivarDistance,5)));
                temp.add(String.valueOf(NumberTool.format(landraceDistance,5)));
                bw.write(String.join("\t", temp));
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
