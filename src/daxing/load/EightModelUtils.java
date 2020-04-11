package daxing.load;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.Triad;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class EightModelUtils {

    public static void merge(String inputDir, String vmapIIGroupFile, String outDir){
        List<File> dirs=IOTool.getVisibleDir(inputDir);
        RowTableTool<String> table=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=table.getHashMap(0, 15);
        Predicate<File> p=f->taxonGroupMap.get(f.getName()).equals("Landrace_Europe") || taxonGroupMap.get(f.getName()).equals("Cultivar");
        List<File> ldCl=dirs.stream().filter(p).collect(Collectors.toList());
        List<File> temp;
        RowTableTool<String> tableTool;
        String outName;
        for (int i = 0; i < ldCl.size(); i++) {
            temp= IOUtils.getVisibleFileListInDir(ldCl.get(i).getAbsolutePath());
            tableTool=new RowTableTool<>(temp.get(0).getAbsolutePath());
            outName= PStringUtils.fastSplit(temp.get(0).getName(), ".").get(1);
            for (int j = 1; j < temp.size(); j++) {
                tableTool.add(new RowTableTool<>(temp.get(j).getAbsolutePath()));
            }
            tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
        }
    }

    public static void retainTriad(String triadFile, String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        IntStream.range(0, files.size()).parallel().forEach(e->retainTriad(triadFile, files.get(e), outDir));
    }

    private static void retainTriad(String triadFile, File inputFile, String outDir){
        Triad triad=new Triad(triadFile);
        RowTableTool<String> table=new RowTableTool<>(inputFile.getAbsolutePath());
        List<String> geneNames=table.getColumn(0);
        String[] threeName;
        int indexABD[]=new int[3];
        String triadID, taxon, temp;
        List<String> lineA, lineB, lineD;
        taxon=PStringUtils.fastSplit(inputFile.getName(), ".").get(0);
        StringBuilder sb;
        BufferedWriter bw;
        try {
            bw=IOTool.getTextGzipWriter(new File(outDir, taxon + ".triad.txt.gz"));
            bw.write("TriadID\tGeneName\tnumSyn\tnumDerivedInSyn\tnumNonsyn\tnumDerivedInNonsyn\tnumHGDeleterious" +
                    "\tnumDerivedInHGDeleterious");
            bw.newLine();
            for (int i = 0; i < triad.getRowNum(); i++) {
                threeName=triad.getTriad(i);
                for (int j = 0; j < indexABD.length; j++) {
                    indexABD[j]=Collections.binarySearch(geneNames, threeName[j]);
                }
                if (indexABD[0] < 0) continue;
                if (indexABD[1] < 0) continue;
                if (indexABD[2] < 0) continue;
                lineA=table.getRow(indexABD[0]);
                lineB=table.getRow(indexABD[1]);
                lineD=table.getRow(indexABD[2]);
                if (lineA.get(1).equals("-1")) continue;
                if (lineB.get(1).equals("-1")) continue;
                if (lineD.get(1).equals("-1")) continue;
                triadID=triad.getTraidID(i);
                for (int j = 0; j < indexABD.length; j++) {
                    temp=String.join("\t", table.getRow(indexABD[j]));
                    sb=new StringBuilder();
                    sb.append(triadID).append("\t");
                    sb.append(temp);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void countEightModel(String inputDir, String vmapIIGroupFile, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String, String> taxonGroupMap=vmapIIGroupTable.getHashMap(0, 15);
        IntStream.range(0, files.size()).parallel().forEach(e->countEightModel(files.get(e), taxonGroupMap, outDir));
    }

    private static void countEightModel(File inputFile, Map<String,String> taxonGroupMap, String outDir){
        String taxonName=PStringUtils.fastSplit(inputFile.getName(), ".").get(0);
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextGzipWriter(new File(outDir, taxonName+".triad.eightModel.txt.gz"))) {
            br.readLine();
            String line;
            String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
            Arrays.sort(models);
            String[] deleterious={"synDerivedCount","nonsynDerivedCount","delDerivedCount"};
            int[][] eightModel=new int[deleterious.length][models.length];
            List<String>[] temp;
            StringBuilder sbDel, sbNonsyn, sbSyn;
            int modelIndexSyn, modelIndexNonsyn, modelIndexDel;
            while ((line=br.readLine())!=null){
                temp=new List[3];
                temp[0]=PStringUtils.fastSplit(line);
                for (int i = 1; i < temp.length; i++) {
                    temp[i]=PStringUtils.fastSplit(br.readLine());
                }
                sbSyn=new StringBuilder();
                sbNonsyn=new StringBuilder();
                sbDel=new StringBuilder();
                sbSyn.append("M");
                sbNonsyn.append("M");
                sbDel.append("M");
                for (int i = 0; i < temp.length; i++) {
                    sbSyn.append(temp[i].get(3));
                    sbNonsyn.append(temp[i].get(5));
                    sbDel.append(temp[i].get(7));
                }
                modelIndexSyn=Arrays.binarySearch(models, sbSyn.toString());
                modelIndexNonsyn=Arrays.binarySearch(models, sbNonsyn.toString());
                modelIndexDel=Arrays.binarySearch(models, sbDel.toString());
                eightModel[0][modelIndexSyn]=eightModel[0][modelIndexSyn]+1;
                eightModel[1][modelIndexNonsyn]=eightModel[1][modelIndexNonsyn]+1;
                eightModel[2][modelIndexDel]=eightModel[2][modelIndexDel]+1;
            }
            bw.write("model\tsynDerivedCount\tnonsynDerivedCount\tdelDerivedCount\ttaxon\tgroup");
            bw.newLine();
            StringBuilder sb;
            for (int i = 0; i < models.length; i++) {
                sb=new StringBuilder();
                sb.append(models[i]).append("\t");
                for (int j = 0; j < eightModel.length; j++) {
                    sb.append(eightModel[j][i]).append("\t");
                }
                sb.append(taxonName).append("\t").append(taxonGroupMap.get(taxonName));
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void mergeModel(String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        RowTableTool<String> table=new RowTableTool<>(files.get(0).getAbsolutePath());
        for (int i = 1; i < files.size(); i++) {
            table.add(new RowTableTool<>(files.get(i).getAbsolutePath()));
        }
        table.write(new File(outDir, "modelMerged.txt.gz"), IOFileFormat.TextGzip);
    }
}
