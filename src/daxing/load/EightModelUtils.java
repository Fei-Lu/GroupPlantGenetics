package daxing.load;

import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.common.Triad;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
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
        taxon=PStringUtils.fastSplit(inputFile.getName(), ".").get(0);
        StringBuilder sb;
        BufferedWriter bw;
        int[] cdsLen;
        try {
            bw=IOTool.getTextGzipWriter(new File(outDir, taxon + ".triad.txt.gz"));
            bw.write("TriadID\tcdsLen\tGeneName\tnumSyn\tnumDerivedInSyn\tnumNonsyn\tnumDerivedInNonsyn" +
                    "\tnumHGDeleterious\tnumDerivedInHGDeleterious\tsubgenome");
            bw.newLine();
            String subgenome, geneName;
            for (int i = 0; i < triad.getRowNum(); i++) {
                threeName=triad.getTriad(i);
                for (int j = 0; j < indexABD.length; j++) {
                    indexABD[j]=Collections.binarySearch(geneNames, threeName[j]);
                }
                if (indexABD[0] < 0) continue;
                if (indexABD[1] < 0) continue;
                if (indexABD[2] < 0) continue;
                triadID=triad.getTraidID(i);
                cdsLen=triad.getCDSLen(i);
                for (int j = 0; j < indexABD.length; j++) {
                    temp=String.join("\t", table.getRow(indexABD[j]));
                    sb=new StringBuilder();
                    sb.append(triadID).append("\t").append(cdsLen[j]).append("\t");
                    sb.append(temp).append("\t");
                    geneName=table.getRow(indexABD[j]).get(0);
                    subgenome=geneName.substring(8,9);
                    sb.append(subgenome);
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

    public static void filter(String triadFileInputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadFileInputDir);
        files.stream().parallel().forEach(e->filter(e, outDir));
    }

    private static void filter(File triadFile, String outDir){
        String fileName=triadFile.getName();
        try (BufferedReader br = IOTool.getReader(triadFile);
             BufferedWriter bw=IOTool.getTextGzipWriter(new File(outDir, fileName))) {
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            String[] lineABD;
            int[] snpNum;
            List<String> temp;
            int triadNum=0;
            int retainedNum=0;
            while ((line=br.readLine())!=null){
                lineABD=new String[3];
                lineABD[0]=line;
                lineABD[1]=br.readLine();
                lineABD[2]=br.readLine();
                triadNum++;
                snpNum=new int[3];
                for (int i = 0; i < lineABD.length; i++) {
                    temp=PStringUtils.fastSplit(lineABD[i]);
                    snpNum[i]=Integer.parseInt(temp.get(3))+Integer.parseInt(temp.get(5));
                }
                if (snpNum[0]==0) continue;
                if (snpNum[1]==0) continue;
                if (snpNum[2]==0) continue;
                retainedNum++;
                for (int i = 0; i < lineABD.length; i++) {
                    bw.write(lineABD[i]);
                    bw.newLine();
                }
            }
            System.out.println(triadFile.getName()+": "+triadNum+" "+retainedNum+" "+(double)retainedNum/triadNum);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void countEightModel(String triadInputDir, String vmapIIGroupFile, int derivedCountThresh, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadInputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String, String> taxonGroupMap=vmapIIGroupTable.getHashMap(0, 15);
        IntStream.range(0, files.size()).parallel().forEach(e->countEightModel(files.get(e), taxonGroupMap,
                derivedCountThresh, outDir));
    }

    public static void countEightModel(File triadInputFile, Map<String,String> taxonGroupMap, int derivedCountThresh,
                                       String outDir){
        String taxonName=PStringUtils.fastSplit(triadInputFile.getName(), ".").get(0);
        try (BufferedReader br = IOTool.getReader(triadInputFile);
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
            double derivedCount;
            int cdsLen, snpNum;
            double geneLoad;
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
                    derivedCount=Double.parseDouble(temp[i].get(4));
                    cdsLen=Integer.parseInt(temp[i].get(1));
                    snpNum=Integer.parseInt(temp[i].get(3))+Integer.parseInt(temp[i].get(5));
                    geneLoad=(1000*10*derivedCount)/(cdsLen*snpNum);
                    if (geneLoad < derivedCountThresh) {
                        sbSyn.append(0);
                    } else {
                        sbSyn.append(1);
                    }
                    derivedCount=Double.parseDouble(temp[i].get(6));
                    geneLoad=(1000*10*derivedCount)/(cdsLen*snpNum);
                    if (geneLoad < derivedCountThresh) {
                        sbNonsyn.append(0);
                    } else {
                        sbNonsyn.append(1);
                    }
                    derivedCount=Double.parseDouble(temp[i].get(8));
                    geneLoad=(1000*10*derivedCount)/(cdsLen*snpNum);
                    if (geneLoad < derivedCountThresh){
                        sbDel.append(0);
                    }else {
                        sbDel.append(1);
                    }
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

    public static void countEightModel(String triadNormalizedInputDir, String vmapIIGroupFile, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadNormalizedInputDir);
        RowTableTool<String> vmapIIGroupTable=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> map=vmapIIGroupTable.getHashMap(0,15);
        IntStream.range(0, files.size()).parallel().forEach(e->countEightModel(files.get(e),map, outDir));
    }

    public static void countEightModel(File triadNormalizedInputFile, Map<String,String> taxonGroupMap, String outDir){
        String outFileName=triadNormalizedInputFile.getName().replaceAll("\\.txt\\.gz$","eightModel.txt.gz");
        String taxonName=PStringUtils.fastSplit(triadNormalizedInputFile.getName(), ".").get(0);
        try (BufferedReader br = IOTool.getReader(triadNormalizedInputFile);
             BufferedWriter bw =IOTool.getTextGzipWriter(new File(outDir, outFileName))) {
            br.readLine();
            //Bluesky.triad.normalized.txt.gz
            //BeiJing8.triad.eightModel.txt.gz
            String line;
            List<String> temp;
            StringBuilder sb;
            String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
            Arrays.sort(models);
            int[][] modelCount=new int[3][8];
            int[] index;
            while((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                index=new int[3];
                index[0]=Arrays.binarySearch(models, temp.get(4));
                index[1]=Arrays.binarySearch(models, temp.get(8));
                index[2]=Arrays.binarySearch(models, temp.get(12));
                modelCount[0][index[0]]+=1;
                modelCount[1][index[1]]+=1;
                modelCount[2][index[2]]+=1;
            }
            bw.write("model\tsynDerivedCount\tnonsynDerivedCount\tdelDerivedCount\ttaxon\tgroup");
            bw.newLine();
            for (int i = 0; i < models.length; i++) {
                sb=new StringBuilder();
                sb.append(models[i]).append("\t");
                for (int j = 0; j < modelCount.length; j++) {
                    sb.append(modelCount[j][i]).append("\t");
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
        Predicate<File> isDir=File::isDirectory;
        List<File> files1=files.stream().filter(isDir.negate()).collect(Collectors.toList());
        RowTableTool<String> table=new RowTableTool<>(files1.get(0).getAbsolutePath());
        for (int i = 1; i < files1.size(); i++) {
            table.add(new RowTableTool<>(files1.get(i).getAbsolutePath()));
        }
        table.write(new File(outDir, "modelMerged.txt.gz"), IOFileFormat.TextGzip);
    }

    public static void normalizedTriad(String inputDir,String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        IntStream.range(0, files.size()).parallel().forEach(f->normalizedTriad(files.get(f), outDir));
    }

    public static void normalizedTriad(File inputFile, String outDir){
        String fileName=inputFile.getName().replaceAll("triad", "triad.normalized");
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextGzipWriter(new File(outDir, fileName))) {
            bw.write("TriadID\tnormalizedNumDerivedInSynA" +
                    "\tnormalizedNumDerivedInSynB\tnormalizedNumDerivedInSynD\tsynRegion" +
                    "\tnormalizedNumDerivedInNonsynA" +
                    "\tnormalizedNumDerivedInNonsynB\tnormalizedNumDerivedInNonsynD\tnonsynRegion" +
                    "\tnormalizedNumDerivedInHGDeleteriousA\tnormalizedNumDerivedInHGDeleteriousB" +
                    "\tnormalizedNumDerivedInHGDeleteriousD\tdelRegion");
            bw.newLine();
            br.readLine();
            String line, triadID, region;
            List<String>[] temp;
            int[] cdsLen;
            int[] snpNum;
            double[] derivedSyn;
            double[] derivedNonsyn;
            double[] deleterious;
            StringBuilder sb;
            double[] normalizedDerivedSyn;
            double[] normalizedDerivedNonSyn;
            double[] normalizedDerivedDel;
            while ((line=br.readLine())!=null){
                temp=new List[3];
                temp[0]=PStringUtils.fastSplit(line);
                temp[1]=PStringUtils.fastSplit(br.readLine());
                temp[2]=PStringUtils.fastSplit(br.readLine());
                cdsLen=new int[3];
                snpNum=new int[3];
                derivedSyn=new double[3];
                derivedNonsyn=new double[3];
                deleterious=new double[3];
                normalizedDerivedSyn=new double[3];
                normalizedDerivedNonSyn=new double[3];
                normalizedDerivedDel=new double[3];
                triadID=temp[0].get(0);
                for (int i = 0; i < cdsLen.length; i++) {
                    cdsLen[i]=Integer.parseInt(temp[i].get(1));
                    snpNum[i]=Integer.parseInt(temp[i].get(3))+Integer.parseInt(temp[i].get(5));
                    derivedSyn[i]=Double.parseDouble(temp[i].get(4));
                    derivedNonsyn[i]=Double.parseDouble(temp[i].get(6));
                    deleterious[i]=Double.parseDouble(temp[i].get(8));
                }
                for (int i = 0; i < 3; i++) {
                    normalizedDerivedSyn[i]=10000*derivedSyn[i]/(cdsLen[i]*snpNum[i]);
                    normalizedDerivedNonSyn[i]=10000*derivedNonsyn[i]/(cdsLen[i]*snpNum[i]);
                    normalizedDerivedDel[i]=10000*deleterious[i]/(cdsLen[i]*snpNum[i]);
                }
                sb=new StringBuilder();
                sb.append(triadID).append("\t");
                for (int i = 0; i < 3; i++) {
                    sb.append(NumberTool.format(normalizedDerivedSyn[i], 5)).append("\t");
                }
                region=Standardization.getNearestPointIndex(normalizedDerivedSyn).getRegion();
                sb.append(region).append("\t");
                for (int i = 0; i < 3; i++) {
                    sb.append(NumberTool.format(normalizedDerivedNonSyn[i], 5)).append("\t");
                }
                region=Standardization.getNearestPointIndex(normalizedDerivedNonSyn).getRegion();
                sb.append(region).append("\t");
                for (int i = 0; i < 3; i++) {
                    sb.append(NumberTool.format(normalizedDerivedDel[i], 5)).append("\t");
                }
                region=Standardization.getNearestPointIndex(normalizedDerivedDel).getRegion();
                sb.append(region);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
