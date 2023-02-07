package daxing.v2.localAncestryInfer.laidp;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.commons.lang3.EnumUtils;
import org.apache.commons.lang3.math.NumberUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class LocalAncestryInferenceStart {

    public static void inferStart(String genotypeDir, String fd_dxyFileDir, String groupInfoFile, int conjunctionNum,
                                  double initializeSwitchCostScore,
                                  String outDir, int maxSolutionCount){
        List<File> genotypeFile = IOTool.getFileListInDirEndsWith(genotypeDir, ".vcf.gz");
        String[] chrs = genotypeFile.stream().map(File::getName).map(s -> s.substring(3,5)).toArray(String[]::new);
        for (int i = 0; i < genotypeFile.size(); i++) {
            inferLocalAncestry(chrs[i], genotypeFile.get(i), new File(groupInfoFile), new File(fd_dxyFileDir),
                    conjunctionNum,initializeSwitchCostScore, new File(outDir), maxSolutionCount);
        }
    }

    public static void inferLocalAncestry(String refChr, File genotypeFile, File groupInfoFile, File fd_dxyFileDir,
                                          int conjunctionNum, double initializeSwitchCostScore,
                                          File outDir, int maxSolutionCount){
        GenotypeTable genoTable = new GenotypeTable(genotypeFile.getAbsolutePath());
        Map<WindowSource.Source, List<String>> srcIndividualMap = getSrcPopMap(groupInfoFile.getAbsolutePath());
        Map<String, WindowSource.Source> taxaSourceMap = getTaxaSourceMap(groupInfoFile.getAbsolutePath());
        Map<String, String> introgressionID2VcfIDMap = RowTableTool.getMap(groupInfoFile.getAbsolutePath(), 1,0);
        List<File> fd_dxyFiles = IOTool.getFileListInDirEndsWith(fd_dxyFileDir.getAbsolutePath(), ".txt.gz");
        boolean containSimulate = fd_dxyFiles.get(0).getName().contains("simulate");
        int taxonIndexLength = containSimulate ? 6 : 7;
        String[] queryTaxa= fd_dxyFiles.stream().map(File::getName).map(s->introgressionID2VcfIDMap.get(s.substring(13,
                13+taxonIndexLength))).toArray(String[]::new);
        String[] outFiles = Arrays.stream(queryTaxa).map(s -> "chr"+refChr+"_"+s+"_LAI.txt").toArray(String[]::new);

//        List<Callable<Integer>> callableTasks = new ArrayList<>();

        for (int i = 0; i < fd_dxyFiles.size(); i++) {
            int k = i;

//            callableTasks.add(()-> inferLocalAncestry3(refChr, genoTable, srcIndividualMap, taxaSourceMap,
//                    fd_dxyFiles.get(k), queryTaxa[k], new File(outDir, outFiles[k]), conjunctionNum,
//                    initializeSwitchCostScore, maxSolutionCount));

            inferLocalAncestry(refChr, genoTable, srcIndividualMap, taxaSourceMap,
                    fd_dxyFiles.get(k), queryTaxa[k], new File(outDir, outFiles[k]), conjunctionNum,
                    initializeSwitchCostScore, maxSolutionCount);

        }
//        ExecutorService executorService = Executors.newFixedThreadPool(32);
//        List<Integer> exitCodes = new ArrayList<>();
//        long start = System.nanoTime();
//        try {
//            List<Future<Integer>> futureList=executorService.invokeAll(callableTasks);
//            executorService.shutdown();
//            executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
//            for (Future<Integer> future : futureList){
//                exitCodes.add(future.get());
//            }
//        } catch (InterruptedException | ExecutionException e) {
//            e.printStackTrace();
//        }
//        List<String> failCommandList= new ArrayList<>();
//        for (int i = 0; i < exitCodes.size(); i++) {
//            if (exitCodes.get(i)!=0){
//                failCommandList.add(fd_dxyFiles.get(i).getName()+ "_"+refChr+" command fail to completed");
//            }
//        }
//        if (failCommandList.size()==0){
//            System.out.println("all commands had completed in "+ Benchmark.getTimeSpanHours(start)+ " hours");
//        }else {
//            System.out.println(failCommandList.size()+ "commands run failed");
//            System.out.println("Total spend "+Benchmark.getTimeSpanHours(start)+ " hours");
//        }
    }

    public static int inferLocalAncestry(String refChr, GenotypeTable genoTable,
                                          Map<WindowSource.Source, List<String>> srcIndividualMap,
                                          Map<String, WindowSource.Source> taxaSourceMap,
                                          File fd_dxyFile, String queryTaxon,
                                          File outFile, int conjunctionNum, double switchCostScore, int maxSolutionCount){
        IndividualSource queryIndividualSource = new IndividualSource(fd_dxyFile.getAbsolutePath(), queryTaxon);
        WindowSource[] toBeInferredWindow = queryIndividualSource.selectCandidateWindow(conjunctionNum);
        List<WindowSource> toBeInferredWindowChrList = new ArrayList<>();
        for (WindowSource windowSource: toBeInferredWindow){
            if (windowSource.getChrRange().getChr().equals(refChr)){
                toBeInferredWindowChrList.add(windowSource);
            }
        }
        EnumSet<WindowSource.Source> sourceEnumSet;
        List<String> srcIndiList;
        int[] srcTaxaIndices;
        int[] siteIndex;
        int startIndex, endIndex, queryTaxonIndex;
        double[][] srcGenotype;
        double[] queryGenotype;
        int startPos, endPos;
        IntList solution;
        EnumMap<Direction, IntList[]> biDirectionCandidateSolution;
        EnumSet<WindowSource.Source> sources;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxon\tChr\tWindowID\tSource\tStart\tEnd");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < toBeInferredWindowChrList.size(); i++) {
                sourceEnumSet = toBeInferredWindowChrList.get(i).getSources();
                srcIndiList = new ArrayList<>();
                for (WindowSource.Source source: sourceEnumSet){
                    srcIndiList.addAll(srcIndividualMap.get(source));
                }

                // src taxa indices;
                srcTaxaIndices = new int[srcIndiList.size()];
                for (int j = 0; j < srcTaxaIndices.length; j++) {
                    srcTaxaIndices[j] = genoTable.getTaxonIndex(srcIndiList.get(j));
                }

                // site start and end index
                siteIndex = new int[2];
                startIndex=genoTable.getSiteIndex(refChr, toBeInferredWindowChrList.get(i).getChrRange().getStart());
                endIndex = genoTable.getSiteIndex(refChr, toBeInferredWindowChrList.get(i).getChrRange().getEnd()-1);
                if (startIndex < -1){
                    startIndex=-startIndex-2;
                }else if (startIndex == -1){
                    startIndex = 0;
                }
                if (endIndex < 0){
                    endIndex=-endIndex-2;
                }
                siteIndex[0]=startIndex;
                siteIndex[1]=endIndex;

                // query taxon index
                queryTaxonIndex = genoTable.getTaxonIndex(queryTaxon);

                // solution of mini distance
                srcGenotype = genoTable.getSrcGenotypeFrom(srcTaxaIndices, siteIndex);
                queryGenotype = genoTable.getQueryGenotypeFrom(queryTaxonIndex, siteIndex);

                biDirectionCandidateSolution = SolutionUtils.getBiDirectionCandidateSourceSolution(srcGenotype,
                        queryGenotype, switchCostScore, srcIndiList, taxaSourceMap, maxSolutionCount);
                solution =SolutionUtils.calculateBreakPoint(biDirectionCandidateSolution);

                // write
                if (solution.size()==0){
                    sb.setLength(0);
                    sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                    sb.append(toBeInferredWindowChrList.get(i).getChrRange().toString()).append("\t");
//                    sb.append("NONE").append("\t");
                    sb.append("4").append("\t");
                    startPos = -1;
                    endPos = -1;
                    sb.append(startPos).append("\t").append(endPos);
                    bw.write(sb.toString());
                    bw.newLine();
                }else {
                    for (int j = 0; j < solution.size(); j=j+3) {
                        sb.setLength(0);
                        sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                        sb.append(toBeInferredWindowChrList.get(i).getChrRange().toString()).append("\t");
                        sources = SourceType.getSourcesFrom(solution.getInt(j));
                        for (WindowSource.Source source1: sources){
//                            sb.append(source1.name()).append(",");
                            sb.append(source1.getIndex()).append(",");
                        }
                        sb.deleteCharAt(sb.length()-1);
                        sb.append("\t");
                        startPos = genoTable.getPosition(solution.getInt(j+1)+startIndex);
                        endPos = genoTable.getPosition(solution.getInt(j+2)+startIndex);
                        sb.append(startPos).append("\t").append(endPos);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return 0;
    }

    public static Map<WindowSource.Source, List<String>> getSrcPopMap(String groupInfoFile){
        Map<WindowSource.Source, List<String>> srcIndividualMap = new HashMap<>();
        EnumSet<WindowSource.Source> srcPopSet = EnumSet.of(WindowSource.Source.WE, WindowSource.Source.DE,
                WindowSource.Source.FTT, WindowSource.Source.AT, WindowSource.Source.NONE);
        for (WindowSource.Source source: srcPopSet){
            srcIndividualMap.put(source, new ArrayList<>());
        }
        WindowSource.Source source;
        try (BufferedReader br = IOTool.getReader(groupInfoFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (NumberUtils.isNumber(temp.get(1))){
                    source = WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp.get(1))).get();
                    if (srcPopSet.contains(source)){
                        srcIndividualMap.get(source).add(temp.get(0));
                    }
                }else if(EnumUtils.isValidEnum(WindowSource.Source.class, temp.get(1))){
                    source = WindowSource.Source.valueOf(temp.get(1));
                    if (srcPopSet.contains(source)){
                        srcIndividualMap.get(source).add(temp.get(0));
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return srcIndividualMap;
    }

    public static Map<String, WindowSource.Source> getTaxaSourceMap(String groupInfoFile){
        Map<String, WindowSource.Source> taxaSourceMap = new HashMap<>();
        EnumSet<WindowSource.Source> srcPopSet = EnumSet.of(WindowSource.Source.WE, WindowSource.Source.DE,
                WindowSource.Source.FTT, WindowSource.Source.AT, WindowSource.Source.NONE);
        WindowSource.Source source;
        try (BufferedReader br = IOTool.getReader(groupInfoFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (NumberUtils.isNumber(temp.get(1))){
                    source = WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp.get(1))).get();
                    if (srcPopSet.contains(source)){
                        taxaSourceMap.put(temp.get(0), source);
                    }
                }
                else if(EnumUtils.isValidEnum(WindowSource.Source.class, temp.get(1))){
                    source = WindowSource.Source.valueOf(temp.get(1));
                    if (srcPopSet.contains(source)){
                        taxaSourceMap.put(temp.get(0), source);
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return taxaSourceMap;
    }


}
