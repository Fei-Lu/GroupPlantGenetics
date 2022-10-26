package daxing.v2.localAncestryInfer;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import org.apache.commons.lang3.EnumUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class LocalAncestryInferenceStart {

    public static void InferLocalAncestry(String refChr, File genotypeFile, File groupInfoFile, File fd_dxyFileDir,
                                          int conjunctionNum, double initializeSwitchCostScore,
                                          File outDir, int maxSolutionCount){
        GenotypeTable genoTable = new GenotypeTable(genotypeFile.getAbsolutePath());
        Map<WindowSource.Source, List<String>> srcIndividualMap = getSrcPopMap(groupInfoFile.getAbsolutePath());
        Map<String, WindowSource.Source> taxaSourceMap = getTaxaSourceMap(groupInfoFile.getAbsolutePath());
        Map<String, String> introgressionID2VcfIDMap = RowTableTool.getMap(groupInfoFile.getAbsolutePath(), 1,0);
        List<File> fd_dxyFiles = IOTool.getFileListInDirEndsWith(fd_dxyFileDir.getAbsolutePath(), ".txt");
        String[] queryTaxa= fd_dxyFiles.stream().map(File::getName).map(s->introgressionID2VcfIDMap.get(s.substring(13,
                20))).toArray(String[]::new);
        String[] outFiles = Arrays.stream(queryTaxa).map(s -> "chr"+refChr+"_"+s+"_LAI.txt").toArray(String[]::new);

        List<Callable<Integer>> callableTasks = new ArrayList<>();
        for (int i = 0; i < fd_dxyFiles.size(); i++) {
            int k = i;
//            callableTasks.add(()-> inferLocalAncestry(refChr, genoTable, srcIndividualMap, taxaSourceMap,
//                    fd_dxyFiles.get(k), queryTaxa[k], new File(outDir, outFiles[k])));

            inferLocalAncestry(refChr, genoTable, srcIndividualMap, taxaSourceMap,
                    fd_dxyFiles.get(k), queryTaxa[k], new File(outDir, outFiles[k]), conjunctionNum,
                    initializeSwitchCostScore, maxSolutionCount);

//            inferLAI(refChr, genoTable, srcIndividualMap, taxaSourceMap, queryTaxa[k], new File(outDir, outFiles[k]),
//                    initializeSwitchCostScore, maxSolutionCount);
        }
//        ExecutorService executorService = Executors.newFixedThreadPool(1);
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
        long start0 =System.nanoTime();
        IndividualSource queryIndividualSource = new IndividualSource(fd_dxyFile.getAbsolutePath(), queryTaxon);
        System.out.println();
        System.out.println("Current taxon and chr: "+queryIndividualSource.getIndividualID()+" "+refChr);
        System.out.println();
        System.out.println("Current option");
        System.out.println("conjunctionNum: "+ conjunctionNum);
        System.out.println("switchCostScore: "+ switchCostScore);
//        System.out.println("maxSolutionCount: "+maxSolutionCount);
//        System.out.println("maxSwitchCostScore: "+maxSwitchCostScore);
        System.out.println();
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
//        List<List<SolutionElement>> solution;
        List<int[]> solution;
        int startPos, endPos;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxon\tChr\tWindowID\tSource\tStart\tEnd\tSolution");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            StringBuilder log = new StringBuilder();
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
                if (startIndex < 0){
                    startIndex=-startIndex-2;
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

                long start = System.nanoTime();
                System.out.println();
                System.out.println("********* Start iteration *********");
                System.out.println(toBeInferredWindowChrList.get(i).getChrRange().toString());
                log.setLength(0);
                log.append("Window sources: ");
                log.append(String.join("\t",sourceEnumSet.stream().map(source -> source.name()).collect(Collectors.toList())));
                System.out.println(log);

                solution = GenotypeTable.getMiniPath22(srcGenotype, queryGenotype,
                switchCostScore, srcIndiList, taxaSourceMap, maxSolutionCount);

                System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds");
                System.out.println("********* End iteration *********");
                System.out.println();
                // write

                WindowSource.Source source;
                if (solution==null){
                    sb.setLength(0);
                    sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                    sb.append(toBeInferredWindowChrList.get(i).getChrRange().toString()).append("\t");
                    sb.append("NONE").append("\t");
                    startPos = -1;
                    endPos = -1;
                    sb.append(startPos).append("\t").append(endPos).append("\t");
                    sb.append(-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }else {
                    for (int j = 0; j < solution.size(); j++) {
                        for (int k = 0; k < solution.get(j).length; k=k+3) {
                            if (solution.get(j)[k]==-1) continue;
                            sb.setLength(0);
                            sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                            sb.append(toBeInferredWindowChrList.get(i).getChrRange().toString()).append("\t");
                            source = WindowSource.Source.getInstanceFromSubNum(solution.get(j)[k]).get();
                            sb.append(source.name()).append("\t");
                            startPos = genoTable.getPosition(solution.get(j)[k+1]+startIndex);
                            endPos = genoTable.getPosition((solution.get(j)[k+2]-1+startIndex));
                            sb.append(startPos).append("\t").append(endPos).append("\t");
                            sb.append(j);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(queryTaxon+" "+refChr+" completed in "+Benchmark.getTimeSpanMinutes(start0)+ " minutes");
        return 0;
    }

    public static void inferLAI(String refChr, GenotypeTable genoTable,
                                Map<WindowSource.Source, List<String>> srcIndividualMap,
                                Map<String, WindowSource.Source> taxaSourceMap, String queryTaxon,
                                File outFile, double switchCostScore, int maxSolutionCount){
        System.out.println();


        int[] siteIndex;
        int startIndex, endIndex, queryTaxonIndex;
        double[][] srcGenotype;
        double[] queryGenotype;
//        BitSet[] srcGenotype;
//        BitSet queryGenotype;
        int startPos;
        int endPos;
//        List<List<SolutionElement>> solution;
        List<int[]> solution;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxon\tChr\tSource\tStart\tEnd\tSolution");
            bw.newLine();
            StringBuilder sb = new StringBuilder();
            WindowSource.Source currentHaplotype;
            int firstPos;
            StringBuilder log = new StringBuilder();
            EnumSet<WindowSource.Source> sourceEnumSet = WindowSource.Source.getABSource();
            List<String> srcIndiList = new ArrayList<>();
            for (WindowSource.Source source: sourceEnumSet){
                srcIndiList.addAll(srcIndividualMap.get(source));
            }

            // src taxa indices;
            int[] srcTaxaIndices = new int[srcIndiList.size()];
            for (int j = 0; j < srcTaxaIndices.length; j++) {
                srcTaxaIndices[j] = genoTable.getTaxonIndex(srcIndiList.get(j));
            }

            // site start and end index
            siteIndex = new int[2];
            startIndex=0;
            endIndex = genoTable.getSiteNumber()-1;
            siteIndex[0]=startIndex;
            siteIndex[1]=endIndex;

            // query taxon index
            queryTaxonIndex = genoTable.getTaxonIndex(queryTaxon);

            // solution of mini distance
            srcGenotype = genoTable.getSrcGenotypeFrom(srcTaxaIndices, siteIndex);
            queryGenotype = genoTable.getQueryGenotypeFrom(queryTaxonIndex, siteIndex);
//              srcGenotype = genoTable.getSrcGenotypeBitSetFrom(srcTaxaIndices, siteIndex);
//              queryGenotype = genoTable.getQueryGenotypeBitSetFrom(queryTaxonIndex, siteIndex);
            System.out.println();
            System.out.println("********* Start iteration *********");
            System.out.println();
            log.setLength(0);
            log.append("Window sources: ");
            log.append(String.join("\t",sourceEnumSet.stream().map(source -> source.name()).collect(Collectors.toList())));
            System.out.println(log);
            solution = GenotypeTable.getMiniPath22(srcGenotype, queryGenotype,
                    switchCostScore, srcIndiList, taxaSourceMap, maxSolutionCount);

            System.out.println("********* End iteration *********");
            System.out.println();
            // write
            if (solution==null){
                sb.setLength(0);
                sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                sb.append("NONE").append("\t");
                startPos = -1;
                endPos = -1;
                sb.append(startPos).append("\t").append(endPos).append("\t");
                sb.append(-1);
                bw.write(sb.toString());
                bw.newLine();
            }else {
                WindowSource.Source source;
                for (int j = 0; j < solution.size(); j++) {
                    for (int k = 0; k < solution.get(j).length; k=k+3) {
                        if (solution.get(j)[k]==-1) continue;
                        sb.setLength(0);
                        sb.append(queryTaxon).append("\t").append(refChr).append("\t");
                        source = WindowSource.Source.getInstanceFromSubNum(solution.get(j)[k]).get();
                        sb.append(source.name()).append("\t");
                        startPos = genoTable.getPosition(solution.get(j)[k+1]+startIndex);
                        endPos = genoTable.getPosition((solution.get(j)[k+2]-1+startIndex));
                        sb.append(startPos).append("\t").append(endPos).append("\t");
                        sb.append(j);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<String> getQuery(String groupInfoFile){
        List<String> res = new ArrayList<>();
        String line;
        List<String> temp;
        try (BufferedReader br = IOTool.getReader(groupInfoFile)) {
            while ((line = br.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                if (temp.get(1).startsWith("LR") || temp.get(1).startsWith("CL")){
                    res.add(temp.get(0));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
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
                if(EnumUtils.isValidEnum(WindowSource.Source.class, temp.get(1))){
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
                if(EnumUtils.isValidEnum(WindowSource.Source.class, temp.get(1))){
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
