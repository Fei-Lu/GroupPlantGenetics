package daxing.v2.localAncestryInfer;

import daxing.common.utiles.IOTool;
import gnu.trove.list.TIntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Future;

public class LocalAncestryInferenceStart {

    public static void InferLocalAncestry(File genotypeFile, File groupInfoFile, List<File> fd_dxyFiles){
        GenotypeTable genoGrid = new GenotypeTable(genotypeFile.getAbsolutePath());
        Map<WindowSource.Source, List<String>> srcIndividualMap = getSrcPopMap(groupInfoFile.getAbsolutePath());
        List<String> queryList = getQuery(groupInfoFile.getAbsolutePath());
        List<Future<List<TIntList>>> callableList = new ArrayList<>();
    }

    public static int inferLocalAncestry(String refChr, GenotypeTable genoTable,
                                                    Map<WindowSource.Source, List<String>> srcIndividualMap,
                                                    Map<String, WindowSource.Source> taxaSourceMap,
                                                    IndividualSource individualSource,
                                         File outFile){
        WindowSource[] toBeInferredWindow = individualSource.selectCandidateWindow(2);
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
        String queryTaxon;
        int startIndex, endIndex, queryTaxonIndex;
        double[][] srcGenotype;
        double[] queryGenotype;
        List<TIntList> solution;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxon\tChr\tStart\tEnd\tSource");
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
                siteIndex[0]=startIndex;
                siteIndex[1]=endIndex;

                // query taxon index
                queryTaxon =  individualSource.getIndividualID();
                queryTaxonIndex = genoTable.getTaxonIndex(queryTaxon);

                srcGenotype = genoTable.getSrcGenotypeFrom(srcTaxaIndices, siteIndex);
                queryGenotype = genoTable.getQueryGenotypeFrom(queryTaxonIndex, siteIndex);
                solution = GenotypeTable.getMiniPath(srcGenotype, queryGenotype);

                // solution + srcIndiList + siteIndex + taxaSourceMap
                // final result

            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return 0;
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
                source = WindowSource.Source.valueOf(temp.get(1));
                if (srcPopSet.contains(source)){
                    srcIndividualMap.get(source).add(temp.get(0));
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return srcIndividualMap;
    }


}
