package daxing;

import daxing.common.chrrange.ChrRange;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.*;
import gnu.trove.list.array.TIntArrayList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntSet;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Start {

    public static void main(String[] args) {

//        String genotypeDir = args[0];
//        String sample2PopFile_AB = args[1];
//        String sample2PopFile_D = args[2];
//        String recombinationMap = args[3];
//        String prunedInSNPFile_AB = args[4];
//        String prunedInSNPFile_D = args[5];
//        String alleleCountDir = args[6];


//        String genotypeDir = "/Users/xudaxing/Desktop/ancestryHMM/test/002_genotype";
//        String sample2PopInfoDir = "/Users/xudaxing/Desktop/ancestryHMM/test/003_sample2Pop";
//        String recombinationMapFile = "/Users/xudaxing/Desktop/ancestryHMM/test/001_geneticMap/iwgsc_refseqv1.0_mapping_data.txt";
//        String prunedInSNPFile_AB = "/Users/xudaxing/Desktop/ancestryHMM/test/004_prunedIn/chrAB_vmap2.1_onlyGenotype_haploid_panelAllUnion_pruned.in";
//        String prunedInSNPFile_D = "/Users/xudaxing/Desktop/ancestryHMM/test/004_prunedIn/chrD_vmap2.1_onlyGenotype_haploid_panelAllUnion_pruned.in";
//        String outDirPanel = "/Users/xudaxing/Desktop/ancestryHMM/test/005_out";
//
//        Panel.preparePanelForAncestryHMM(genotypeDir, sample2PopInfoDir, recombinationMapFile, prunedInSNPFile_AB,
//                prunedInSNPFile_D, outDirPanel);
////
//        String genotypeFile = "/Users/xudaxing/Desktop/ABBA/001_gentotype/chr2A_vmap2.1_onlyGenotype_haploid_imputation.vcf";
//        String fd_dxyFile = "/Users/xudaxing/Desktop/ABBA/002_dxy_fd";
//        String groupByPop2IndividualFile="/Users/xudaxing/Desktop/ABBA/groupByPop2Indi_indianDwarfToNONE_LRClose2IndianDwarf.txt";
//        String outDir="/Users/xudaxing/Desktop/ABBA/003_outDir";
//////
//////
//        int conjunctionNum=3; // 2
//        double switchCostScore=1.5; // 1.5
//        int maxSolutionCount=20; // 100
//        int threadNum= 10;
////
////        String genotypeFile = args[0];
////        String fd_dxyFile = args[1];
////        String groupByPop2IndividualFile=args[2];
////        String outDir=args[3];
////
////        int conjunctionNum=Integer.parseInt(args[4]); // 2
////        double switchCostScore=Double.parseDouble(args[5]); // 1.5
////        int maxSolutionCount=Integer.parseInt(args[6]); // 100
//////        int threadNum= Integer.parseInt(args[7]);
////
////
//        long start = System.nanoTime();
//        LocalAncestryInferenceStart.inferLocalAncestry("2A", new File(genotypeFile),
//                new File(groupByPop2IndividualFile), new File(fd_dxyFile), conjunctionNum, switchCostScore,
//                new File(outDir), maxSolutionCount);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds");

        String genotypeDir = args[0];
        String fd_dxyFileDir = args[1];
        String groupInfoFile = args[2];
        String outDir = args[3];
        int conjunctionNum = Integer.parseInt(args[4]);
        double initializeSwitchCostScore  =Double.parseDouble(args[5]);
        int maxSolutionCount = Integer.parseInt(args[6]);

//        String genotypeDir = "/Users/xudaxing/Desktop/ABBA/001_gentotype";
//        String fd_dxyFileDir = "/Users/xudaxing/Desktop/ABBA/002_dxy_fd";
//        String groupInfoFile = "/Users/xudaxing/Desktop/ABBA/groupByPop2Indi_indianDwarfToNONE_LRClose2IndianDwarf.txt";
//        String outDir = "/Users/xudaxing/Desktop/ABBA/003_outDir";
//        int conjunctionNum = 2;
//        double initializeSwitchCostScore  =1.5;
//        int maxSolutionCount = 32;
////        int threadNum = Integer.parseInt(args[7]);


        LocalAncestryInferenceStart.inferStart(genotypeDir, fd_dxyFileDir, groupInfoFile,conjunctionNum,
                initializeSwitchCostScore,outDir,maxSolutionCount);

//        List<SolutionElement> solutionElementList = new ArrayList<>();
//        solutionElementList.add(new SolutionElement(WindowSource.Source.WE, 0, 1));
//        solutionElementList.add(new SolutionElement(WindowSource.Source.DE, 2, 3));
//        solutionElementList.add(new SolutionElement(WindowSource.Source.WE, 4, 5));
//
//        List<SolutionElement> copy = new ArrayList<>();
//        Iterator<SolutionElement> iterator = solutionElementList.listIterator();
//        while (iterator.hasNext()){
//            copy.add(iterator.next().clone());
//        }
//        System.out.println();

//        double switchCostScore= 1.5;
//        double[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        double[] queryGenotype = {1,1,0,0,0,1,0,0,1,1,1,1};
//        double[][] miniCost = SolutionUtils.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);
//        IntSet[][] candidateSolution = SolutionUtils.getCandidateSolution2(miniCost, 1.5, );
//        System.out.println();

//        List<String> a = WheatLineage.abLineage();
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < a.size(); i++) {
//            sb.append(a.get(i)).append(" ");
//        }
//        System.out.println(sb.toString());

//        String genotypeDir = args[0];
//        String fd_dxyFileDir = args[1];
//        String groupInfoFile = args[2];
//        int conjunctionNum = Integer.parseInt(args[3]);
//        double initializeSwitchCostScore = Double.parseDouble(args[4]);
//        String outDir = args[5];
//        int maxSolutionCount = Integer.parseInt(args[6]);
//        LocalAncestryInferenceStart.inferStart(genotypeDir, fd_dxyFileDir, groupInfoFile, conjunctionNum,
//                initializeSwitchCostScore, outDir, maxSolutionCount);

    }

    public static void test(String introgressionMapFile, String infoSiteFile, String outFile){
        try (BufferedReader brInfoSite = IOTool.getReader(infoSiteFile);
             BufferedReader brMap = IOTool.getReader(introgressionMapFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;

            EnumMap<InfoSite, TIntArrayList> infoSite2PostionMap = new EnumMap<>(InfoSite.class);
            for (InfoSite infoSite: InfoSite.values()){
                infoSite2PostionMap.put(infoSite, new TIntArrayList());
            }
            InfoSite infoSite;
            brInfoSite.readLine();
            while ((line=brInfoSite.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                infoSite = InfoSite.valueOf(temp.get(2));
                infoSite2PostionMap.get(infoSite).add(Integer.parseInt(temp.get(1)));
            }
            for (Map.Entry<InfoSite, TIntArrayList> entry: infoSite2PostionMap.entrySet()){
                entry.getValue().sort();
            }

            int start, end, indexStart= -1, indexEnd = -1, hit, count;
            StringBuilder sb = new StringBuilder();
            TIntArrayList positionList;
            line=brMap.readLine();
            sb.setLength(0);
            sb.append(line).append("\t");
            for (int i = 0; i < InfoSite.values().length; i++) {
                sb.append(InfoSite.values()[i]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            while ((line=brMap.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                start = Integer.parseInt(temp.get(1));
                end = Integer.parseInt(temp.get(2));
                sb.setLength(0);
                sb.append(String.join("\t", temp)).append("\t");
                for (int i = 0; i < InfoSite.values().length; i++) {
                    positionList = infoSite2PostionMap.get(InfoSite.values()[i]);
                    hit = positionList.binarySearch(start);
                    if (hit == -1){
                        indexStart = 0;
                    }else if (hit < -1){
                        indexStart = -hit -1;
                        indexStart = indexStart > positionList.size() ? positionList.size() : indexStart;
                    }else {
                        indexStart = hit;
                    }

                    hit = positionList.binarySearch(end);
                    if (hit == -1){
                        indexEnd = 0;
                    }else if (hit < -1){
                        indexEnd = -hit - 1;
                        indexEnd = indexEnd > positionList.size() ? positionList.size() : indexEnd;
                    }else {
                        indexEnd = hit;
                    }
                    count = indexEnd - indexStart;
                    sb.append(count).append("\t");
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

    enum InfoSite{
        Info_IndianDwarf(0),
        Info_WE(1), Info_DE(2), Info_FTT(3);

        int index;
        InfoSite(int index) {
            this.index=index;
        }

        public int getIndex() {
            return index;
        }
    }
    public static int find(int[][] A) {
        int[][] solution = new int[A.length][A.length];

        solution[0][0] = A[0][0];
        // fill the first row
        for (int i = 1; i < A.length; i++) {
            solution[0][i] = A[0][i] + solution[0][i - 1];
        }

        // fill the first column
        for (int i = 1; i < A.length; i++) {
            solution[i][0] = A[i][0] + solution[i - 1][0];
        }

        // path will be either from top or left, choose which ever is minimum
        for (int i = 1; i < A.length; i++) {
            for (int j = 1; j < A.length; j++) {
                solution[i][j] = A[i][j]
                        + Math.min(solution[i - 1][j], solution[i][j - 1]);
            }
        }
        return solution[A.length - 1][A.length - 1];
    }

}