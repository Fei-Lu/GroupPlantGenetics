package daxing;

import com.google.common.collect.Ordering;
import daxing.common.chrrange.ChrRange;
import daxing.common.genotype.GenoGrid;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.GenotypeTable;
import daxing.v2.localAncestryInfer.IndividualSource;
import daxing.v2.localAncestryInfer.WindowSource;
import gnu.trove.list.TIntList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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

        String genotypeFile = "/Users/xudaxing/Desktop/ABBA/chr2A_vmap2.1_test.recode.vcf";
        String fd_dxyFile = "/Users/xudaxing/Desktop/ABBA/chr2A_vmap2.1_test_LR0002_fd_dxy.txt";
        String groupByPop2IndividualFile="";
        String outFile = "";
        int[] taxaIndex={0,1,2,3};
        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
        System.out.println();
//        IndividualSource individualSource = new IndividualSource(fd_dxyFile, "LR_0002");
//        WindowSource[] chrRangeList=individualSource.selectCandidateWindow(3);
//        System.out.println();

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