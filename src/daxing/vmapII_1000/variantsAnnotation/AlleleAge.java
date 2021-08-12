package daxing.vmapII_1000.variantsAnnotation;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class AlleleAge {

    public static void interSection(String geneSiteAnnoDBDir,String ancestraSiteDir, String outDir){
        List<File> geneSiteFiles = IOTool.getFileListInDirEndsWith(geneSiteAnnoDBDir, "gz");
        List<File> ancestralFiles = IOTool.getFileListInDirEndsWith(ancestraSiteDir, "txt");
        String[] outNames =
                geneSiteFiles.stream().map(File::getName).map(s -> s.replaceAll("txt.gz","ancestralSite.txt")).toArray(String[]::new);
        IntStream.range(0, geneSiteFiles.size()).forEach(e->{
            Set<String> geneSites = RowTableTool.getColumnSet(geneSiteFiles.get(e).getAbsolutePath(), 2);
            Set<String> ancestralSites = RowTableTool.getColumnSet(ancestralFiles.get(e).getAbsolutePath(), 0);
            List<String> ancestralSiteList= new ArrayList<>(ancestralSites);
            Collections.sort(ancestralSiteList, Comparator.comparingInt(value -> Integer.parseInt(value)));
            try (BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                for (String str : ancestralSiteList) {
                    if (!geneSites.contains(str)) continue;
                    bw.write(str);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void buildSH_alleleAge(String inputFile_bin, String positionsDir_42,
                                         String maxIndexPerChrFile,
                                         String alleleAgeSH){
        List<File> binFiles = IOTool.getFileListInDirEndsWith(inputFile_bin, "bin");
        List<File> positionFiles = IOTool.getFileListInDirEndsWith(positionsDir_42, "txt");
        Map<String, String> chrMaxIndexMap = RowTableTool.getMap(maxIndexPerChrFile, 0, 1);
        int chrID, maxIndex;
        StringBuilder sb = new StringBuilder();
        String outFileName, outFileNamePrefix;
        try (BufferedWriter bw = IOTool.getWriter(alleleAgeSH)) {
            for (int i = 0; i < binFiles.size(); i++) {
                chrID = Integer.parseInt(binFiles.get(i).getName().substring(3,6));
                maxIndex = Integer.parseInt(chrMaxIndexMap.get(String.valueOf(chrID)));
                for (int j = 0; j < maxIndex + 1; j++) {
                    outFileName = positionFiles.get(i).getName().replaceAll(".txt", "_"+ PStringUtils.getNDigitNumber(3, j));
                    outFileNamePrefix = outFileName+".run1";
                    sb.setLength(0);
                    sb.append("geva_v1beta -t 5 -i 002_convertFormat/").append(binFiles.get(i).getName()).append(" -o" +
                            " ");
                    sb.append("004_alleleAge/").append(outFileNamePrefix);
                    sb.append(" --positions 003_positionsByChrID/004_interSection_position_geneSiteAnnoDB_200Sites/");
                    sb.append(outFileName+".txt ");
                    sb.append("--maxConcordant 500 --maxDiscordant 500 --Ne 10000 --mut 6.5e-9 --hmm ");
                    sb.append("geva/hmm/hmm_initial_probs.txt ");
                    sb.append("geva/hmm/hmm_emission_probs.txt");
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
