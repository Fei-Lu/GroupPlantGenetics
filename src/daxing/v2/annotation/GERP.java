package daxing.v2.annotation;

import com.ibm.icu.text.NumberFormat;
import daxing.common.factors.FunctionalElement;
import daxing.common.utiles.IOTool;
import daxing.common.wheat.PGF;
import daxing.common.table.RowTableTool;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;

/**
 * GERP 16 species
 * Each line consists of the neutral rate N (from step 2) and RS score S (from step 4), separated by a tab character
 */
public class GERP {

    public static void split_toChrID(String gerpSortedBedDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(gerpSortedBedDir, "gz");
        Map<Integer, BufferedWriter> chrBWMap = new HashMap<>();
        BufferedWriter bw = null;
        String chr;
        int[] chrIDArray;
        File outFile;
        try {
            for (int i = 0; i < files.size(); i++){
                chr = files.get(i).getName().substring(0,2);
                chrIDArray= RefV1Utils.getChrIDs();
                for (int e : chrIDArray){
                    outFile = new File(outDir, files.get(i).getName().replaceAll(chr,
                            "chr"+PStringUtils.getNDigitNumber(3
                            , e)));
                    bw = IOTool.getWriter(outFile);
                    bw.write("ChrID\tPos\tRSScore");
                    bw.newLine();
                    chrBWMap.put(e, bw);
                }
            }
            BufferedReader br;
            String line, refChr;
            List<String> temp;
            int chrIDStart, chrIDEnd, posStart, posEnd, refStartPos, refEndPos ;
            double gerpScore;
            List<Integer> chrIDList;
            StringBuilder sb = new StringBuilder();
            for (File file: files){
                br = IOTool.getReader(file);
                chrIDList = new ArrayList<>();
                while ((line = br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    refChr = temp.get(0).substring(3,5);
                    refStartPos = Integer.parseInt(temp.get(1))+1;
                    refEndPos = Integer.parseInt(temp.get(2));
                    gerpScore = Double.parseDouble(temp.get(3));
                    if (gerpScore != 0){
                        chrIDEnd = RefV1Utils.getChrID(refChr, refEndPos);
                        posEnd = RefV1Utils.getPosOnChrID(refChr, refEndPos);
                        chrIDList.add(chrIDEnd);
                        bw = chrBWMap.get(chrIDEnd);
                        sb.setLength(0);
                        sb.append(chrIDEnd).append("\t").append(posEnd).append("\t");
                        sb.append(gerpScore);
                        bw.write(sb.toString());
                        bw.newLine();
                    }else {
                        for (int i = refStartPos; i < refEndPos+1; i++) {
                            posStart = RefV1Utils.getPosOnChrID(refChr,i);
                            chrIDStart = RefV1Utils.getChrID(refChr, i);
                            bw = chrBWMap.get(chrIDStart);
                            sb.setLength(0);
                            sb.append(chrIDStart).append("\t").append(posStart).append("\t");
                            sb.append(gerpScore);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                br.close();
            }
            for (Map.Entry<Integer, BufferedWriter> entry : chrBWMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void calculateConstraintRegionSize(String inputGERPDir, String outFile){
        List<File> files= IOTool.getVisibleDir(inputGERPDir);
        TDoubleArrayList gerpRes=new TDoubleArrayList();
        for (File file: files){
            gerpRes.add(calculateGERPGreaterThanZeroCount(file.getAbsolutePath()));
        }
        double[] gerpArray= gerpRes.toArray();
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setGroupingUsed(false);
        numberFormat.setMinimumFractionDigits(5);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < gerpArray.length; i++) {
                sb.setLength(0);
                sb.append("chr").append(PStringUtils.getNDigitNumber(3, i+1)).append("\t").append(gerpArray[i]);
                sb.append("\t").append(RefV1Utils.getChrIDLength(i+1)).append("\t");
                sb.append(numberFormat.format(gerpArray[i]/RefV1Utils.getChrIDLength(i+1)));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int calculateGERPGreaterThanZeroCount(String gerpFile){
        int count=0;
        try (BufferedReader br = IOTool.getReader(gerpFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if (Double.parseDouble(temp.get(2)) > 0){
                    count++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return count;
    }

    public static void addGERPGeneFuture(String gerpInputDir, String outDir, String pgfFile, String nonoverlapGeneFile){
        PGF pgf = new PGF(pgfFile);
        System.out.println(pgf.getGeneNumber());
        RowTableTool<String> table= new RowTableTool<>(nonoverlapGeneFile);
        Predicate<List<String>> uniqueGeneP= l->l.get(3).equals("1");
        table.removeIf(uniqueGeneP.negate());
        Set<String> uniqueGeneSet= new HashSet<>(table.getColumn(0));
        Predicate<PGF.Gene> duplicatedGeneP= gene -> !uniqueGeneSet.contains(gene.getGeneName());
        pgf.removeIf(duplicatedGeneP);
        pgf.sortGeneByGeneRange();
        System.out.println(pgf.getGeneNumber());
        List<File> files = IOTool.getFileListInDirEndsWith(gerpInputDir, ".gz");
        String[] outNames =
                files.stream().map(File::getName).map(s -> s.replaceAll(".sort.gz",".GeneFuture.txt.gz")).toArray(String[]::new);
        BufferedWriter bw;
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            int chrID, pos, geneIndex, longestTranscriptIndex;
            int utr5Index, utr3Index, cdsIndex;
            FunctionalElement functionalElement;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < files.size(); i++) {
                br = IOTool.getReader(files.get(i));
                br.readLine();
                bw = IOTool.getWriter(new File(outDir, outNames[i]));
                bw.write("ChrID\tPos\tGERP\tGeneFuture");
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    chrID = Integer.parseInt(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    geneIndex = pgf.getGeneIndex(chrID, pos);
                    if (geneIndex < 0){
                        functionalElement = FunctionalElement.INTERGENIC;
                        sb.setLength(0);
                        sb.append(String.join("\t", temp)).append("\t").append(functionalElement.name());
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    longestTranscriptIndex = pgf.getLongestTranscriptIndex(geneIndex);
                    utr5Index = pgf.getGene(geneIndex).getTs().get(longestTranscriptIndex).get5UTRIndex(chrID, pos);
                    utr3Index = pgf.getGene(geneIndex).getTs().get(longestTranscriptIndex).get3UTRIndex(chrID, pos);
                    cdsIndex = pgf.getGene(geneIndex).getTs().get(longestTranscriptIndex).getCDSIndex(chrID, pos);
                   if (utr5Index >=0){
                        functionalElement = FunctionalElement.UTR5;
                    }else if (utr3Index >=0){
                        functionalElement = FunctionalElement.UTR3;
                    }else if (cdsIndex >=0){
                        functionalElement= FunctionalElement.CDS;
                    }else {
                        functionalElement= FunctionalElement.INTRON;
                    }
                    sb.setLength(0);
                    sb.append(String.join("\t", temp)).append("\t").append(functionalElement.name());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
