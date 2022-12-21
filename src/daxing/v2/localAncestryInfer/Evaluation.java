package daxing.v2.localAncestryInfer;

import daxing.common.chrrange.ChrRange;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class Evaluation {

    String[] taxa;
    List<ChrRange>[] introgressionRange;
    IntList[] donors;
    GenotypeTable genotypeTable;

    public Evaluation(String simulatedTractDir, GenotypeTable genotypeTable){
        List<File> simulatedFiles = IOTool.getFileListInDirEndsWith(simulatedTractDir, ".txt");
        String[] taxa = simulatedFiles.stream().map(File::getName).map(s -> s.substring(10, 12)).toArray(String[]::new);
        List<ChrRange>[] introgressionRange = new List[taxa.length];
        for (int i = 0; i < introgressionRange.length; i++) {
            introgressionRange[i] = new ArrayList<>();
        }
        IntList[] donors  = new IntList[taxa.length];
        for (int i = 0; i < donors.length; i++) {
            donors[i] = new IntArrayList();
        }
        try {
            BufferedReader br;
            String line;
            List<String> temp;
            ChrRange chrRange;
            for (int i = 0; i < simulatedFiles.size(); i++) {
                br = IOTool.getReader(simulatedFiles.get(i));
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    chrRange = new ChrRange("1", Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)));
                    introgressionRange[i].add(chrRange);
                    donors[i].add(0);
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        this.introgressionRange=introgressionRange;
        this.donors=donors;
        this.taxa=taxa;
        this.genotypeTable=genotypeTable;
    }

    public GenotypeTable getGenotypeTable() {
        return genotypeTable;
    }

    public List<ChrRange>[] getIntrogressionRange() {
        return introgressionRange;
    }

    public IntList[] getDonors() {
        return donors;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public int getVariantsNum(){
        return this.getGenotypeTable().getSiteNumber();
    }

    /**
     *
     * @return actual value, dim1 is taxa, dim2 is variants source
     *         variants source equal WindowSource.Source.getIndex()
     */
    public int[][] getActualValue(){
        int[][] actualDonors = new int[this.getTaxa().length][];
        for (int i = 0; i < actualDonors.length; i++) {
            actualDonors[i] = new int[this.getVariantsNum()];
            Arrays.fill(actualDonors[i], 4);
        }
        int start, end, startHit, endHit, startIndex, endIndex;
        WindowSource.Source source;
        List<ChrRange>[] chrRanges = this.getIntrogressionRange();
        IntList[] donors = this.getDonors();
        for (int i = 0; i < chrRanges.length; i++) {
            for (int j = 0; j < chrRanges[i].size(); j++) {
                start = chrRanges[i].get(j).getStart();
                end = chrRanges[i].get(j).getEnd();
                source = WindowSource.Source.getInstanceFromSubNum(donors[i].getInt(j)).get();
                startHit = this.getGenotypeTable().getSiteIndex("1", start);
                endHit = this.getGenotypeTable().getSiteIndex("1", end);
                startIndex = startHit < 0 ? -startHit-1 : startHit;
                endIndex = endHit < 0 ? -endHit-1 : endHit;
                Arrays.fill(actualDonors[i], startIndex, endIndex, source.getIndex());
            }
        }
        return actualDonors;
    }

    /**
     *
     * @param laidpResDir
     * @return predicted value, dim1 is taxa, dim2 is variants source
     * variants source equal WindowSource.Source.getIndex()
     */
    public int[][] getPredictedValue(String laidpResDir){
        List<File> files = IOTool.getFileListInDirEndsWith(laidpResDir, ".txt");
        String[] taxa = files.stream().map(File::getName).map(s -> s.substring(9, 11)).toArray(String[]::new);
        int[][] predictedDonors = new int[taxa.length][];
        for (int i = 0; i < predictedDonors.length; i++) {
            predictedDonors[i] = new int[this.getVariantsNum()];
            Arrays.fill(predictedDonors[i], 4);
        }
        try {
            BufferedReader br;
            String line;
            List<String> temp;
            WindowSource.Source source;
            int startIndex, endIndex;
            for (int i = 0; i < files.size(); i++) {
                br = IOTool.getReader(files.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    if (temp.get(3).contains(",")) continue;
                    if (temp.get(3).equals("4")) continue;
                    source = WindowSource.Source.getInstanceFromSubNum(Integer.parseInt(temp.get(3))).get();
                    startIndex = this.getGenotypeTable().getSiteIndex("1", Integer.parseInt(temp.get(4)));
                    endIndex = this.getGenotypeTable().getSiteIndex("1", Integer.parseInt(temp.get(5)));
                    assert startIndex >=0 : "make sure start position is in genotype file";
                    assert endIndex >=0 : "make sure end position is in genotype file";
                    Arrays.fill(predictedDonors[i], startIndex, endIndex+1, source.getIndex());
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return predictedDonors;
    }

    public double[][] calculateAccuracyRecallPrecision(int[][] predictedValue, int[][] actualValue){
        double[][] accuracy_recall_precision = new double[5][];
        for (int i = 0; i < accuracy_recall_precision.length; i++) {
            accuracy_recall_precision[i] = new double[this.getTaxa().length];
            Arrays.fill(accuracy_recall_precision[i], -1);
        }
        for (int i = 0; i < predictedValue.length; i++) {
            int count_truePositive=0;
            int count_falseNegative=0;
            int count_falsePositive=0;
            int count_trueNegative=0;
            for (int j = 0; j < this.getVariantsNum(); j++) {
                if ((predictedValue[i][j]!=4) && (actualValue[i][j] == predictedValue[i][j])){
                    count_truePositive++;
                }else if ((predictedValue[i][j]==4) && (actualValue[i][j]!=4)){
                    count_falseNegative++;
                }else if ((predictedValue[i][j]!=4) && (actualValue[i][j]==4)){
                    count_falsePositive++;
                }else if ((predictedValue[i][j]==4) && (actualValue[i][j]==4)){
                    count_trueNegative++;
                }
            }
            accuracy_recall_precision[0][i] = ((double)(count_truePositive+count_trueNegative))/this.getVariantsNum();
            accuracy_recall_precision[1][i] = ((double)count_truePositive)/(count_truePositive+count_falseNegative);
            accuracy_recall_precision[2][i] = ((double)count_truePositive)/(count_truePositive+count_falsePositive);
            accuracy_recall_precision[3][i] = ((double)count_falsePositive)/(count_falsePositive+count_trueNegative);
            accuracy_recall_precision[4][i] = ((double)count_trueNegative)/(count_falsePositive+count_trueNegative);

        }
        return accuracy_recall_precision;
    }

    /**
     *
     * @param predictedValue
     * @param actualValue
     * @return contingency table, dim1 is [True positive, False negative, False positive, True negative], dim2 is taxa
     */
    public int[][] calculateContingencyTable(int[][] predictedValue, int[][] actualValue){
        int[][] contingency = new int[4][];
        for (int i = 0; i < contingency.length; i++) {
            contingency[i] = new int[this.getTaxa().length];
            Arrays.fill(contingency[i], -1);
        }
        for (int i = 0; i < predictedValue.length; i++) {
            int count_truePositive=0;
            int count_falseNegative=0;
            int count_falsePositive=0;
            int count_trueNegative=0;
            for (int j = 0; j < this.getVariantsNum(); j++) {
                if ((predictedValue[i][j]!=4) && (actualValue[i][j] == predictedValue[i][j])){
                    count_truePositive++;
                }else if ((predictedValue[i][j]==4) && (actualValue[i][j]!=4)){
                    count_falseNegative++;
                }else if ((predictedValue[i][j]!=4) && (actualValue[i][j]==4)){
                    count_falsePositive++;
                }else if ((predictedValue[i][j]==4) && (actualValue[i][j]==4)){
                    count_trueNegative++;
                }
            }
            contingency[0][i] = count_truePositive;
            contingency[1][i] = count_falseNegative;
            contingency[2][i] = count_falsePositive;
            contingency[3][i] = count_trueNegative;

        }
        return contingency;
    }

    public static void write_accuracy_recall_precision(String simulatedTractDir, String genotypeTableFile, String laidpResDir,
                                                       String outFile_accuracy){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeTableFile);
        Evaluation evaluation = new Evaluation(simulatedTractDir, genotypeTable);
        int[][] predictedValue = evaluation.getPredictedValue(laidpResDir);
        int[][] actualValue = evaluation.getActualValue();
        double[][] accuracyRecallPrecision = evaluation.calculateAccuracyRecallPrecision(predictedValue, actualValue);
        try (BufferedWriter bw = IOTool.getWriter(outFile_accuracy)) {
            StringBuilder sb = new StringBuilder();
            sb.append("Taxa\tAccuracy\tRecall\tPrecision\tFalsePositiveRate\tSpecificity");
            bw.write(sb.toString());
            bw.newLine();
            NumberFormat numberFormat = NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            String[] taxa = evaluation.getTaxa();
            for (int i = 0; i < accuracyRecallPrecision[0].length; i++) {
                sb.setLength(0);
                sb.append(taxa[i]).append("\t").append(numberFormat.format(accuracyRecallPrecision[0][i])).append("\t");
                sb.append(numberFormat.format(accuracyRecallPrecision[1][i])).append("\t");
                sb.append(numberFormat.format(accuracyRecallPrecision[2][i])).append("\t");
                sb.append(numberFormat.format(accuracyRecallPrecision[3][i])).append("\t");
                sb.append(numberFormat.format(accuracyRecallPrecision[4][i]));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void write_ContingencyTable(String simulatedTractDir, String genotypeTableFile, String laidpResDir,
                                              String outFile_ContingencyTable){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeTableFile);
        Evaluation evaluation = new Evaluation(simulatedTractDir, genotypeTable);
        int[][] predictedValue = evaluation.getPredictedValue(laidpResDir);
        int[][] actualValue = evaluation.getActualValue();
        int[][] contingencyTable = evaluation.calculateContingencyTable(predictedValue, actualValue);
        try (BufferedWriter bw = IOTool.getWriter(outFile_ContingencyTable)) {
            StringBuilder sb = new StringBuilder();
            sb.append("Taxa\tTruePositive\tFalseNegative\tFalsePositive\tTrueNegative");
            bw.write(sb.toString());
            bw.newLine();
            NumberFormat numberFormat = NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            String[] taxa = evaluation.getTaxa();
            for (int i = 0; i < contingencyTable[0].length; i++) {
                sb.setLength(0);
                sb.append(taxa[i]).append("\t").append(numberFormat.format(contingencyTable[0][i])).append("\t");
                sb.append(numberFormat.format(contingencyTable[1][i])).append("\t");
                sb.append(numberFormat.format(contingencyTable[2][i])).append("\t");
                sb.append(numberFormat.format(contingencyTable[3][i]));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
