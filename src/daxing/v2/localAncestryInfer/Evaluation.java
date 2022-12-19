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
     * True positive (TP)
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
                    if (temp.get(3).equals("NONE")) continue;
                    source = WindowSource.Source.valueOf(temp.get(3));
                    startIndex = this.getGenotypeTable().getSiteIndex("1", Integer.parseInt(temp.get(4)));
                    endIndex = this.getGenotypeTable().getSiteIndex("1", Integer.parseInt(temp.get(5)));
                    assert startIndex >=0 : "make sure start position is in genotype file";
                    assert endIndex >=0 : "make sure end position is in genotype file";
                    if (startIndex > endIndex){
                        System.out.println();
                    }
                    Arrays.fill(predictedDonors[i], startIndex, endIndex+1, source.getIndex());
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return predictedDonors;
    }

    public double[] calculateAccuracy(int[][] predictedValue, int[][] actualValue){
        double[] accuracy = new double[this.getTaxa().length];
        Arrays.fill(accuracy, -1);
        for (int i = 0; i < predictedValue.length; i++) {
            int count = 0;
            for (int j = 0; j < this.getVariantsNum(); j++) {
                if (predictedValue[i][j]!=actualValue[i][j]) continue;
                count++;
            }
            accuracy[i] = ((double)count)/this.getVariantsNum();
        }
        return accuracy;
    }

    public static void write_accuracy(String simulatedTractDir, String genotypeTableFile, String laidpResDir,
                                      String outFile_accuracy){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeTableFile);
        Evaluation evaluation = new Evaluation(simulatedTractDir, genotypeTable);
        int[][] predictedValue = evaluation.getPredictedValue(laidpResDir);
        int[][] actualValue = evaluation.getActualValue();
        double[] accuracy = evaluation.calculateAccuracy(predictedValue, actualValue);
        try (BufferedWriter bw = IOTool.getWriter(outFile_accuracy)) {
            StringBuilder sb = new StringBuilder();
            NumberFormat numberFormat = NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            String[] taxa = evaluation.getTaxa();
            for (int i = 0; i < accuracy.length; i++) {
                sb.setLength(0);
                sb.append(taxa[i]).append("\t").append(numberFormat.format(accuracy[i]));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
