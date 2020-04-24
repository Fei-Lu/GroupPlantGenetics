package daxing.load.neutralSiteLoad;

import daxing.common.IOTool;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class IndividualChrLoad {

    String taxonName;
    String[] geneNames;
    int[] derivedCounts;

    public IndividualChrLoad(String taxonName,String[] geneNames){
        this.taxonName=taxonName;
        this.geneNames=geneNames;
        this.derivedCounts=new int[geneNames.length];
    }

    public void addGeneDerivedCount(int geneIndex, int derivedCount){
        this.derivedCounts[geneIndex]=this.derivedCounts[geneIndex]+derivedCount;
    }

    public void write(File parentDir, String outFile){
        try (BufferedWriter bw = IOTool.getTextWriter(new File(parentDir, outFile))) {
            bw.write("GeneName\tDerivedCountInGeneLocal");
            bw.newLine();
            StringBuilder sb;
            for (int i = 0; i < this.geneNames.length; i++) {
                sb=new StringBuilder();
                sb.append(geneNames[i]).append("\t").append(derivedCounts[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
