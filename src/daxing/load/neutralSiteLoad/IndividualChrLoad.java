package daxing.load.neutralSiteLoad;

import daxing.common.IOTool;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class IndividualChrLoad {

    String taxonName;
    String[] geneNames;
    int[][] derivedCounts;

    public IndividualChrLoad(String taxonName,String[] geneNames){
        this.taxonName=taxonName;
        this.geneNames=geneNames;
        this.derivedCounts=new int[geneNames.length][];
        for (int i = 0; i < derivedCounts.length; i++) {
            derivedCounts[i]=new int[2];
        }
    }

    public void addGeneDerivedCount(int geneIndex, int[] derivedCount){
        this.derivedCounts[geneIndex][0]=this.derivedCounts[geneIndex][0]+derivedCount[0];
        this.derivedCounts[geneIndex][1]=this.derivedCounts[geneIndex][1]+derivedCount[1];
    }

    public void write(File parentDir, String outFile){
        try (BufferedWriter bw = IOTool.getTextGzipWriter(new File(parentDir, outFile))) {
            bw.write("GeneName\tnumGeneLocal\tnumDerivedInGeneLocal");
            bw.newLine();
            StringBuilder sb;
            for (int i = 0; i < this.geneNames.length; i++) {
                sb=new StringBuilder();
                sb.append(geneNames[i]).append("\t").append(derivedCounts[i][0]).append("\t").append(derivedCounts[i][1]);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
