package daxing.load.ancestralSite;

import daxing.common.IOTool;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public class IndividualChrLoad{

    String taxonName;
    int chr;
    GeneLoad[] geneLoads;

    public IndividualChrLoad(String taxonName, String[] geneNames, int chr){
        this.taxonName=taxonName;
        this.chr=chr;
        List<GeneLoad> geneLoads=new ArrayList<>(2000);
        for (String geneName : geneNames) {
            geneLoads.add(new GeneLoad(geneName));
        }
        Collections.sort(geneLoads);
        this.geneLoads=new GeneLoad[geneLoads.size()];
        for (int i = 0; i < this.geneLoads.length; i++) {
            this.geneLoads[i]=geneLoads.get(i);
        }
    }

    public int getGeneIndex(String geneName){
        return Arrays.binarySearch(geneLoads, new GeneLoad(geneName));
    }

    public void addGenotype(String geneName, byte[] siteGenotype){
        int geneIndex=getGeneIndex(geneName);
        geneLoads[geneIndex].addGenotype(siteGenotype);
    }

    public void addGenotype(String geneName, byte[] siteGenotype, double corrRatio){
        int geneIndex=getGeneIndex(geneName);
        geneLoads[geneIndex].addGenotype(siteGenotype, corrRatio);
    }

    public String[] getGeneNames() {
        String[] geneName=new String[this.geneLoads.length];
        for (int i = 0; i < this.geneLoads.length; i++) {
            geneName[i]=this.geneLoads[i].getGeneName();
        }
        return geneName;
    }

    public int getChr() {
        return chr;
    }

    public String getTaxonName() {
        return taxonName;
    }

    public void write(String outDir){
        int chr=this.getChr();
        String taxonName=this.getTaxonName();
        File outFile=new File(outDir, "chr"+ PStringUtils.getNDigitNumber(3, chr)+"."+taxonName+".txt.gz");
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("geneName\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn\tnumDerivedInNonsyn\t" +
                    "numHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious\tnumHeterInHGDeleterious");
            bw.newLine();
            for (GeneLoad geneLoad : this.geneLoads) {
                bw.write(geneLoad.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
