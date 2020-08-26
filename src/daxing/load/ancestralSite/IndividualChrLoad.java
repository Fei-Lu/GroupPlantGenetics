package daxing.load.ancestralSite;

import daxing.common.IOTool;
import gnu.trove.list.array.TIntArrayList;
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
        for (int i = 0; i < geneNames.length; i++) {
            geneLoads.add(new GeneLoad(geneNames[i]));
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

    public String[] getGeneNames() {
        String[] geneName=new String[this.geneLoads.length];
        for (int i = 0; i < this.geneLoads.length; i++) {
            geneName[i]=this.geneLoads[i].getGeneName();
        }
        return geneName;
    }

    public TIntArrayList getNumHGDeleterious() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getHGDeleteriousNum());
        }
        return res;
    }

    public int getChr() {
        return chr;
    }

    public String getTaxonName() {
        return taxonName;
    }

    public TIntArrayList getNumDerivedInHGDeleterious() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getHGDeleteriousDerivedNum());
        }
        return res;
    }

    public TIntArrayList getNumDerivedInNonsyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getNonsynDerivedNum());
        }
        return res;
    }

    public TIntArrayList getNumDerivedInSyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getSynDerivedNum());
        }
        return res;
    }

    public TIntArrayList getNumNonsyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getNonsynNum());
        }
        return res;
    }

    public TIntArrayList getNumSyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getSynNum());
        }
        return res;
    }

    public TIntArrayList getSynHeterSitesNum(){
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getSynHeterSitesNum());
        }
        return res;
    }

    public TIntArrayList getNonsynHeterSitesNum(){
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getNonsynHeterSitesNum());
        }
        return res;
    }

    public TIntArrayList getHGDeleteriousHeterSitesNum(){
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            res.add(geneLoads[i].getHGDeleteriousHeterSitesNum());
        }
        return res;
    }

    public void write(String outDir){
        int chr=this.getChr();
        String taxonName=this.getTaxonName();
        File outFile=new File(outDir, "chr"+ PStringUtils.getNDigitNumber(3, chr)+"."+taxonName+".txt.gz");
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("geneName\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn\tnumDerivedInNonsyn\t" +
                    "numHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious\tnumHeterInHGDeleterious");
            bw.newLine();
            StringBuilder sb;
            String[] geneNames=this.getGeneNames();
            TIntArrayList numSyn=this.getNumSyn();
            TIntArrayList numDerivedInSyn=this.getNumDerivedInSyn();
            TIntArrayList numHeterInSyn=this.getSynHeterSitesNum();
            TIntArrayList numNonsyn=this.getNumNonsyn();
            TIntArrayList numDerivedInNonsyn=this.getNumDerivedInNonsyn();
            TIntArrayList numHeterInNonsyn=this.getNonsynHeterSitesNum();
            TIntArrayList numHGDeleterious=this.getNumHGDeleterious();
            TIntArrayList numDerivedInHGDeleterious=this.getNumDerivedInHGDeleterious();
            TIntArrayList numHeterInHGDeleterious=this.getHGDeleteriousHeterSitesNum();
            for (int i = 0; i < geneNames.length; i++) {
                sb=new StringBuilder();
                sb.append(geneNames[i]).append("\t").append(numSyn.get(i)).append("\t");
                sb.append(numDerivedInSyn.get(i)).append("\t");
                sb.append(numHeterInSyn.get(i)).append("\t");
                sb.append(numNonsyn.get(i)).append("\t");
                sb.append(numDerivedInNonsyn.get(i)).append("\t");
                sb.append(numHeterInNonsyn.get(i)).append("\t");
                sb.append(numHGDeleterious.get(i)).append("\t");
                sb.append(numDerivedInHGDeleterious.get(i)).append("\t");
                sb.append(numHeterInHGDeleterious.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
