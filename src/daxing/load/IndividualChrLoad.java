package daxing.load;

import daxing.common.IOTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;

public class IndividualChrLoad {

    String taxonName;
    String[] geneNames;
    int chr;
    GeneLoad[] geneLoads;

    public IndividualChrLoad(String taxonName, String[] geneNames, int chr){
        this.taxonName=taxonName;
        this.geneNames=geneNames;
        this.chr=chr;
        this.geneLoads=new GeneLoad[geneNames.length];
    }

    public void addGeneLoad(int geneIndex, GeneLoad geneLoad){
        geneLoads[geneIndex]=geneLoad;
    }

    public String[] getGeneNames() {
        return geneNames;
    }

    public TIntArrayList getNumHGDeleterious() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getHGDeleteriousNum());
            }
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
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getHGDeleteriousDerivedNum());
            }
        }
        return res;
    }

    public TIntArrayList getNumDerivedInNonsyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getNonsynDerivedNum());
            }
        }
        return res;
    }

    public TIntArrayList getNumDerivedInSyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getSynDerivedNum());
            }
        }
        return res;
    }

    public TIntArrayList getNumNonsyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getNonsynNum());
            }
        }
        return res;
    }

    public TIntArrayList getNumSyn() {
        TIntArrayList res=new TIntArrayList();
        for (int i = 0; i < geneLoads.length; i++) {
            if (geneLoads[i]==null){
                res.add(-1);
            }else {
                res.add(geneLoads[i].getSynNum());
            }
        }
        return res;
    }

    public void write(String outDir){
        int chr=this.getChr();
        String taxonName=this.getTaxonName();
        File outFile=new File(outDir, "chr"+ PStringUtils.getNDigitNumber(3, chr)+"."+taxonName+".txt.gz");
        try (BufferedWriter bw = IOTool.getTextGzipWriter(outFile)) {
            bw.write("transcriptName\tnumSyn\tnumDerivedInSyn\tnumNonsyn\tnumDerivedInNonsyn" +
                    "\tnumHGDeleterious\tnumDerivedInHGDeleterious");
            bw.newLine();
            StringBuilder sb;
            String[] geneNames=this.getGeneNames();
            TIntArrayList numSyn=this.getNumSyn();
            TIntArrayList numDerivedInSyn=this.getNumDerivedInSyn();
            TIntArrayList numNonsyn=this.getNumNonsyn();
            TIntArrayList numDerivedInNonsyn=this.getNumDerivedInNonsyn();
            TIntArrayList numHGDeleterious=this.getNumHGDeleterious();
            TIntArrayList numDerivedInHGDeleterious=this.getNumDerivedInHGDeleterious();
            for (int i = 0; i < geneNames.length; i++) {
                sb=new StringBuilder();
                sb.append(geneNames[i]).append("\t").append(numSyn.get(i)).append("\t");
                sb.append(numDerivedInSyn.get(i)).append("\t");
                sb.append(numNonsyn.get(i)).append("\t");
                sb.append(numDerivedInNonsyn.get(i)).append("\t");
                sb.append(numHGDeleterious.get(i)).append("\t");
                sb.append(numDerivedInHGDeleterious.get(i));
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
