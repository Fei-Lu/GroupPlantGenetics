/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import analysis.maizeRNASeq.*;
import format.table.RowTable;
import format.tree.Newick;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOUtils;

/**
 *
 * @author feilu
 */
public class TaxaDiversitycp {
    
    public TaxaDiversitycp () {
        //this.selectedTaxaFromDistanceMatrix();
        //this.selectedTaxaFromDistanceMatrix2();
        //this.selectedTaxaFromNewickTree();
        this.mdsOnSelectedTaxa();
        //this.colorSelectedTaxaOnTree();
    }
    
    public void colorSelectedTaxaOnTree () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/source/Tree_chr10000sites_432taxa.xml";
        String outfileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/Tree_chr10000sites_432taxa_colorSel.xml";
        String taxaFileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/selectedTaxa.txt";
        String addS = "<property ref=\"style:font_color\" datatype=\"xsd:token\" applies_to=\"node\">#ff0000</property>";
        RowTable<String> t = new RowTable(taxaFileS);
        List<String> taxaList = t.getColumn(0);
        Collections.sort(taxaList);
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.contains("<name>")) {
                    String currentName = temp.split(">")[1];
                    currentName = currentName.split("<")[0];
                    bw.write(temp);
                    bw.newLine();
                    bw.write(br.readLine());
                    bw.newLine();
                    if (Collections.binarySearch(taxaList, currentName) < 0) {
                        
                    }
                    else {
                        bw.write(addS);
                        bw.newLine();
                    } 
                }
                else {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();

        }
        catch (Exception e) {
            e.printStackTrace();
            
        }
    }
    
    public void mdsOnSelectedTaxa () {
        String infileS = "/Users/Aoyue/Documents/tree/selectDiversity/Fei_result/selectedTaxa.txt";
        String mdsInfileS = "/Users/Aoyue/Documents/mds/851oriMds.txt";
        String outfileS = "/Users/Aoyue/Documents/mds/300oriMds.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> taxaList = t.getColumn(0);
        Collections.sort(taxaList);
        try {
            BufferedReader br = IOUtils.getTextReader(mdsInfileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa\tPC1\tPC2\tPC3");
            bw.newLine();
            String temp = br.readLine();
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split("\t");
                if (Collections.binarySearch(taxaList, tem[0]) < 0) continue;
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void selectedTaxaFromNewickTree () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/source/Tree_chr10000sites_432taxa.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/selectedTaxa.txt";
        String[] vip = {"Jing2416", "mPH6WC", "Huangzaosi", "B73", "Mo17Ht"};     
        int n = 300;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = br.readLine();
            br.close();
            //temp = "((raccoon:19.19959,bear:6.80041):0.84600,((sea_lion:11.99700, seal:12.00300):7.52973,((monkey:100.85930,cat:47.14069):20.59201, weasel:18.87953):2.09460):3.87382,dog:25.46154)";
            Newick nwk = new Newick (temp);
            List<String> taxaList = nwk.selectTaxaWithMaxDiversity(n);
            Set<String> taxaSet = new HashSet(taxaList);
            for (int i = 0; i < vip.length; i++) {
                taxaSet.add(vip[i]);
            }
            taxaList = new ArrayList(taxaSet);
            Collections.sort(taxaList);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < taxaList.size(); i++) {
                bw.write(taxaList.get(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    public void selectedTaxaFromDistanceMatrix2 () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/source/Matrix_chr10000sites_432taxa.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/selectedTaxa.txt";
        int selectedNum = 100;
        double[] dis = null;
        String[] taxaNames = null;
        int taxaNum = 0;
        TaxonDistance[] td = null;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                break;
            }
            taxaNum = Integer.parseInt(temp);
            dis = new double[taxaNum];
            taxaNames = new String[taxaNum];
            td = new TaxonDistance[taxaNum];
            for (int i = 0; i < taxaNum; i++) {
                String[] tem = br.readLine().split("\t");
                taxaNames[i] = tem[0];
                for (int j = 0; j < taxaNum; j++) {
                    dis[i]+= Double.parseDouble(tem[j+1]);
                }
                td[i] = new TaxonDistance(tem[0], dis[i]);
            }
            br.close();
            Arrays.sort(td, Collections.reverseOrder());
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < selectedNum; i++) {
                bw.write(td[i].taxon);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    class TaxonDistance implements Comparable<TaxonDistance>{
        String taxon = null;
        double dis = 0;
        public TaxonDistance (String taxon, double dis) {
            this.taxon = taxon;
            this.dis = dis;
        }
        @Override
        public int compareTo(TaxonDistance o) {
            if (this.dis < o.dis) return -1;
            else if (this.dis > o.dis) return 1;
            return 0;
        }
    }
    
    //Doesn't work
    public void selectedTaxaFromDistanceMatrix () {
        String infileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/source/Matrix_chr10000sites_432taxa.txt";
        String outfileS = "/Users/feilu/Documents/analysisL/production/maizeRNA/sampleSelection/selectedTaxa.txt";
        int selectedNum = 100;
        double[][] dis = null;
        String[] taxaNames = null;
        int taxaNum = 0;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                break;
            }
            taxaNum = Integer.parseInt(temp);
            dis = new double[taxaNum][taxaNum];
            taxaNames = new String[taxaNum];
            for (int i = 0; i < taxaNum; i++) {
                String[] tem = br.readLine().split("\t");
                taxaNames[i] = tem[0];
                for (int j = 0; j < taxaNum; j++) {
                    dis[i][j] = Double.parseDouble(tem[j+1]);
                }
            }
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        List<TwoTaxaDistance> dList = new ArrayList<>();
        for (int i = 0; i < taxaNum-1; i++) {
            for (int j = i+1; j < taxaNum; j++) {
                TwoTaxaDistance d = new TwoTaxaDistance(taxaNames[i], taxaNames[j], dis[i][j]);
                dList.add(d);
            }
        }
        Collections.sort(dList);
        Set<String> removedSet = new HashSet<String>();
        while (this.getTaxaNum(dList)>selectedNum) {
            String removingTaxon = null;
            for (int i = 0; i < dList.size(); i++) {
                TwoTaxaDistance current = dList.get(i);
                if (current.taxon1 != null && current.taxon2 != null) {
                    removingTaxon = current.taxon1;
                    current.taxon1 = null;
                    removedSet.add(current.taxon1);
                    break;
                }
            }
            
            for (int i = 0; i < dList.size(); i++) {
                TwoTaxaDistance current = dList.get(i);
                if (removingTaxon.equals(current.taxon1)) current.taxon1 = null;
                if (removingTaxon.equals(current.taxon2)) current.taxon2 = null;
            }
        }
        Set<String> taxaSet = new HashSet();
        for (int i = 0; i < dList.size(); i++) {
            TwoTaxaDistance current = dList.get(i);
            if (current.taxon1 != null) taxaSet.add(current.taxon1);
            if (current.taxon2 != null) taxaSet.add(current.taxon2);
        }
        String[] remainingTaxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(remainingTaxa);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Taxa");
            bw.newLine();
            for (int i = 0; i < remainingTaxa.length; i++) {
                bw.write(remainingTaxa[i]);
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
    private int getTaxaNum (List<TwoTaxaDistance> dList) {
        Set<String> taxaSet = new HashSet();
        for (int i = 0; i < dList.size(); i++) {
            TwoTaxaDistance current = dList.get(i);
            if (current.taxon1 != null) taxaSet.add(current.taxon1);
            if (current.taxon2 != null) taxaSet.add(current.taxon2);
        }
        return taxaSet.size();
    }
    
    class TwoTaxaDistance implements Comparable<TwoTaxaDistance> {
        String taxon1 = null;
        String taxon2 = null;
        double distance;
        boolean ifProsessed = false;
        boolean taxon1Removed = false;
        boolean taxon2Removed = false;
        
        public TwoTaxaDistance (String taxon1, String taxon2, double distance) {
            this.taxon1 = taxon1;
            this.taxon2 = taxon2;
            this.distance = distance;
        }

        @Override
        public int compareTo(TwoTaxaDistance o) {
            if (this.distance < o.distance) return -1;
            else if (this.distance > o.distance) return 1;
            return 0;
        }
        
        public boolean isThereNullTaxon () {
            if (taxon1 == null || taxon2 == null) return true;
            return false;
        }
        
        public void setTaxon1Removed () {
            taxon1Removed = true;
        }
        
        public void setTaxon2Removed () {
            taxon1Removed = true;
        }
        
    }
}
