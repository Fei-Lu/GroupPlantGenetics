package xiaohan.utils;

import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class geneUpstreamSnp {

    public geneUpstreamSnp(Object arg) {
        this.getSnpGeneMap((String) arg);
        this.getGenes((String) arg);
        this.getSites((String) arg);
    }

    public static HashMap<String, ArrayList<String>> getSnpGeneMap(String infile) {
        HashMap<String, ArrayList<String>> siteGeneMap = new HashMap<>();
        HashSet<String> siteSet = new HashSet<>();
        ArrayList<String> newlist = new ArrayList<>();
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infile);
        try {
            String site = null;
            while ((temp = br1.readLine()) != null) {
                temps = temp.split("\t");
                site = temps[0];
                if(!siteSet.contains(site)) {
                    siteSet.add(site);
                }
            }
            String[] sitestr = siteSet.toArray(new String[siteSet.size()]);
            Arrays.sort(sitestr);
            for (int i = 0; i < sitestr.length; i++) {
                siteGeneMap.put(sitestr[i], newlist);
            }
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                site = temps[0];
                ArrayList<String> list1 = siteGeneMap.get(site);
                String gene = temps[1].split(";")[0].split("=")[1];
                if(!list1.contains(gene)) {
                    list1.add(gene);
                }
                siteGeneMap.put(site,list1);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return siteGeneMap;
    }

    public static String[] getGenes(String infile) {
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextReader(infile);
        HashSet<String> geneSet = new HashSet<>();
        String[] genes = null;
        try {
            String gene = null;
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                gene = temps[1].split(";")[0].split("=")[1];
                geneSet.add(gene);
            }
            genes = geneSet.toArray(new String[geneSet.size()]);
            Arrays.sort(genes);
        } catch (Exception e) {
            e.printStackTrace();
        }
        return genes;
    }

    public static HashSet<String> getSites(String infile) {
        HashMap<String, ArrayList<String>> siteGeneMap = new HashMap<>();
        HashSet<String> siteSet = new HashSet<>();
        ArrayList<String> newlist = new ArrayList<>();
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infile);
        try {
            String site = null;
            String site1 = null;
            int countline = 0;
            while ((temp = br1.readLine()) != null) {
                temps = temp.split("\t");
                if(countline == 0) {
                    site = temps[0];
                    site1 = temps[0];
                    siteSet.add(site);
                    continue;
                }
                site = temps[0];
                if(!site.equals(site1)){
                    siteSet.add(site);
                }
                site1 = site;
                countline ++;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return siteSet;
    }

    public static void main(String[] args) {
        new geneUpstreamSnp(args[0]);
    }
}
