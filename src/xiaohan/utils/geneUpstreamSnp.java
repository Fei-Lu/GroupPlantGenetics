package xiaohan.utils;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;

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

    public static Multimap<String, String> getSnpGeneMap(String infile) {
        HashMap<String, ArrayList<String>> siteGeneMap = new HashMap<>();
        Multimap<String, String> newmap = ArrayListMultimap.create();
        HashSet<String> siteSet = new HashSet<>();
        ArrayList<String> newlist = new ArrayList<>();
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        BufferedReader br1 = IOUtils.getTextGzipReader(infile);
//        BufferedReader br = null;
//        BufferedReader br1 = null;
//        if (infile.endsWith("gz")) {
//            br = IOUtils.getTextGzipReader(infile);
//            br1 = IOUtils.getTextGzipReader(infile);
//        }
//        br = IOUtils.getTextReader(infile);
//        br1 = IOUtils.getTextReader(infile);
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                String site = temps[0];
                String gene = temps[1].split(";")[0].split("=")[1];
                newmap.put(site, gene);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return newmap;
    }

    public static String[] getGenes(String infile) {
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextGzipReader(infile);
//        BufferedReader br = null;
//        if (infile.endsWith("gz")) {
//            br = IOUtils.getTextGzipReader(infile);
//        }
//        br = IOUtils.getTextReader(infile);
        HashSet<String> geneSet = new HashSet<>();
        String[] genes = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                String gene = temps[1].split(";")[0].split("=")[1];
                if(!geneSet.contains(gene)){
                    geneSet.add(gene);
                }
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
//        BufferedReader br = null;
//        BufferedReader br1 = null;
//        if (infile.endsWith("gz")) {
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        BufferedReader br1 = IOUtils.getTextGzipReader(infile);
//        }
//        br = IOUtils.getTextReader(infile);
//        br1 = IOUtils.getTextReader(infile);
        try {
            String site = null;
            String site1 = null;
            int countline = 0;
            while ((temp = br1.readLine()) != null) {
                temps = temp.split("\t");
                if (countline == 0) {
                    site = temps[0];
                    site1 = temps[0];
                    siteSet.add(site);
                    continue;
                }
                site = temps[0];
                if (!site.equals(site1)) {
                    siteSet.add(site);
                }
                site1 = site;
                countline++;
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
