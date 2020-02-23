package xiaohan.rareallele;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.Range;
import pgl.infra.range.RangeValStr;
import pgl.infra.utils.PStringUtils;
import smile.stat.Stat;
import xujun.analysis.rnaseq.GeneFeature;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Test {

    public Test() throws IOException {
        //this.findTaxon();
        //this.findSNPnumber();
        //this.findFalseSample();
        //this.callposition();
        //this.findDifference();
    }

    public static void mkshell(File file){
        File[] fl = file.listFiles();
        for(File f:fl) {
            if(f.isDirectory()){
                mkshell(f);
            }
            if(f.isFile()){
                System.out.println(f);
            }
        }
    }

    public void callposition(){
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/wheat_v1.1_Lulab.gtf";
        String geneNameS = null;
        int gfIndex=0;
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            Set<String> geneSet = new HashSet<String>();
            Set<Integer> chrSet = new HashSet<>();
            String[] tem = null;
            String geneName = null;
            String[] geneNames = null;
            HashMap<String, Integer> geneChrMap = new HashMap();
            HashMap<String, Integer> geneMinMap = new HashMap();
            HashMap<String, Integer> geneMaxMap = new HashMap();
            HashMap<String, Byte> geneStrandMap = new HashMap();
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                if (!tem[2].startsWith("exon")) continue;
                String[] te = tem[8].split(";");
                geneName = te[1].split("\"")[1].toString();
                if (!geneSet.contains(geneName)) {
                    geneMinMap.put(geneName, Integer.MAX_VALUE);
                    geneMaxMap.put(geneName, Integer.MIN_VALUE);
                    geneSet.add(geneName);
                }
                //geneSet.add(geneName);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                if (geneMinMap.get(geneName) > min) {
                    geneMinMap.put(geneName, min);
                }
                if (geneMaxMap.get(geneName) < max) {
                    geneMaxMap.put(geneName, max);
                }
                int chr = Integer.parseInt(tem[0]);
                chrSet.add(chr);
                geneChrMap.put(geneName, chr);
                if (tem[6].startsWith("-")) geneStrandMap.put(geneName, (byte) 1);
                else geneStrandMap.put(geneName, (byte) 0);
            }
            geneNames = geneSet.toArray(new String[geneSet.size()]);
            Arrays.sort(geneNames);
            for (int i =0 ; i< geneNames.length;i++) {
                String geneS = geneNames[i];
                System.out.println(geneChrMap.get(geneS)+"\t"+geneMinMap.get(geneS)+"\t"+geneMaxMap.get(geneS)+"\t"+geneS);
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            System.out.println(geneNameS+"\t"+gfIndex);
        }
    }

    public void findDifference() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/phenotype.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedReader br1 = IOUtils.getTextReader(infileS);
        try {
            String temp = br.readLine();
            String temp1 = br1.readLine();
            temp = br.readLine();
            temp1 =br1.readLine();
            while((temp = br.readLine()) != null) {
                String[] temps = temp.split("\t");
                String[] temps1 =temp1.split("\t");
                String chr =temps[0];
                String chr1 =temps1[0];
                String largerNumber = temps[1];
                String smallNumber = temps1[1];
                int larger = Integer.parseInt(largerNumber);
                int small = Integer.parseInt(smallNumber);
                if (chr.equals(chr1) && larger < small){
                    System.out.println(temps1[0]+"\t"+temps1[1]);
                }
                temp1 = br1.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void findFalseSample(){
        String infileS = "/Users/yxh/Documents/eQTL/003experiment/浓度汇总/RNAconcentration.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String[] temps = null ;
        String temp = null ;
        HashSet<String> NameSet = new HashSet<String>() ;
        HashMap<String, Double> A280Map = new HashMap<>();
        HashMap<String, Double> A230Map = new HashMap<>();
        HashMap<String, Double> ConcentrationMap = new HashMap<>();
        try{
            while((temp = br.readLine()) != null){
                if(temp.startsWith("#")){
                    continue;
                }
                temps = temp.split("\t");
                String Name = temps[0];
                Double A280 = Double.valueOf(temps[1]);
                Double A230 = Double.valueOf(temps[2]);
                Double Concentration = Double.valueOf(temps[3]);
                NameSet.add(Name);
                A280Map.put(Name,A280);
                A230Map.put(Name,A230);
                ConcentrationMap.put(Name,Concentration);
                //A280不合格，A230合格
                /*if((A280 > 2.3 || A280 < 1.6) && (A230<2.5 && A230 > 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }
                //A280合格，A230不合格
                if((A280 < 2.3 && A280 > 1.6) && (A230 > 2.5 || A230 < 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                //A280，A230都不合格
                /*if((A280 > 2.3 || A280 < 1.6) && (A230 > 2.5 || A230 < 1.6)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                //浓度不合格
                /*if(Concentration < 30){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                /*if((A280 > 2.3 || A280 < 1.6) || (A230 > 2.5 || A230 < 1.6) || (Concentration < 30)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }*/
                if((A280 > 2.3 || A280 < 1.6) || (A230 > 2.5 || A230 < 1.6) || (Concentration < 30)){
                System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }
                
            }
            
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    public void findTaxon(){
    }

    public void findSNPnumber() throws IOException {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/output/genotype.vcf";
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        int count = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                } else {
                    count = count++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
            {
            }
        }
    }

    public static void main (String[] args) throws IOException {
        new Test();
        //String path = "/Users/yxh/Desktop/ANNO_ANCBJ170529_PM-ANCBJ170529-14_2020-02-19";
        //File file = new File(path);
        //mkshell(file);
    }
}

