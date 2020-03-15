package xiaohan.rareallele;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.Range;
import pgl.infra.range.RangeValStr;
import pgl.infra.utils.PStringUtils;
import smile.stat.Stat;
import xujun.analysis.rnaseq.GeneFeature;

import java.io.BufferedReader;
import java.io.BufferedWriter;
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
        //this.vcffiltering();
        //this.vcfmerge();
        this.addinfo();
    }
    public void addinfo(){
        String infile = "/data2/xiaohan/nonfiltering/getsub/all.nonfilterDS.vcf";
        String output = "/data2/xiaohan/nonfiltering/getsub";
        String infor = "/data2/xiaohan/nonfiltering/getsub/all.nonfilter.vcf";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infor);
        BufferedWriter bw = IOUtils.getTextWriter(new File(output,"all.non-DS.vcf").getAbsolutePath());
        String[] temps = null;
        String temp = null;
        String[] tems = null;
        String[] temps1 = null;
        String temp1 = null;
//        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
//        String[] names = name.split(",");
        try{
            for(int i =0;i<42;i++){
                temp1 = br1.readLine();
            }
            while ((temp = br.readLine()) != null) {
                temp1 = br1.readLine();
                temps1 = temp1.split("\t");
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    bw.write("#");
                    for (int i = 0; i < 2; i++) {
                        tems = temps[i].split("]");
                        bw.write(tems[1] + "\t");
                    }
                    bw.write("ID"+ "\t");
                    for (int i = 2; i < 4; i++) {
                        tems = temps[i].split("]");
                        bw.write(tems[1] + "\t");
                    }
                    bw.write(temps1[5]+"\t"+temps1[6]+"\t"+temps1[7]+"\t"+temps1[8]+"\t");
                    for (int i = 4; i < 100; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    bw.newLine();
                    continue;
                } else {
                    temps = temp.split("\t");
                    temps1 = temp1.split("\t");
                    for (int i = 0; i < 2; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    bw.write("snp_" + temps[1] + "\t");
                    for (int i = 2; i < 4; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    bw.write(temps1[5]+"\t"+temps1[6]+"\t"+temps1[7]+"\t"+"DS"+"\t");
                    for (int i = 4; i < 100; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    bw.newLine();
                    continue;
                }
            }
            bw.flush();bw.close();
            br.close();
        }
        catch (Exception e ){
            e.printStackTrace();
        }
    }

    public void vcfmerge(){
        StringBuilder sb = new StringBuilder();
        sb.append("nohup vcf-concat ");
        for( int i = 0;i<45;i++){
            sb.append("/data3/wgs/vcf/GATK/vmap3/1.SNP/"+i+".snp.vcf.gz ");
        }
        sb.append("| bgzip -c > all.non-filtering.vcf.gz &");
        System.out.println(sb.toString());
    }
    public void vcffiltering(){
        StringBuilder sb = new StringBuilder();
        for(int i = 0 ; i < 45;i++){
//            sb.append("nohup vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/");
//            sb.append(i+".snp.vcf.gz --maf 0 --max-maf 0.05 --out ");
//            sb.append("/data2/xiaohan/SNP/");
//            sb.append(i+".snp.maf005 --recode "+" || ");
            sb.append("bgzip "+i+".snp.maf005.recode.vcf"+"\n");
//            sb.append("log"+i+".txt 2>&1 &"+"\n");
        }
        String words = sb.toString();
        System.out.println(words);
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
//        String path = "/data2/junxu/SiPASData/data_P101SC18112845-01-F004-B4-21/1.rawdata";
//        File file = new File(path);
//        mkshell(file);
    }
}

