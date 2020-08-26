package xiaohan.wheatRNAseq;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.math.BigDecimal;
import java.util.*;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * @author yxh
 */
public class test {

    public test() throws IOException {
//        this.mkPosGeneMap();
//        this.writecode();
//        this.decimals();
//        this.clonefiles();
//        this.posAllele();
//        this.NewFile();
//        this.readFile();
//        this.countTable();
//        this.forfun();
//        this.test();
        this.test1();
//        this.refGene();
    }

    /**
     * @throws IOException
     */

    public void refGene() {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "refegeneWheat.txt").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String mRNA = null;
        int exonNumber = 0;
        int CDS = 1;
        StringBuilder sb1 = new StringBuilder();
        StringBuilder sb2 = new StringBuilder();
        String start = null;
        String end = null;
        int count = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##") && !temp.equals("###")) continue;
                if (!temp.startsWith("###")) {
                    temps = temp.split("\t");
                    if (temps[2].equals("gene")) {
                        String name = temps[8].split(";")[0].replace("ID=", "");
                        bw.write(name + "\t" + name + "\t" + temps[0] + "\t" + temps[7] + "\t");
                        continue;
                    }
                    if (temps[2].equals("mRNA") && mRNA == null) {
                        bw.write(temps[3] + "\t" + temps[4] + "\t");
                        mRNA = "mRNA";
                        continue;
                    }
                    if (temps[2].equals("CDS") & CDS == 1) {
                        start = temps[3];
                        end = temps[4];
                        CDS++;
                        continue;
                    }
                    if (temps[2].equals("CDS") & CDS > 1) {
                        start = start;
                        end = temps[4];
                        CDS++;
                        continue;
                    }
                    if (temps[2].equals("exon")) {
                        exonNumber++;
                        sb1.append(temps[3] + ",");
                        sb2.append(temps[4] + ",");
                        continue;
                    }
                    continue;
                } else if (temp.startsWith("###")) {
                    System.out.println(temp);
                    count++;
                    System.out.println(count);
                    bw.write(start + "\t" + end + "\t" + exonNumber + "\t" + sb1.toString() + "\t" + sb2.toString());
                    bw.newLine();
                    mRNA = null;
                    exonNumber = 0;
                    CDS = 1;
                    start = null;
                    end = null;
                    sb1.replace(0, sb1.length(), "");
                    sb2.replace(0, sb2.length(), "");
                    continue;
                } else {
                    continue;
                }

            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void test1() {
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            System.out.println("nohup rvtest --inVcf /data2/xiaohan/genotype_root/genotype_rootMaf005/87B18.chr" + chr + ".maf005.vcf.gz  --pheno /data1/home/xiaohan/rareallele/RVtest/phenoDir/S7.ped --geneFile /data1/home/xiaohan/rareallele/RVtest/geneFile/refegeneWheat.txt  --kernel skat,kbac --out chr" + chr + " &");
        }
    }

    public void test() {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/tall.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/";
        String Samplename = "B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E237,B18-E242,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
        String[] names = Samplename.split(",");
        HashMap nametall = new HashMap();
        String temp = null;
        String[] temps = null;
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "S7.txt").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                nametall.put(temps[0], temps[1]);
            }
            for (int i = 0; i < names.length; i++) {
                bw.write(names[i] + "\t" + nametall.get(names[i]));
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void forfun() {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/SNPsite/allSNPsite.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/SNPsite/";
        String temp = null;
        String[] temps = null;
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "SNPsite.txt").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                bw.write("SNP_" + temps[0] + "_" + temps[1] + "\t" + temps[0] + "\t" + temps[1]);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void countTable() {
        String outputDirS = "/data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/temp/";
//        String[] subDirS = {};
        String[] subDirS = {"subFastqs", "sams", "geneCount", "countTable"};
        String geneAnnotationFileS = "/data1/home/xiaohan/rareallele/SiPASpipeline/reference/ERCC92.gtf";
        List<String> nameList = new ArrayList<>();
        List<String> fileList = new ArrayList<>();
        String subCountDirS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        File[] fs = new File(subCountDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
        List<File> fList = Arrays.asList(fs);
        for (int i = 0; i < fList.size(); i++) {
            fileList.add(fList.get(i).getName().replace("Count.txt", ""));
        }
        Collections.sort(fileList);
        int geneNumber = 0;
        String geneName = null;
        ArrayList<String> geneList = new ArrayList();
        try {
            BufferedReader br = IOUtils.getTextReader(geneAnnotationFileS);
            String temp = null;
            String[] tem = null;
            while ((temp = br.readLine()) != null) {
                List<String> tList = PStringUtils.fastSplit(temp);
                tem = tList.toArray(new String[tList.size()]);
                if (tem[2].startsWith("CDS")) continue;
                if (tem[2].startsWith("exon")) {
                    String[] te = tem[8].split(";");
                    geneName = te[1].split("\"")[1].substring(0, te[1].split("\"")[1].length());
                    if (!(geneList.contains(geneName))) {
                        geneList.add(geneName);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        int[][] count = new int[geneList.size()][fList.size()];
        fList.stream().forEach(f -> {
            String temp = null;
            String[] tem = null;
            try {
                BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    List<String> tList = PStringUtils.fastSplit(temp);
                    tem = tList.toArray(new String[tList.size()]);
//                      if(tem[0].startsWith("TraesCS")){
                    if (!tem[0].startsWith("__")) {
                        if (!nameList.contains(tem[0])) {
                            nameList.add(tem[0]);
                        }
                        int index = nameList.indexOf(tem[0]);
                        count[index][fileList.indexOf(f.getName().replace("Count.txt", ""))] = Integer.parseInt(tem[1]);
                    }
                }

            } catch (Exception ex) {
                System.out.println(tem[0] + "\t1234");
                ex.printStackTrace();

            }
        });
        File subDir = new File(outputDirS, subDirS[3]);
        String outputFileS = new File(subDir, "countResult.txt").getAbsolutePath();
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputFileS).getAbsolutePath());
            sb.append("Gene" + "\t");
            for (int i = 0; i < fileList.size(); i++) {
                sb.append(fileList.get(i) + "\t");
            }
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                sb = new StringBuilder();
                for (int j = 0; j < fileList.size(); j++) {
                    if (j == 0) {
                        sb.append(nameList.get(i) + "\t");
                    }
                    sb.append(count[i][j] + "\t");
                }
                bw.write(sb.toString());
                bw.newLine();
            }

            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void readFile() {
        String infile = "/Users/yxh/Desktop/Untitled.rtf";
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        try {
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void NewFile() {
        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/hapscanner/outputDir";
        for (int i = 0; i < 42; i++) {
            int chrNumber = i + 1;
            new File(outputDirS, "output_chr" + chrNumber).mkdir();

        }
    }

    public void posAllele() {
        String inputDir = "/data2/junxu/genotype";
        String outputDir1 = "/data1/home/xiaohan/rareallele/Hapscanner/inputfile/pos";
        String outputDir2 = "/data1/home/xiaohan/rareallele/Hapscanner/inputfile/posAllele";
        String temp = null;
        String[] temps = null;
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = xiaohan.rareallele.IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split(".")[0];
            nameSet.add(Name);
            System.out.println(Name);
        }
        try {
            BufferedReader[] br = new BufferedReader[42];
            BufferedWriter[] bw1 = new BufferedWriter[42];
            BufferedWriter[] bw2 = new BufferedWriter[42];
            for (int i = 0; i < br.length; i++) {
                int chrNumber = i + 1;
                br[i] = IOUtils.getTextGzipReader(new File(inputDir, chrNumber + ".346.B18.recode.vcf.gz").getAbsolutePath());
            }
            for (int i = 0; i < br.length; i++) {
                int chrNumber = i + 1;
                bw1[i] = IOUtils.getTextWriter(new File(outputDir1, "pos_hapscanner_chr" + chrNumber + ".txt").getAbsolutePath());
                bw2[i] = IOUtils.getTextWriter(new File(outputDir2, "posAllele_hapscanner_chr" + chrNumber + ".txt").getAbsolutePath());
                bw2[i].write("Chr\tPos\tRef\tAlt(maximum 2 alternative alleles, which is seperated by \",\", e.g. A,C)");
                bw2[i].newLine();
            }
            for (int i = 0; i < br.length; i++) {
                while ((temp = br[i].readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    bw1[i].write(temps[0] + "\t" + temps[1]);
                    bw1[i].newLine();
                    bw2[i].write(temps[0] + "\t" + temps[1] + "\t" + temps[3] + "\t" + temps[4]);
                    bw2[i].newLine();
                }
            }
            for (int i = 0; i < br.length; i++) {
                br[i].close();
                bw1[i].flush();
                bw1[i].close();
                bw2[i].flush();
                bw2[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void clonefiles() {
//        String inputfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/hapscanner/parameter/parameters_hapScanner_chr1.txt";
        String inputfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/hapscanner/taxaRef/taxaRefBAM_hapscanner_chr1.txt";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/hapscanner/parameter";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/hapscanner/taxaRef";
        String temp = null;
        String prefix = "/data1/home/xiaohan/rareallele/Hapscanner/";
        BufferedReader br = IOUtils.getTextReader(inputfile);
        try {
            BufferedWriter[] bw = new BufferedWriter[41];
            for (int i = 0; i < 41; i++) {
                int chrNumber = i + 2;
                bw[i] = IOUtils.getTextWriter(new File(outputDir, "taxaRefBAM_hapscanner_chr" + chrNumber + ".txt").getAbsolutePath());
            }
            while ((temp = br.readLine()) != null) {
//                if (temp.startsWith(prefix)) {
                for (int i = 0; i < bw.length; i++) {
                    int chrNumber = i + 2;
                    System.out.print(chrNumber);
                    String temp1 = temp.replace("chr1", "chr" + chrNumber);
                    System.out.println(temp1);
                    bw[i].write(temp1);
                    bw[i].newLine();
                }
//                } else if (temp.startsWith("1")){
//                    for(int i = 0;i< bw.length;i++){
//                        System.out.println(temp);
//                        int chrNumber = i+2;
//                        String temp1 = temp.replace("1",chrNumber+"");
//                        bw[i].write(String.valueOf(temp1));
////                        bw[i].write("10");
//                        bw[i].newLine();
            }
//                } 
//                else {
//                    for (int i = 0; i < bw.length; i++) {
//                        bw[i].write(temp);
//                        bw[i].newLine();
//                    }
//                }
//            }
            br.close();
            for (int i = 0; i < bw.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void decimals() {
        double MAF = 0.0000;
        int site0 = 1;
        int sitesum = 150;
        MAF = new BigDecimal((float) site0 / sitesum).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue();
        System.out.print(MAF);

    }

    public void writecode() {
//        String[] FileNames = {"0k_1k", "1k_2k", "2k_3k", "3k_4k", "4k_5k", "5k_6k", "6k_7k", "7k_8k", "8k_9k", "9k_10k"};
//        for (int i = 0; i < FileNames.length; i++) {
//            int DiscontrolS = Integer.parseInt(FileNames[i].split("_")[0].replace("k", "")) * 1000;
//            int DiscontrolE = Integer.parseInt(FileNames[i].split("_")[1].replace("k", "")) * 1000;
//            System.out.print(DiscontrolS + "_" + DiscontrolE);
//            System.out.println("this is calculate " + DiscontrolS + "k to " + DiscontrolE + "k rare allele count……");
//        }
//        for(int i = 18;i<45;i++){
//            System.out.print("nohup vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/");
//            System.out.print(i);
//            System.out.print(".snp.vcf.gz --maf 0 --max-maf 0.05 --out ");
//            System.out.print(i+".snp.maf005 --recode && tabix -p "+i+".snp.maf005.recode.vcf.gz >log1.txt 2>&1 &");
//        }
//          for(int i = 1;i<97;i++){
//              System.out.println("/data1/home/xiaohan/rareallele/SiPASpipeline/S1-7/sams/PL-BC"+i+"Aligned.out.bam.gz");
//          }
//        String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
//        //建立对应的name位点表
//        String[] Samplenames = name.split(",");
//        for (int i = 0; i < Samplenames.length; i++) {
//            System.out.println(Samplenames[i]);
//        }
//        for (int i = 1; i < 9; i++) {
//            int chrNumber = i + 1;
//            System.out.println("nohup java -Xmx10g -jar TIGER.jar -a HapScanner -p /data1/home/xiaohan/rareallele/Hapscanner/inputfile/parameter/parameters_hapScanner_chr" + chrNumber + ".txt > log.txt &");
//            System.out.println("nohup vcftools --vcf /data1/home/xiaohan/rareallele/Hapscanner/outputDir/VCF/chr0"+i+".vcf --max-missing 0.1 --out /data1/home/xiaohan/rareallele/Hapscanner/outputDir/RNAgenotype/chr0"+i+".vcf &" );
//            System.out.println("nohup vcftools --gzvcf /data2/junxu/genotype/"+i+".346.B18.recode.vcf.gz --maf 0 --max-maf 0.05 --out /data2/xiaohan/genotypeMaf005/346B18.chr00"+i+".maf005 --recode && bgzip /data2/xiaohan/genotypeMaf005/346B18.chr00"+i+".maf005.recode.vcf &");
        HashSet p1 = new HashSet();
        HashSet p2 = new HashSet();
        p1.add("qitiandasheng");
        p1.add("sunwukong");
        p2.add("sunwukong");
        p2.add("bimawen");
        HashSet set = new HashSet();
        set.add(p1);
        p1.retainAll(p2);
        System.out.println(p1.isEmpty());
        System.out.print(p1);
        System.out.print(set);
//        }

    }

    public void mkPosGeneMap() {
        String geneNameS = null;
        int gfIndex = 0;
        ArrayList<String> geneNameList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/refer/wheat_v1.1_Lulab.gtf");
            String temp = null;
            Set<String> geneSet = new HashSet<String>();
            Set<Integer> chrSet = new HashSet<>();
            String[] tem = null;
            String geneName = null;
            HashMap<String, Integer> geneChrMap = new HashMap();
            HashMap<String, Integer> geneMinMap = new HashMap();
            HashMap<String, Integer> geneMaxMap = new HashMap();
            HashMap<String, Byte> geneStrandMap = new HashMap();
            while ((temp = br.readLine()) != null) {
                tem = temp.split("\t");
                if (!tem[2].startsWith("exon")) continue;
                String[] te = tem[8].split(";");
                geneName = te[1].split("\"")[1];
                if (!geneSet.contains(geneName)) {
                    geneMinMap.put(geneName, Integer.MAX_VALUE);
                    geneMaxMap.put(geneName, Integer.MIN_VALUE);
                }
                geneSet.add(geneName);
                int min = Integer.parseInt(tem[3]);
                int max = Integer.parseInt(tem[4]);
                if (geneMinMap.get(geneName) > min) geneMinMap.put(geneName, min);
                if (geneMaxMap.get(geneName) < max) geneMaxMap.put(geneName, max);
                int chr = Integer.parseInt(tem[0]);
                chrSet.add(chr);
                geneChrMap.put(geneName, chr);
                if (tem[6].startsWith("-")) geneStrandMap.put(geneName, (byte) 1);
                else geneStrandMap.put(geneName, (byte) 0);
            }

            for (String s : geneSet) {
                geneNameList.add(s);
            }
            Collections.sort(geneNameList);
            for (String s : geneNameList) {
                System.out.println(geneChrMap.get(s) + "\t" + geneMinMap.get(s) + "\t" + geneMaxMap.get(s) + "\t" + s);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            System.out.println(geneNameS + "\t" + gfIndex);
        }

    }

    public void findNumber() throws IOException {
        String infile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/ZX-ZX.txt";
        HashMap<String, String> trans = new HashMap<String, String>();
        String table = "ZX393\tZX1497\n" +
                "ZX394\tZX1498\n" +
                "ZX395\tZX1499\n" +
                "ZX396\tZX1500\n" +
                "ZX397\tZX1501\n" +
                "ZX398\tZX1502\n" +
                "ZX399\tZX1503\n" +
                "ZX400\tZX1504\n" +
                "ZX401\tZX1505\n" +
                "ZX402\tZX1506\n" +
                "ZX403\tZX1507\n" +
                "ZX404\tZX1508\n" +
                "ZX405\tZX1509\n" +
                "ZX406\tZX1510\n" +
                "ZX407\tZX1511\n" +
                "ZX408\tZX1512\n" +
                "ZX409\tZX1513\n" +
                "ZX410\tZX1514\n" +
                "ZX411\tZX1515\n" +
                "ZX412\tZX1516\n" +
                "ZX413\tZX1517\n" +
                "ZX414\tZX1518\n" +
                "ZX415\tZX1519\n" +
                "ZX416\tZX1520\n" +
                "ZX417\tZX1521\n" +
                "ZX418\tZX1522\n" +
                "ZX419\tZX1523\n" +
                "ZX420\tZX1524\n" +
                "ZX421\tZX1525\n" +
                "ZX422\tZX1526\n" +
                "ZX423\tZX1527\n" +
                "ZX424\tZX1528\n" +
                "ZX425\tZX1529\n" +
                "ZX426\tZX1530\n" +
                "ZX427\tZX1531\n" +
                "ZX428\tZX1532\n" +
                "ZX429\tZX1533\n" +
                "ZX430\tZX1534\n" +
                "ZX431\tZX1535\n" +
                "ZX432\tZX1536\n" +
                "ZX433\tZX1537\n" +
                "ZX434\tZX1538\n" +
                "ZX435\tZX1539\n" +
                "ZX436\tZX1540\n" +
                "ZX437\tZX1541\n" +
                "ZX438\tZX1542\n" +
                "ZX439\tZX1543\n" +
                "ZX440\tZX1544\n" +
                "ZX441\tZX1545\n" +
                "ZX442\tZX1546\n" +
                "ZX443\tZX1547\n" +
                "ZX444\tZX1548\n" +
                "ZX445\tZX1549\n" +
                "ZX446\tZX1550\n" +
                "ZX447\tZX1551\n" +
                "ZX448\tZX1552\n" +
                "ZX449\tZX1553\n" +
                "ZX450\tZX1554\n" +
                "ZX451\tZX1555\n" +
                "ZX452\tZX1556\n" +
                "ZX453\tZX1557\n" +
                "ZX454\tZX1558\n" +
                "ZX455\tZX1559\n" +
                "ZX456\tZX1560\n" +
                "ZX457\tZX1561\n" +
                "ZX458\tZX1562\n" +
                "ZX459\tZX1563\n" +
                "ZX460\tZX1564\n" +
                "ZX461\tZX1565\n" +
                "ZX462\tZX1566\n" +
                "ZX463\tZX1567\n" +
                "ZX464\tZX1568\n" +
                "ZX465\tZX1569\n" +
                "ZX466\tZX1570\n" +
                "ZX467\tZX1571\n" +
                "ZX468\tZX1572\n" +
                "ZX469\tZX1573\n" +
                "ZX470\tZX1574\n" +
                "ZX471\tZX1575\n" +
                "ZX472\tZX1576\n" +
                "ZX473\tZX1577\n" +
                "ZX474\tZX1578\n" +
                "ZX475\tZX1579\n" +
                "ZX476\tZX1580\n" +
                "ZX477\tZX1581\n" +
                "ZX478\tZX1582\n" +
                "ZX479\tZX1583\n" +
                "ZX480\tZX1584\n" +
                "ZX481\tZX1585\n" +
                "ZX482\tZX1586\n" +
                "ZX483\tZX1587\n" +
                "ZX484\tZX1588\n" +
                "ZX485\tZX1589\n" +
                "ZX486\tZX1590\n" +
                "ZX487\tZX1591\n" +
                "ZX488\tZX1592\n" +
                "ZX489\tZX1593\n" +
                "ZX490\tZX1594\n" +
                "ZX491\tZX1595\n" +
                "ZX492\tZX1596\n" +
                "ZX493\tZX1597\n" +
                "ZX494\tZX1598\n" +
                "ZX495\tZX1599\n" +
                "ZX496\tZX1600\n" +
                "ZX497\tZX1601\n" +
                "ZX498\tZX1602\n" +
                "ZX499\tZX1603\n" +
                "ZX500\tZX1604\n" +
                "ZX501\tZX1605\n" +
                "ZX502\tZX1606\n" +
                "ZX503\tZX1607\n" +
                "ZX504\tZX1608\n" +
                "ZX505\tZX1609\n" +
                "ZX506\tZX1610\n" +
                "ZX507\tZX1611\n" +
                "ZX508\tZX1612\n" +
                "ZX509\tZX1613\n" +
                "ZX510\tZX1614\n" +
                "ZX511\tZX1615\n" +
                "ZX512\tZX1616\n" +
                "ZX513\tZX1617\n" +
                "ZX514\tZX1618\n" +
                "ZX515\tZX1619\n" +
                "ZX516\tZX1620\n" +
                "ZX517\tZX1621\n" +
                "ZX518\tZX1622\n" +
                "ZX519\tZX1623\n" +
                "ZX520\tZX1624\n" +
                "ZX521\tZX1625\n" +
                "ZX522\tZX1626\n" +
                "ZX523\tZX1627\n" +
                "ZX524\tZX1628\n" +
                "ZX525\tZX1629\n" +
                "ZX526\tZX1630\n" +
                "ZX527\tZX1631\n" +
                "ZX528\tZX1632\n" +
                "ZX529\tZX1633\n" +
                "ZX530\tZX1634\n" +
                "ZX531\tZX1635\n" +
                "ZX532\tZX1636\n" +
                "ZX533\tZX1637\n" +
                "ZX534\tZX1638\n" +
                "ZX535\tZX1639\n" +
                "ZX536\tZX1640\n" +
                "ZX537\tZX1641\n" +
                "ZX538\tZX1642\n" +
                "ZX539\tZX1643\n" +
                "ZX540\tZX1644\n" +
                "ZX541\tZX1645\n" +
                "ZX542\tZX1646\n" +
                "ZX543\tZX1647\n" +
                "ZX544\tZX1648\n" +
                "ZX545\tZX1649\n" +
                "ZX546\tZX1650\n" +
                "ZX547\tZX1651\n" +
                "ZX548\tZX1652\n" +
                "ZX549\tZX1653\n" +
                "ZX550\tZX1654\n" +
                "ZX551\tZX1655\n" +
                "ZX552\tZX1656\n" +
                "ZX553\tZX1657\n" +
                "ZX554\tZX1658\n" +
                "ZX555\tZX1659\n" +
                "ZX556\tZX1660\n" +
                "ZX557\tZX1661\n" +
                "ZX558\tZX1662\n" +
                "ZX559\tZX1663\n" +
                "ZX560\tZX1664\n" +
                "ZX561\tZX1665\n" +
                "ZX562\tZX1666\n" +
                "ZX563\tZX1667\n" +
                "ZX564\tZX1668\n" +
                "ZX565\tZX1669\n" +
                "ZX566\tZX1670\n" +
                "ZX567\tZX1671\n" +
                "ZX568\tZX1672\n" +
                "ZX569\tZX1673\n" +
                "ZX570\tZX1674\n" +
                "ZX571\tZX1675\n" +
                "ZX572\tZX1676\n" +
                "ZX573\tZX1677\n" +
                "ZX574\tZX1678\n" +
                "ZX575\tZX1679\n" +
                "ZX576\tZX1680\n" +
                "ZX577\tZX1681\n" +
                "ZX578\tZX1682\n" +
                "ZX579\tZX1683\n" +
                "ZX580\tZX1684\n" +
                "ZX581\tZX1685\n" +
                "ZX582\tZX1686\n" +
                "ZX583\tZX1687\n" +
                "ZX584\tZX1688\n" +
                "ZX585\tZX1689\n" +
                "ZX586\tZX1690\n" +
                "ZX587\tZX1691\n" +
                "ZX588\tZX1692\n" +
                "ZX589\tZX1693\n" +
                "ZX590\tZX1694\n" +
                "ZX591\tZX1695\n" +
                "ZX592\tZX1696\n" +
                "ZX593\tZX1697\n" +
                "ZX594\tZX1698\n" +
                "ZX595\tZX1699\n" +
                "ZX596\tZX1700\n" +
                "ZX597\tZX1701\n" +
                "ZX598\tZX1702\n" +
                "ZX599\tZX1703\n" +
                "ZX600\tZX1704\n" +
                "ZX601\tZX1705\n" +
                "ZX602\tZX1706\n" +
                "ZX603\tZX1707\n" +
                "ZX604\tZX1708\n" +
                "ZX605\tZX1709\n" +
                "ZX606\tZX1710\n" +
                "ZX607\tZX1711\n" +
                "ZX608\tZX1712\n" +
                "ZX609\tZX1713\n" +
                "ZX610\tZX1714\n" +
                "ZX611\tZX1715\n" +
                "ZX612\tZX1716\n" +
                "ZX613\tZX1717\n" +
                "ZX614\tZX1718\n" +
                "ZX615\tZX1719\n" +
                "ZX616\tZX1720\n" +
                "ZX617\tZX1721\n" +
                "ZX618\tZX1722\n" +
                "ZX619\tZX1723\n" +
                "ZX620\tZX1724\n" +
                "ZX621\tZX1725\n" +
                "ZX622\tZX1726\n" +
                "ZX623\tZX1727\n" +
                "ZX624\tZX1728\n" +
                "ZX625\tZX1729\n" +
                "ZX626\tZX1730\n" +
                "ZX627\tZX1731\n" +
                "ZX628\tZX1732\n" +
                "ZX629\tZX1733\n" +
                "ZX630\tZX1734\n" +
                "ZX631\tZX1735\n" +
                "ZX632\tZX1736\n" +
                "ZX633\tZX1737\n" +
                "ZX634\tZX1738\n" +
                "ZX635\tZX1739\n" +
                "ZX636\tZX1740\n" +
                "ZX637\tZX1741\n" +
                "ZX638\tZX1742\n" +
                "ZX639\tZX1743\n" +
                "ZX640\tZX1744\n" +
                "ZX641\tZX1745\n" +
                "ZX642\tZX1746\n" +
                "ZX643\tZX1747\n" +
                "ZX644\tZX1748\n" +
                "ZX645\tZX1749\n" +
                "ZX646\tZX1750\n" +
                "ZX647\tZX1751\n" +
                "ZX648\tZX1752\n" +
                "ZX649\tZX1753\n" +
                "ZX650\tZX1754\n" +
                "ZX651\tZX1755\n" +
                "ZX652\tZX1756\n" +
                "ZX653\tZX1757\n" +
                "ZX654\tZX1758\n" +
                "ZX655\tZX1759\n" +
                "ZX656\tZX1760\n" +
                "ZX657\tZX1761\n" +
                "ZX658\tZX1762\n" +
                "ZX659\tZX1763\n" +
                "ZX660\tZX1764\n" +
                "ZX661\tZX1765\n" +
                "ZX662\tZX1766\n" +
                "ZX663\tZX1767\n" +
                "ZX664\tZX1768\n" +
                "ZX665\tZX1769\n" +
                "ZX666\tZX1770\n" +
                "ZX667\tZX1771\n" +
                "ZX668\tZX1772\n" +
                "ZX669\tZX1773\n" +
                "ZX670\tZX1774\n" +
                "ZX671\tZX1775\n" +
                "ZX672\tZX1776\n" +
                "ZX673\tZX1777\n" +
                "ZX674\tZX1778\n" +
                "ZX675\tZX1779\n" +
                "ZX676\tZX1780\n" +
                "ZX677\tZX1781\n" +
                "ZX678\tZX1782\n" +
                "ZX679\tZX1783\n" +
                "ZX680\tZX1784\n" +
                "ZX681\tZX1785\n" +
                "ZX682\tZX1786\n" +
                "ZX683\tZX1787\n" +
                "ZX684\tZX1788\n" +
                "ZX685\tZX1789\n" +
                "ZX686\tZX1790\n" +
                "ZX687\tZX1791\n" +
                "ZX688\tZX1792\n" +
                "ZX689\tZX1793\n" +
                "ZX690\tZX1794\n" +
                "ZX691\tZX1795\n" +
                "ZX692\tZX1796\n" +
                "ZX693\tZX1797\n" +
                "ZX694\tZX1798\n" +
                "ZX695\tZX1799\n" +
                "ZX696\tZX1800\n" +
                "ZX697\tZX1801\n" +
                "ZX698\tZX1802\n" +
                "ZX699\tZX1803\n" +
                "ZX700\tZX1804\n" +
                "ZX701\tZX1805\n" +
                "ZX702\tZX1806\n" +
                "ZX703\tZX1807\n" +
                "ZX704\tZX1808\n" +
                "ZX705\tZX1809\n" +
                "ZX706\tZX1810\n" +
                "ZX707\tZX1811\n" +
                "ZX708\tZX1812\n" +
                "ZX709\tZX1813\n" +
                "ZX710\tZX1814\n" +
                "ZX711\tZX1815\n" +
                "ZX712\tZX1816\n" +
                "ZX713\tZX1817\n" +
                "ZX714\tZX1818\n" +
                "ZX715\tZX1819\n" +
                "ZX716\tZX1820\n" +
                "ZX717\tZX1821\n" +
                "ZX718\tZX1822\n" +
                "ZX719\tZX1823\n" +
                "ZX720\tZX1824\n" +
                "ZX721\tZX1825\n" +
                "ZX722\tZX1826\n" +
                "ZX723\tZX1827\n" +
                "ZX724\tZX1828\n" +
                "ZX725\tZX1829\n" +
                "ZX726\tZX1830\n" +
                "ZX727\tZX1831\n" +
                "ZX728\tZX1832\n" +
                "ZX729\tZX1833\n" +
                "ZX730\tZX1834\n" +
                "ZX731\tZX1835\n" +
                "ZX732\tZX1836\n" +
                "ZX733\tZX1837\n" +
                "ZX734\tZX1838\n" +
                "ZX735\tZX1839\n" +
                "ZX736\tZX1840\n" +
                "ZX737\tZX1841\n" +
                "ZX738\tZX1842\n" +
                "ZX739\tZX1843\n" +
                "ZX740\tZX1844\n" +
                "ZX741\tZX1845\n" +
                "ZX742\tZX1846\n" +
                "ZX743\tZX1847\n" +
                "ZX744\tZX1848\n" +
                "ZX745\tZX1849\n" +
                "ZX746\tZX1850\n" +
                "ZX747\tZX1851\n" +
                "ZX748\tZX1852\n" +
                "ZX749\tZX1853\n" +
                "ZX750\tZX1854\n" +
                "ZX751\tZX1855\n" +
                "ZX752\tZX1856\n" +
                "ZX753\tZX1857\n" +
                "ZX754\tZX1858\n" +
                "ZX755\tZX1859\n" +
                "ZX756\tZX1860\n" +
                "ZX757\tZX1861\n" +
                "ZX758\tZX1862\n" +
                "ZX759\tZX1863\n" +
                "ZX760\tZX1864\n" +
                "ZX761\tZX1865\n" +
                "ZX762\tZX1866\n" +
                "ZX763\tZX1867\n" +
                "ZX764\tZX1868\n" +
                "ZX765\tZX1869\n" +
                "ZX766\tZX1870\n" +
                "ZX767\tZX1871\n" +
                "ZX768\tZX1872\n" +
                "ZX769\tZX1873\n" +
                "ZX770\tZX1874\n" +
                "ZX771\tZX1875\n" +
                "ZX772\tZX1876\n" +
                "ZX773\tZX1877\n" +
                "ZX774\tZX1878\n" +
                "ZX775\tZX1879\n" +
                "ZX776\tZX1880\n" +
                "ZX777\tZX1881\n" +
                "ZX778\tZX1882\n" +
                "ZX779\tZX1883\n" +
                "ZX780\tZX1884\n" +
                "ZX781\tZX1885\n" +
                "ZX782\tZX1886";
        String[] temps = table.split("\n");
        String[] temp = null;
        for (int i = 1; i < temps.length; i++) {
            temp = temps[i].split("\t");
            trans.put(temp[0], temp[1]);
        }
        try {
            BufferedReader br = IOUtils.getTextReader(infile1);
            for (int i = 1; i <= 96; i++) {
                String count = br.readLine();
                System.out.println(count + "\t");
                System.out.println(trans.get(count) + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    public static void main(String args[]) throws IOException {
        new test();
    }
}

