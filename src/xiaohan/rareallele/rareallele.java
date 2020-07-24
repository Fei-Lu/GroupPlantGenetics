package xiaohan.rareallele;

import pgl.infra.table.RowTable;
import pgl.infra.utils.PStringUtils;

import java.io.*;
import java.util.*;

import static java.lang.Integer.parseInt;

public class rareallele {
    public rareallele() {
//        String infileS ="snp1";
//        String inputDir = "/data2/xiaohan/sub/";
//        String outputDir = "/data2/xiaohan/GT/";
//        String infileS = "36.snp.maf005.recode";
//        String inputDir = "/data2/xiaohan/SNP/";
//        String outputDir = "/data1/home/xiaohan/rareallele/fastQTL/";
//        String infileS = "vcf";
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/";
//        this.GetCandidateVariants();
//        this.parseFqByIndex();
//        this.read10lines();
//        this.ifanyindex();
//        this.parseFqByIndexAndBarcode();
        //this.getExampleVCF();
        //this.rankGenes();
//        this.getTransNumber();
//        this.getVCFposition(infileS, inputDir);
//        this.getsubVCF(infileS,inputDir,outputDir);
        //this.checklines();
        //this.countlines();
        //this.test();
//        this.changeSampleName();
//        this.getGTvcf2();
//        this.getGTvcf(infileS, outputDir);
//        this.changeName();
//        this.get5kSNPcount();
        this.getupstreamSNPcount();
//        this.expressionRank();
//        this.SNPcount();
//        this.SNPTable();
//        this.rankcorrelation();
//        this.CalculateeGenes();
//        this.SplitPhenoBychr();
//        this.code();
//        this.countTable();
//        this.is3M();
//        this.candidateSNP();
//        this.candidate();
//        this.writesetFile();
//        this.extractTop();
        this.extractHomologousGene();
//        this.extractsubgenome();
//        this.extractTopPheno();
    }

    public void extractTopPheno() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3";
        String infileS = "countResult3_DESeq2.txt";
//        String outfileS = "Top5000_S3.txt";
        String outfileS = "Top5000-10000_S3.txt";
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7";
//        String infileS = "countResult7_DESeq2.txt";
//        String outfileS = "Top5000_S7.txt";
//        String outfileS = "Top5000-10000_S7.txt";
//        String infor = "Top_5000_median_exp_genes.txt";
        String infor = "Top_5001_to_10000_median_exp_genes.txt";
        try {
            String temp = null;
            String[] temps = null;
            HashSet GeneSet = new HashSet();
            HashMap<String, Integer> geneMap = new HashMap();
            BufferedReader br = IOUtils.getTextReader(new File(inputDir, infor).getAbsolutePath());
            System.out.println(new File(inputDir, infor).getAbsolutePath());
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("x")) continue;
                temps = temp.split("\t");
                System.out.println(temps[1]);
                GeneSet.add(temps[1]);
                geneMap.put(temps[1], Integer.parseInt(temps[0]));
                continue;
            }
            String temp1 = null;
            String[] temps1 = null;
            String[] rank = new String[5000];
            BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, infileS).getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, outfileS).getAbsolutePath());
            while ((temp1 = br1.readLine()) != null) {
                temps1 = temp1.split("\t");
                if (temp1.startsWith("g")) {
                    bw.write(temp1);
                    bw.newLine();
                }
                if (GeneSet.contains(temps1[0])) {
                    String gene = temps1[0];
                    int index = geneMap.get(gene);
                    rank[index - 1] = temp1;
                }
            }
            for (int i = 0; i < rank.length; i++) {
                bw.write(rank[i]);
                bw.newLine();
            }
            br1.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void extractTop() {
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
//        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/expressionTable";
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable";
        String infor = "Top_5000_median_exp_genes.txt";
//        String infor = "Top_5001_to_10000_median_exp_genes.txt";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                HashSet GeneSet = new HashSet();
                HashMap<String, Integer> geneMap = new HashMap();
                BufferedReader br = IOUtils.getTextReader(new File(inforDir, infor).getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("x")) continue;
                    temps = temp.split("\t");
                    GeneSet.add(temps[1]);
                    geneMap.put(temps[1], Integer.parseInt(temps[0]));
                    continue;
                }
                String temp1 = null;
                String[] temps1 = null;
                String[] rank = new String[5000];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedReader br2 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
//                    BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, "Top5000_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw = IOUtils.getTextWriter(new File(inputDir, "Top5000-10000_" +f+ "_count.txt").getAbsolutePath());
                    while ((temp1 = br2.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (temps1[0].startsWith("gene")) {
                            bw.write(temp1);
                            bw.newLine();
                        }
                    }
                    while ((temp1 = br1.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (GeneSet.contains(temps1[0])) {
                            String gene = temps1[0];
                            int index = geneMap.get(gene);
                            rank[index - 1] = temp1;
                        }
                    }
                    for (int i = 0; i < rank.length; i++) {
                        bw.write(rank[i]);
                        bw.newLine();
                    }
                    br1.close();
                    br2.close();
                    bw.flush();
                    bw.close();
                    br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }

    public void extractHomologousGene() {
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data";
        String infor = "TheABD.txt";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                HashSet GeneSetA = new HashSet();
                HashSet GeneSetB = new HashSet();
                HashSet GeneSetD = new HashSet();
                HashMap<String, Integer> genesubAMap = new HashMap();
                HashMap<String, Integer> genesubBMap = new HashMap();
                HashMap<String, Integer> genesubDMap = new HashMap();
                BufferedReader br = IOUtils.getTextReader(new File(inforDir, infor).getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    GeneSetA.add(temps[1]);
                    GeneSetB.add(temps[2]);
                    GeneSetD.add(temps[3]);
                    genesubAMap.put(temps[1], Integer.parseInt(temps[0]));
                    genesubBMap.put(temps[2], Integer.parseInt(temps[0]));
                    genesubDMap.put(temps[3], Integer.parseInt(temps[0]));
                    continue;
                }
                String temp1 = null;
                String[] temps1 = null;
                String[] rankA = new String[17108];
                String[] rankB = new String[17108];
                String[] rankD = new String[17108];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedReader br2 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw1 = IOUtils.getTextWriter(new File(inputDir, "Ahomo_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw2 = IOUtils.getTextWriter(new File(inputDir, "Bhomo_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw3 = IOUtils.getTextWriter(new File(inputDir, "Dhomo_" +f+ "_count.txt").getAbsolutePath());
                    while ((temp1 = br2.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (temps1[0].startsWith("gene")) {
                            bw1.write(temp1);
                            bw1.newLine();
                            bw2.write(temp1);
                            bw2.newLine();
                            bw3.write(temp1);
                            bw3.newLine();
                        }
                    }
                    while ((temp1 = br1.readLine()) != null) {
                        temps1 = temp1.split("\t");
                        if (GeneSetA.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexA = genesubAMap.get(gene);
                            rankA[indexA - 1] = temp1;
                        }
                        if (GeneSetB.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexB = genesubBMap.get(gene);
                            rankB[indexB - 1] = temp1;
                        }
                        if (GeneSetD.contains(temps1[0])) {
                            String gene = temps1[0];
                            int indexD = genesubDMap.get(gene);
                            rankD[indexD - 1] = temp1;
                        }
                    }
                    for (int i = 0; i < rankA.length; i++) {
                        bw1.write(rankA[i]);
                        bw1.newLine();
                    }
                    for (int i = 0; i < rankB.length; i++) {
                        bw2.write(rankB[i]);
                        bw2.newLine();
                    }
                    for (int i = 0; i < rankD.length; i++) {
                        bw3.write(rankD[i]);
                        bw3.newLine();
                    }
                    br1.close();
                    br2.close();
                    bw1.flush();bw2.flush();bw3.flush();
                    bw1.close();bw2.close();bw3.close();
                    br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
    
    public void extractsubgenome() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/snp_count/major_SNP_site";
//        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count/major_SNP_site";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2] ;
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                String[] rank = new String[5000];
                    BufferedReader br1 = IOUtils.getTextReader(new File(inputDir, "all_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw1 = IOUtils.getTextWriter(new File(inputDir, "A_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw2 = IOUtils.getTextWriter(new File(inputDir, "B_" +f+ "_count.txt").getAbsolutePath());
                    BufferedWriter bw3 = IOUtils.getTextWriter(new File(inputDir, "D_" +f+ "_count.txt").getAbsolutePath());
                    temp = br1.readLine();
                    bw1.write(temp);
                    bw2.write(temp);
                    bw3.write(temp);
                    bw1.newLine();
                    bw2.newLine();
                    bw3.newLine();
                    while((temp = br1.readLine())!= null){
                        if(!temp.startsWith("T"))continue;
                        temps = temp.split("\t");
                        String geneName = temps[0];
                        String subgenome = geneName.toString().substring(8, 9);
                        if(subgenome.equals("A")){
                            bw1.write(temp);
                            bw1.newLine();
                        }
                        if(subgenome.equals("B")){
                            bw2.write(temp);
                            bw2.newLine();
                        }
                        if(subgenome.equals("D")){
                            bw3.write(temp);
                            bw3.newLine();
                        }
                    }
                    br1.close();
                    bw1.flush();bw1.close();
                    bw2.flush();bw2.close();
                    bw3.flush();bw3.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void writesetFile() {
        String inputDir = "/data1/home/xiaohan/rareallele/rvtest/output/different";
        String outputDir = "/data1/home/xiaohan/rareallele/rvtest/output/different";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".SingleWald.assoc");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            System.out.println(fs[i]);
            String Name = fs[i].getName().split("\\.")[0];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                int start = 0;
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, f + ".SingleWald.assoc").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, f + "correlation").getAbsolutePath());
                if (f.endsWith("A")) {
                    start = 268932112;
                }
                if (f.endsWith("B")) {
                    start = 405524690;
                }
                if (f.endsWith("D")) {
                    start = 62563846;
                }
                String temp = null;
                String[] temps = null;
                int begin = start - 1000000;
                bw.write("CHR" + "\t" + "BP" + "\t" + "SNP" + "\t" + "P");
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("C")) continue;
                    temps = temp.split("\t");
                    int position = Integer.parseInt(temps[1]);
                    if (position > begin && position < start) {
                        bw.write(temps[0] + "\t" + position + "\t" + "snp_" + position + "\t" + temps[8]);
                        bw.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    public void candidate() {
        String input = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S7/output";
        String output = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S7/output/e";
        String infor = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression";
        File[] fs = new File(input).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".Skat.assoc");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            System.out.println(fs[i].getName().split("\\.")[0]);
            String Name = fs[i].getName().split("\\.")[0].split("t")[1];
            nameSet.add(Name);
//            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(input, "out" + f + ".Skat.assoc").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(output, "chr" + f + "_Skat_assoc.txt").getAbsolutePath());
//                BufferedWriter bw1 = IOUtils.getTextWriter(new File(output,"chr"+f+"count.txt").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextReader(new File(infor, "S7expression" + f + ".bed").getAbsolutePath());
                String temp1 = null;
                String[] temps1 = null;
                String temp = null;
                String[] temps = null;
                HashMap<String, Integer> StartMap = new HashMap<>();
                while ((temp1 = br1.readLine()) != null) {
                    if (temp1.startsWith("#")) continue;
                    temps1 = temp1.split("\t");
                    StartMap.put(temps1[3], Integer.parseInt(temps1[1]));
                }
                bw.write("CHR" + "\t" + "BP" + "\t" + "SNP" + "\t" + "P");
                bw.newLine();
//                bw1.write("Distance" +"\t" + "Count");
//                bw1.newLine();
                int[] countTable = new int[1000];
                for (int i = 0; i < 1000; i++) {
                    countTable[i] = 0;
                }
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("R")) continue;
                    temps = temp.split("\t");
                    if (temps[6].equals("NA")) continue;
                    String geneName = temps[0].split("_")[0];
                    int start = StartMap.get(geneName);
                    String position1 = temps[1].split("-")[0].split(":")[1];
                    String position2 = temps[1].split("-")[1];
                    int position = (Integer.parseInt(position1) + Integer.parseInt(position2)) / 2;
                    int distance = start - position;
                    bw.write(f + "\t");
                    bw.write(position + "\t" + "snp_" + position + "\t" + temps[6]);
                    bw.newLine();
//                    System.out.println(Double.parseDouble(temps[6]));
//                    if(Double.parseDouble(temps[6])<0.01){
//                        int distanceNumber = distance / 1000;
//                        countTable[distanceNumber] = countTable[distanceNumber] + Integer.parseInt(temps[3]);
//                    }
                }
//                for(int i = 0 ; i <1000;i++){
//                    bw1.write(i + "\t" + countTable[i]);
//                    bw1.newLine();
//                }
//                bw1.flush();bw1.close();
                br.close();
                br1.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void candidateSNP() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S3/output";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S3/output/e";
        String infor = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/splitexpression";
//        String inputDir =  "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S7/output" ;
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/S7/output/candidate" ;
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".SingleWald.assoc");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            System.out.println(fs[i].getName().split("\\.")[0]);
            String Name = fs[i].getName().split("\\.")[0].split("t")[1];
            nameSet.add(Name);
//            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br1 = IOUtils.getTextReader(new File(infor, "S3expression" + f + ".bed").getAbsolutePath());
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, "out" + f + ".SingleWald.assoc").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "chr" + f + "_single_assoc.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDir, "chr" + f + "single.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                HashMap<String, Integer> StartMap = new HashMap<>();
                while ((temp1 = br1.readLine()) != null) {
                    if (temp1.startsWith("#")) continue;
                    temps1 = temp1.split("\t");
                    StartMap.put(temps1[3], Integer.parseInt(temps1[1]));
                }
                bw.write("CHR" + "\t" + "BP" + "\t" + "SNP" + "\t" + "P" + "\t" + "Beta");
                bw.newLine();
                bw1.write("Distance" + "Effect");
                bw1.newLine();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("R")) continue;
                    temps = temp.split("\t");
                    if (temps[9].equals("NA") || temps[7].equals("NA")) continue;
                    bw.write(temps[1] + "\t" + temps[2] + "\t");
                    String SNPname = "snp" + "_" + temps[1] + "_" + temps[2] + "_" + temps[0];
                    bw.write(SNPname + "\t");
                    bw.write(temps[9] + "\t");
                    bw.write(temps[7]);
                    bw.newLine();
                    double pvalue = Double.parseDouble(temps[9]);
                    System.out.println(pvalue);
                    if (pvalue < 0.01) {
                        String geneName = temps[0].split("_")[0];
                        int start = StartMap.get(geneName);
                        System.out.println(start);
                        int distance = start - Integer.parseInt(temps[2]);
                        bw1.write(distance + "\t" + temps[7]);
                        bw1.newLine();
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    public void is3M() {
        String infile1 = "/Users/yxh/Documents/eQTL/009ERCC test/SiPAS-Truseq_1M/SiPAS3M/logamR1.txt";
        String infile2 = "/Users/yxh/Documents/eQTL/009ERCC test/SiPAS-Truseq_1M/SiPAS3M/logpmR1.txt";
        BufferedReader br1 = IOUtils.getTextReader(infile1);
        BufferedReader br2 = IOUtils.getTextReader(infile2);
        String temp1 = null;
        String temp2 = null;
        try {
            int count = 1;
            while ((temp1 = br1.readLine()) != null) {
                temp2 = br2.readLine();
                if (Integer.valueOf(temp1).equals(Integer.valueOf(temp2))) {
                    System.out.println(count);
                }
                count++;
            }
            br1.close();
            br2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void countTable() {
        String outputDir = "/Users/yxh/Documents/eQTL/009ERCC test/SiPAS-Truseq_1M/SiPAS3M";
        List<String> plateList = new ArrayList<>();
        List<File>[] fileList = new List[96];//
        for (int i = 0; i < 96; i++) {//
            fileList[i] = new ArrayList<File>();
        }
//        String subCountDirS = new File ("/Users/xujun/Desktop/eQTL/total/geneCount").getAbsolutePath();
        String subCountDirS = new File(outputDir, "geneCount").getAbsolutePath();
        File[] fs = new File(subCountDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "Count.txt");
        List<File> fList = Arrays.asList(fs);
        for (int i = 0; i < fList.size(); i++) {
            String plate = fList.get(i).getName().split("_")[1];
            if (!plateList.contains(plate)) {
                plateList.add(plate);
            }
            fileList[plateList.indexOf(plate)].add(fList.get(i));
        }
        int geneNumber = 0;
        String geneName = null;
        Set<String> geneSet = new HashSet<String>();
        for (int i = 0; i < fileList.length; i++) {
            Collections.sort(fileList[i]);
        }
        int geneNmuber = 0;
        StringBuilder wc = new StringBuilder();
        wc.append("wc -l ").append(fileList[0].get(0));
        String command = wc.toString();
        System.out.println(command);
        try {
            Runtime rt = Runtime.getRuntime();
            Process p = rt.exec(command);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                geneNumber = Integer.valueOf(temp.split(" ")[0]) - 5;
            }
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < fileList.length; i++) {
            Set<String> nameSet = new HashSet<String>();
            List<String> nameList = new ArrayList<>();
            int[][] count = new int[geneNumber][fileList[i].size()];
            List<File> list = fileList[i];
            list.stream().forEach(f -> {
                String temp = null;
                String[] tem = null;
                try {
                    BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
                    int rowN = 0;
                    while ((temp = br.readLine()) != null) {
                        List<String> tList = PStringUtils.fastSplit(temp);
                        tem = tList.toArray(new String[tList.size()]);
                        if (!tem[0].startsWith("__")) {
                            if (!nameSet.contains(tem[0])) {
                                nameList.add(tem[0]);
                            }
                            nameSet.add(tem[0]);
//                          int index1=nameList.indexOf(tem[0]);
                            count[rowN][list.indexOf(f)] = Integer.parseInt(tem[1]);
                            rowN++;
                        }
                    }

                } catch (Exception ex) {
                    System.out.println(tem[0] + "\t" + geneSet.size() + "\t1234");
                    ex.printStackTrace();

                }
            });
            File subDir = new File(outputDir, "countTable");
            String outputFileS = new File(subDir, fileList[i].get(0).getName().split("_")[1] + "_countResult.txt").getAbsolutePath();
            try {
                StringBuilder sb = new StringBuilder();
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputFileS).getAbsolutePath());
                sb.append("Gene" + "\t");
                for (int j = 0; j < fileList[i].size(); j++) {
                    sb.append(fileList[i].get(j).getName().replaceAll("_Count.txt", "") + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
                for (int k = 0; k < count.length; k++) {
                    sb = new StringBuilder();
                    for (int j = 0; j < fileList[i].size(); j++) {
                        if (j == 0) {
                            sb.append(nameList.get(k) + "\t");
                        }
                        sb.append(count[k][j] + "\t");
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                }

                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }

    public void code() {

        for (int i = 1; i < 49; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append("nohup /data1/home/xiaohan/myprogram/STAR-2.6.0c/bin/Linux_x86_64/STAR --runThreadN 160 --genomeDir /data1/home/xiaohan/rareallele/SiPASpipeline/ERCC/ampm/output/starLibERCC ");
            sb.append("--readFilesIn /data1/home/xiaohan/rareallele/SiPASpipeline/ST-subfq/SiPAS/");
            sb.append(i);
            sb.append("_am_R1.fq /data1/home/xiaohan/rareallele/SiPASpipeline/ST-subfq/SiPAS/");
            sb.append(i);
            sb.append("_am_R2.fq ");
            sb.append("--outFileNamePrefix /data1/home/xiaohan/rareallele/SiPASpipeline/ST-subfq/outputSiPAS/bams/");
            sb.append(i);
            sb.append("_am_ --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated  --outSAMtype SAM --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outFilterMatchNmin 80  &");
            System.out.println(sb.toString());
        }
    }

    public void GetCandidateVariants() {
        //information = nominal results;
        String information = "";
        //infilesS GT vcf
        String infileVCF = "";
        String outputDir = "";
        HashSet<String> VariantSet = new HashSet<>();
        BufferedReader br1 = IOUtils.getTextReader(information);
        String temp = null;
        String[] temps = null;
        BufferedReader br2 = IOUtils.getTextReader(infileVCF);
        String temp2 = null;
        String[] temps2 = null;
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "candidate.vcf").getAbsolutePath());
        try {
            while ((temp = br1.readLine()) != null) {
                temps = temp.split("\\s");
                if (!VariantSet.contains(temps[1])) {
                    VariantSet.add(temps[1]);
                }
            }
            String title = br2.readLine();
            bw.write(title);
            while ((temp2 = br2.readLine()) != null) {
                temps2 = temp2.split("\t");
                if (VariantSet.contains(temps2[2])) {
                    bw.write(temp2);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void SplitPhenoBychr() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/expressionboxcox3.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S3/splitexpression";
//        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionboxcox7.txt";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression";
        try {
            BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(infileS);
            BufferedWriter[] bw = new BufferedWriter[43];
            BufferedWriter[] bw1 = new BufferedWriter[43];
            for (int i = 0; i < 43; i++) {
//                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".bed").getAbsolutePath());
//                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".txt").getAbsolutePath());
                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S3expression" + i + ".bed").getAbsolutePath());
                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S3expression" + i + ".txt").getAbsolutePath());
            }
            String temp = null;
            String[] temps = null;
            temp = br.readLine();
            for (int i = 0; i < 43; i++) {
                bw[i].write(temp);
                bw[i].newLine();
                bw1[i].write(temp);
                bw1[i].newLine();
            }
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                int count = Integer.parseInt(temps[0]);
                bw[count].write(temp);
                bw[count].newLine();
                bw1[count].write(temp);
                bw1[count].newLine();
            }
            for (int i = 0; i < 43; i++) {
                bw[i].flush();
                bw[i].close();
                bw1[i].flush();
                bw1[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void CalculateeGenes() {
        String inputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6";
        String outputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/candidate";
//        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/eGene/nominals1/test/";
//        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/eGene/nominals1/candidate";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName());
            System.out.println(fs[i].getName());
        }
        nameSet.stream().forEach(f -> {
            BufferedWriter[] bw = new BufferedWriter[5];
            for (int i = 1; i < 5; i++) {
                bw[i] = IOUtils.getTextGzipWriter(new File(outputDirS, "candidateGenes" + i + ".txt.gz").getAbsolutePath());
            }
            String infile1 = new File(inputDirS, f).getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infile1);
            String temp = null;
            String[] temps = null;
            int count = 0;
            try {
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\\s");
//                System.out.println(temps[10]);
                    Float p = Float.parseFloat(temps[4]);
                    if (p <= 0.00001) {
                        bw[1].write(temp);
                        bw[1].newLine();
                        count++;
                    }
                    if (p <= 0.0001) {
                        bw[2].write(temp);
                        bw[2].newLine();
                        count++;
                    }
                    if (p <= 0.001) {
                        bw[3].write(temp);
                        bw[3].newLine();
                        count++;
                    }
                    if (p <= 0.01) {
                        bw[4].write(temp);
                        bw[4].newLine();
                        count++;
                    }
                }
                for (int i = 1; i < 5; i++) {
                    bw[i].flush();
                    bw[i].close();
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void ifanyindex() {
        String inputfileS = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/input/6/L006_R1_001.fastq.gz";
        BufferedReader br = IOUtils.getTextGzipReader(inputfileS);
        String temp = null;
        String index = "ATTATA";
        String currentindex = null;
        try {
            while ((temp = br.readLine()) != null) {
                currentindex = temp.split(":")[9].substring(1, 6);
                if (currentindex.equals(index)) {
                    System.out.println("TRUE");
                    br.readLine();
                    br.readLine();
                    br.readLine();
                }
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void read10lines() {
        String infileS = "L006_R2_001";
        String inputDir = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/input/6/";
        String outputDir = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/input/6";
        BufferedReader br = IOUtils.getTextGzipReader(inputDir + infileS + ".fastq.gz");
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, infileS + "-10lines.fq").getAbsolutePath());
        String temp = null;
        int count = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if (count < 10) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                count++;
                if (count > 10) {
                    break;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void parseFqByIndexAndBarcode() {
        HashMap barcodeStrain = new HashMap();
        HashMap indexStrain = new HashMap();
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rnaseq/RNA-seq-coleoptile-20180315-1.txt");
        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/Book1.txt");

        List<String> barcodeList = new ArrayList<String>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            barcodeStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0));
            barcodeList.add(t.getCell(i, 2).substring(1));
        }
        String inputDirS = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/input/6/";
//        String inputDirS = "/data1/home/xiaohan/coleoptile/unsplitedtest";
        String outputDirS = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/output/";
//        String outputDirS = "/data1/home/xiaohan/coleoptile/splitedtest";
        File[] fs = new File(inputDirS).listFiles();
        fs = pgl.infra.utils.IOUtils.listFilesEndsWith(fs, ".fastq.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                String infile1 = new File(inputDirS, p + "_R1_001.fastq.gz").getAbsolutePath();
                String infile2 = new File(inputDirS, p + "_R2_001.fastq.gz").getAbsolutePath();
                BufferedWriter[] bw = new BufferedWriter[barcodeList.size()];
                BufferedWriter[] bw1 = new BufferedWriter[barcodeList.size()];
                for (int i = 0; i < barcodeList.size(); i++) {
                    bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i)) + "_R1.fq").getAbsolutePath());
                    bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDirS, barcodeStrain.get(barcodeList.get(i)) + "_R2.fq").getAbsolutePath());
                }
                BufferedReader br = pgl.infra.utils.IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = pgl.infra.utils.IOUtils.getTextGzipReader(infile2);
                int pos = -1;
                String temp = null;
                String seq = null;
                String index = null;
                String currentBarcode = null;
                while ((temp = br.readLine()) != null) {
                    index = temp.split(":")[9].substring(1, 6);
                    if (indexStrain.get(index) != null) {
                        seq = br.readLine();
                        currentBarcode = seq.substring(1, 8);
                        if (barcodeStrain.get(currentBarcode) != null) {
                            pos = barcodeList.indexOf(currentBarcode);
                            bw[pos].write(temp);
                            bw[pos].newLine();
                            bw[pos].write(seq);
                            bw[pos].newLine();
                            bw[pos].write(br.readLine());
                            bw[pos].newLine();
                            bw[pos].write(br.readLine());
                            bw[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                        } else {
                            br.readLine();
                            br.readLine();
                            br1.readLine();
                            br1.readLine();
                            br1.readLine();
                            br1.readLine();
                        }
                    } else {
                        br.readLine();
                        br.readLine();
                        br.readLine();
                        br1.readLine();
                        br1.readLine();
                        br1.readLine();
                        br1.readLine();
                    }

                }
                for (int i = 0; i < barcodeList.size(); i++) {
                    bw[i].flush();
                    bw[i].close();
                    bw1[i].flush();
                    bw1[i].close();
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    private void parseFqByIndex() {
        HashMap barcodeStrain = new HashMap();
        HashMap indexStrain = new HashMap();
        long startTimePoint = System.nanoTime();
//        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rnaseq/RNA-seq-coleoptile-20180315-1.txt");
        RowTable<String> t = new RowTable<>("/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest/Book1.txt");

        List<String> IndexList = new ArrayList<String>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            indexStrain.put(t.getCell(i, 1).substring(1), t.getCell(i, 0));
            //barcodeStrain.put(t.getCell(i, 2).substring(1), t.getCell(i, 0));
            IndexList.add(t.getCell(i, 1).substring(1));
        }
        String inputDirS = "/data2/junxu/SiPASData/three-lanes-firsttime-copydata0129/P101SC18112845-01-F004-WSW50-0129-weifen";
//        String inputDirS = "/data1/home/xiaohan/coleoptile/unsplitedtest";
        String outputDirS = "/data1/home/xiaohan/rareallele/SiPASpipeline/ampmtest";
//        String outputDirS = "/data1/home/xiaohan/coleoptile/splitedtest";
        File[] fs = new File(inputDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".fastq.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String name = fs[i].getName().split("_")[0] + "_" + fs[i].getName().split("_")[1] + "_" + fs[i].getName().split("_")[2];
            nameSet.add(name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                String infile1 = new File(inputDirS, p + "_R1_001.fastq.gz").getAbsolutePath();
                String infile2 = new File(inputDirS, p + "_R2_001.fastq.gz").getAbsolutePath();
                BufferedWriter[] bw = new BufferedWriter[IndexList.size()];
                BufferedWriter[] bw1 = new BufferedWriter[IndexList.size()];
                for (int i = 0; i < IndexList.size(); i++) {
                    bw[i] = IOUtils.getTextWriter(new File(outputDirS, IndexList.get(i) + "_R1.fastq.fq").getAbsolutePath());
                    bw1[i] = IOUtils.getTextWriter(new File(outputDirS, IndexList.get(i) + "_R2.fastq.fq").getAbsolutePath());
                }
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                int pos = -1;
                String temp = null;
                String seq = null;
                String index = null;
                String currentindex = null;
                //String currentBarcode = null;
                while ((temp = br.readLine()) != null) {
                    index = temp.split(":")[9].substring(1, 6);
                    if (indexStrain.get(index) != null) {
                        currentindex = index;
                        if (barcodeStrain.get(currentindex) != null) {
                            pos = IndexList.indexOf(currentindex);
                            bw[pos].write(temp);
                            bw[pos].newLine();
                            bw[pos].write(seq);
                            bw[pos].newLine();
                            bw[pos].write(br.readLine());
                            bw[pos].newLine();
                            bw[pos].write(br.readLine());
                            bw[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                            bw1[pos].write(br1.readLine());
                            bw1[pos].newLine();
                        } else {
                            br.readLine();
                            br.readLine();
                            br1.readLine();
                            br1.readLine();
                            br1.readLine();
                            br1.readLine();
                        }
                    } else {
                        br.readLine();
                        br.readLine();
                        br.readLine();
                        br1.readLine();
                        br1.readLine();
                        br1.readLine();
                        br1.readLine();
                    }

                }
                for (int i = 0; i < IndexList.size(); i++) {
                    bw[i].flush();
                    bw[i].close();
                    bw1[i].flush();
                    bw1[i].close();
                }
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    public void rankcorrelation() {
        String exprfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/exprfile.txt";
        String snpfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/SNPcountTable.txt";
        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        BufferedReader brexpr = IOUtils.getTextReader(exprfile);
        BufferedReader brsnp = IOUtils.getTextReader(snpfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, "rank.txt").getAbsolutePath());
        String tempexpr = null;
        String tempsnp = null;
        String[] tempe = null;
        String[] temps = null;
        String Sample = "B18-E002\tB18-E007\tB18-E008\tB18-E010\tB18-E011\tB18-E014\tB18-E016\tB18-E018\tB18-E023\tB18-E024\tB18-E029\tB18-E032\tB18-E035\tB18-E038\tB18-E043\tB18-E045\tB18-E046\tB18-E049\tB18-E051\tB18-E062\tB18-E065\tB18-E070\tB18-E072\tB18-E074\tB18-E081\tB18-E082\tB18-E087\tB18-E089\tB18-E097\tB18-E099\tB18-E115\tB18-E118\tB18-E124\tB18-E127\tB18-E134\tB18-E138\tB18-E139\tB18-E141\tB18-E152\tB18-E166\tB18-E170\tB18-E180\tB18-E184\tB18-E185\tB18-E188\tB18-E199\tB18-E203\tB18-E204\tB18-E205\tB18-E210\tB18-E214\tB18-E215\tB18-E218\tB18-E219\tB18-E228\tB18-E233\tB18-E236\tB18-E237\tB18-E242\tB18-E244\tB18-E245\tB18-E251\tB18-E252\tB18-E253\tB18-E256\tB18-E262\tB18-E265\tB18-E267\tB18-E270\tB18-E271\tB18-E273\tB18-E277\tB18-E280\tB18-E286\tB18-E288\tB18-E289\tB18-E290\tB18-E298\tB18-E299\tB18-E305\tB18-E306\tB18-E312\tB18-E316\tB18-E318\tB18-E320\tB18-E324\tB18-E330\tB18-E332\tB18-E335\tB18-E337\tB18-E346\tB18-E347\tB18-E348\tB18-E355\tB18-E356\tB18-E357";
        String[] SampleName = Sample.split("\t");
        HashMap<String, Integer> exprRank = new HashMap<>();
        try {
            while ((tempsnp = brsnp.readLine()) != null) {
                if (tempsnp.startsWith("#")) {
                    continue;
                }
                temps = tempsnp.split("\t");
                bw.write(temps[0] + "\t");
                for (int i = 1; i < temps.length; i++) {
                    exprRank.put(SampleName[i - 1], Integer.parseInt(temps[i]));
                }
                tempexpr = brexpr.readLine();
                tempe = tempexpr.split("\t");
                for (int j = 1; j < tempe.length; j++) {
                    bw.write(exprRank.get(tempe[j]) + "\t");
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            brexpr.close();
            brsnp.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void SNPTable() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/SNPcount.txt";
        String Table = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/Top_512__expressed_genes_median_counts.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        Set<String> SNP = new HashSet<>();
        String geneName = null;
        String temp = null;
        String[] SNPtemp = null;
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "SNPcountTable.txt").getAbsolutePath());
        BufferedReader br = IOUtils.getTextReader(infileS);
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp.toString());
                    bw.newLine();
                    continue;
                }
                SNP.add(temp);
            }
            SNPtemp = SNP.toArray(new String[SNP.size()]);
            RowTable<String> t = new RowTable<>(Table);
            for (int i = 0; i < t.getRowNumber(); i++) {
                geneName = t.getCell(i, 0);
                System.out.println(geneName);
                for (int j = 0; j < SNPtemp.length; j++) {
                    if (SNPtemp[j].startsWith(geneName)) {
                        bw.write(SNPtemp[j].toString());
                        bw.newLine();
                        continue;
                    }
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void SNPcount() {
        String VCFfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/5kSNP.vcf";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/";
        String expressionfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/0.05/phenotypes36.txt";
        BufferedReader brVCF = IOUtils.getTextReader(VCFfile);
        BufferedReader brexpr = IOUtils.getTextReader(expressionfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "SNPcount.txt").getAbsolutePath());
        String tempVCF = null;
        String[] VCF = null;
        String Sample = "B18-E002\tB18-E007\tB18-E008\tB18-E010\tB18-E011\tB18-E014\tB18-E016\tB18-E018\tB18-E023\tB18-E024\tB18-E029\tB18-E032\tB18-E035\tB18-E038\tB18-E043\tB18-E045\tB18-E046\tB18-E049\tB18-E051\tB18-E062\tB18-E065\tB18-E070\tB18-E072\tB18-E074\tB18-E081\tB18-E082\tB18-E087\tB18-E089\tB18-E097\tB18-E099\tB18-E115\tB18-E118\tB18-E124\tB18-E127\tB18-E134\tB18-E138\tB18-E139\tB18-E141\tB18-E152\tB18-E166\tB18-E170\tB18-E180\tB18-E184\tB18-E185\tB18-E188\tB18-E199\tB18-E203\tB18-E204\tB18-E205\tB18-E210\tB18-E214\tB18-E215\tB18-E218\tB18-E219\tB18-E228\tB18-E233\tB18-E236\tB18-E237\tB18-E242\tB18-E244\tB18-E245\tB18-E251\tB18-E252\tB18-E253\tB18-E256\tB18-E262\tB18-E265\tB18-E267\tB18-E270\tB18-E271\tB18-E273\tB18-E277\tB18-E280\tB18-E286\tB18-E288\tB18-E289\tB18-E290\tB18-E298\tB18-E299\tB18-E305\tB18-E306\tB18-E312\tB18-E316\tB18-E318\tB18-E320\tB18-E324\tB18-E330\tB18-E332\tB18-E335\tB18-E337\tB18-E346\tB18-E347\tB18-E348\tB18-E355\tB18-E356\tB18-E357";
        String[] SampleName = null;
        SampleName = Sample.split("\t");
        String tempexpr = null;
        String[] expr = null;
        HashMap<String, String> TSSMap = new HashMap<>();
        Set<String> TSSset = new HashSet<>();
        String[] TSS = null;
        try {
            while ((tempexpr = brexpr.readLine()) != null) {
                if (tempexpr.startsWith("#")) continue;
                expr = tempexpr.split("\t");
//                System.out.println(expr[1]);
                TSSset.add(expr[1]);
                TSSMap.put(expr[1], expr[3]);
            }
            TSS = TSSset.toArray(new String[TSSset.size()]);
            Arrays.sort(TSS);
            int[][] SNPcount = new int[TSS.length][SampleName.length];

            while ((tempVCF = brVCF.readLine()) != null) {
                if (tempVCF.startsWith("##")) continue;
                if (tempVCF.startsWith("#")) {
                    bw.write("SNPcount" + "\t");
                    bw.write(Sample);
                    bw.newLine();
                    continue;
                }
                VCF = tempVCF.split("\t");
                for (int i = 3; i < SampleName.length + 3; i++) {
                    for (int j = 0; j < TSS.length; j++) {
                        int StartPoint = Integer.parseInt(TSS[j]);
                        int distance = StartPoint - Integer.parseInt(VCF[1]);
//                            if(distance < 5120 && distance > 0 && VCF[i].endsWith("1")){
//                                SNPcount[j][i - 3] = SNPcount[j][i-3] + 1 ;
//                            }
                        if (distance < 5120 && distance > 0 && VCF[i].equals("0/1")) {
                            SNPcount[j][i - 3] = SNPcount[j][i - 3] + 1;
                        }
                        if (distance < 5120 && distance > 0 && VCF[i].equals("1/1")) {
                            SNPcount[j][i - 3] = SNPcount[j][i - 3] + 2;
                        }
                    }
                }
            }
            for (int m = 0; m < TSS.length; m++) {
                bw.write(TSSMap.get(TSS[m]) + "\t");
                for (int n = 0; n < SampleName.length; n++) {
                    bw.write(SNPcount[m][n] + "\t");
                }
                bw.newLine();
            }
            bw.flush();
            bw.close();
            brVCF.close();
            brexpr.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void expressionRank() {
        String infile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression/S7expression2.bed";
        String[] expression = null;
        String[] name = null;
        double[] expressionlevel = new double[96];
        BufferedReader br = IOUtils.getTextReader(infile);
        String expr = null;
        HashMap<Double, String> Namewithexpression = new HashMap<>();
        try {
            while ((expr = br.readLine()) != null) {
                if (expr.startsWith("#")) {
                    name = expr.split("\t");
                    continue;
                }
                expression = expr.split("\t");
                for (int i = 4; i < expressionlevel.length + 4; i++) {
                    expressionlevel[i - 4] = Double.parseDouble(expression[i]);
                    Namewithexpression.put(expressionlevel[i - 4], name[i]);
                    System.out.println(name[i]);
                }
                Arrays.sort(expressionlevel);
                for (int i = 1; i < expressionlevel.length + 1; i++) {
                    System.out.println(Namewithexpression.get(expression[i]));
                }
                continue;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void get5kSNPcount() {
        //DS
//        String VCFfileDir = "/data2/xiaohan/addinfo/S7";
        //GT
//        String VCFfileDir = "/data2/xiaohan/GT/S3";
//        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S3";
        String VCFfileDir = "/data2/xiaohan/GT/S7";
        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S7";
//        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S3";
        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S7";
        File[] fs = new File(VCFfileDir).listFiles();
        fs = pgl.infra.utils.IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String name = fs[i].getName().split("\\.")[0].split("r")[1].split("G")[0];
            nameSet.add(name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                System.out.println("this is running :" + p);
                String[] TSS = null;
                BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "chr" + p + "GT.vcf.gz").getAbsolutePath());
//                BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "genotypes" + p + ".vcf.gz").getAbsolutePath());
//                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S3expression" + p + ".bed.gz").getAbsolutePath());
                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S7expression" + p + ".bed.gz").getAbsolutePath());
                BufferedWriter[] bwS = new BufferedWriter[12];
                bwS[0] = IOUtils.getTextWriter(new File(outputDir, p + "_5k_all_count.txt").getAbsolutePath());
                bwS[11] = IOUtils.getTextWriter(new File(outputDir, p + "_10k_all_count.txt").getAbsolutePath());
                for (int i = 1; i < 11; i++) {
                    int count = i - 1;
                    bwS[i] = IOUtils.getTextWriter(new File(outputDir, p + "_" + count + "k_" + i + "k_count.txt").getAbsolutePath());
                }
                String tempVCF = null;
                String[] VCF = null;
                String tempexpr = null;
                String[] expr = null;
                HashMap<String, String> TSSMap = new HashMap<>();
                Set<String> TSSset = new HashSet<>();
//                S3
//                String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//                S7
                String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
                String[] names = name.split(",");
                System.out.println("this is reading expression file:" + p);
                while ((tempexpr = brexpr.readLine()) != null) {
                    if (tempexpr.startsWith("#")) continue;
                    expr = tempexpr.split("\t");
                    TSSMap.put(expr[3], expr[1]);
                    TSSset.add(expr[3]);
                }
                TSS = TSSset.toArray(new String[TSSset.size()]);
                Arrays.sort(TSS);
                int[][] count5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5k[i][j] = 0;
                    }
                }
                int[][] count10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count10k[i][j] = 0;
                    }
                }
                int[][] count0_1k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count0_1k[i][j] = 0;
                    }
                }
                int[][] count1_2k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count1_2k[i][j] = 0;
                    }
                }
                int[][] count2_3k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count2_3k[i][j] = 0;
                    }
                }
                int[][] count3_4k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count3_4k[i][j] = 0;
                    }
                }
                int[][] count4_5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count4_5k[i][j] = 0;
                    }
                }
                int[][] count5_6k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5_6k[i][j] = 0;
                    }
                }
                int[][] count6_7k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count6_7k[i][j] = 0;
                    }
                }
                int[][] count7_8k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count7_8k[i][j] = 0;
                    }
                }
                int[][] count8_9k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count8_9k[i][j] = 0;
                    }
                }
                int[][] count9_10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count9_10k[i][j] = 0;
                    }
                }
                System.out.println("complete create int array: " + p);
                while ((tempVCF = brVCF.readLine()) != null) {
                    if (tempVCF.startsWith("#")) {
                        VCF = tempVCF.split("\t");
                        for (int i = 0; i < 12; i++) {
                            bwS[i].write(VCF[0] + "\t");
                            for (int m = 9; m < names.length+9; m++) {
                                bwS[i].write(VCF[m] + "\t");
                            }
                            bwS[i].newLine();
                        }
                        continue;
                    }
                    VCF = tempVCF.split("\t");
                    int snpsite = parseInt(VCF[1]);
                    for (int i = 0; i < TSS.length; i++) {
                        int startsite = Integer.parseInt(TSSMap.get(TSS[i]));
                        int distance = (int) (startsite - snpsite);
                        if (distance > 0 && distance < 5000) {
                            for (int m = 9; m < VCF.length; m++) {
//                                if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count5k[i][m - 9]++;
                                }
//                                if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count5k[i][m - 9]++;
                                    count5k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 0 && distance < 10000) {
                            for (int m = 9; m < VCF.length; m++) {
//                                if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count10k[i][m - 9]++;
                                }
//                                if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count10k[i][m - 9]++;
                                    count10k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 0 && distance < 1000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count0_1k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count0_1k[i][m - 9]++;
                                    count0_1k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 1000 && distance < 2000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count1_2k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count1_2k[i][m - 9]++;
                                    count1_2k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 2000 && distance < 3000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count2_3k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count2_3k[i][m - 9]++;
                                    count2_3k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 3000 && distance < 4000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count3_4k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count3_4k[i][m - 9]++;
                                    count3_4k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 4000 && distance < 5000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count4_5k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count4_5k[i][m - 9]++;
                                    count4_5k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 5000 && distance < 6000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count5_6k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count5_6k[i][m - 9]++;
                                    count5_6k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 6000 && distance < 7000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count6_7k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count6_7k[i][m - 9]++;
                                    count6_7k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 7000 && distance < 8000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count7_8k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count7_8k[i][m - 9]++;
                                    count7_8k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 8000 && distance < 9000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count8_9k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count8_9k[i][m - 9]++;
                                    count8_9k[i][m - 9]++;
                                }
                            }
                        }
                        if (distance > 9000 && distance < 10000) {
                            for (int m = 9; m < VCF.length; m++) {
                                //if (VCF[m].equals("1.0") || VCF[m].equals("-1.0")) {
                                if (VCF[m].equals("0/1") || VCF[m].equals("0/2") || VCF[m].equals("0/3")) {
                                    count9_10k[i][m - 9]++;
                                }
                                //if (VCF[m].equals("2.0")) {
                                if (!VCF[m].equals("0/1") && !VCF[m].equals("0/2") && !VCF[m].equals("0/3") && !VCF[m].equals("0/0")) {
                                    count9_10k[i][m - 9]++;
                                    count9_10k[i][m - 9]++;
                                }
                            }
                        }
                        continue;
                    }
                }
                System.out.println("complete making int array : " + p);
                for (int i = 0; i < TSS.length; i++) {
                    bwS[0].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[0].write(count5k[i][j] + "\t");
                    }
                    bwS[0].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[11].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[11].write(count10k[i][j] + "\t");
                    }
                    bwS[11].newLine();
                }
                System.out.println("complete writing file : " + p + "_10k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[1].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[1].write(count0_1k[i][j] + "\t");
                    }
                    bwS[1].newLine();

                }
                System.out.println("complete writing file : " + p + "_0k_1k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[2].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[2].write(count1_2k[i][j] + "\t");
                    }
                    bwS[2].newLine();
                }
                System.out.println("complete writing file : " + p + "_1k_2k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[3].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[3].write(count2_3k[i][j] + "\t");
                    }
                    bwS[3].newLine();
                }
                System.out.println("complete writing file : " + p + "_2k_3k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[4].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[4].write(count3_4k[i][j] + "\t");
                    }
                    bwS[4].newLine();
                }
                System.out.println("complete writing file : " + p + "_3k_4k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[5].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[5].write(count4_5k[i][j] + "\t");
                    }
                    bwS[5].newLine();
                }
                System.out.println("complete writing file : " + p + "_4k_5k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[6].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[6].write(count5_6k[i][j] + "\t");
                    }
                    bwS[6].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_6k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[7].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[7].write(count6_7k[i][j] + "\t");
                    }
                    bwS[7].newLine();
                }
                System.out.println("complete writing file : " + p + "_6k_7k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[8].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[8].write(count7_8k[i][j] + "\t");
                    }
                    bwS[8].newLine();
                }
                System.out.println("complete writing file : " + p + "_7k_8k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[9].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[9].write(count8_9k[i][j] + "\t");
                    }
                    bwS[9].newLine();
                }
                System.out.println("complete writing file : " + p + "_8k_9k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[10].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[10].write(count9_10k[i][j] + "\t");
                    }
                    bwS[10].newLine();
                }
                System.out.println("complete writing file : " + p + "_9k_10k_count.txt");
                System.out.println("complete writing file : " + p);
                for (int i = 0; i < 12; i++) {
                    bwS[i].flush();
                    bwS[i].close();
                }
                brexpr.close();
                brVCF.close();
            } catch (
                    Exception e) {
                e.printStackTrace();
            }
        });
    }

    
    public void getupstreamSNPcount() {
        //GT
        String VCFfileDir = "/data2/xiaohan/GT/S3";
        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S3";
//        String VCFfileDir = "/data2/xiaohan/GT/S7";
//        String expressionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S7";
        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S3/minor";
//        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/S7/minor";
        File[] fs = new File(VCFfileDir).listFiles();
        fs = pgl.infra.utils.IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        HashSet<String> nameSet = new HashSet<String>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String name = fs[i].getName().split("\\.")[0].split("r")[1].split("G")[0];
            nameSet.add(name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                System.out.println("this is running :" + p);
                String[] TSS = null;
                BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "chr" + p + "GT.vcf.gz").getAbsolutePath());
                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S3expression" + p + ".bed.gz").getAbsolutePath());
//                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(expressionfileDir, "S7expression" + p + ".bed.gz").getAbsolutePath());
                BufferedWriter[] bwS = new BufferedWriter[12];
                bwS[0] = IOUtils.getTextWriter(new File(outputDir, p + "_5k_all_count.txt").getAbsolutePath());
                bwS[11] = IOUtils.getTextWriter(new File(outputDir, p + "_10k_all_count.txt").getAbsolutePath());
                for (int i = 1; i < 11; i++) {
                    int count = i - 1;
                    bwS[i] = IOUtils.getTextWriter(new File(outputDir, p + "_" + count + "k_" + i + "k_count.txt").getAbsolutePath());
                }
                String tempVCF = null;
                String[] VCF = null;
                String tempexpr = null;
                String[] expr = null;
                HashMap<String, String> TSSMap = new HashMap<>();
                Set<String> TSSset = new HashSet<>();
//                S3
                String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//                S7
//                String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
                String[] names = name.split(",");
                System.out.println("this is reading expression file:" + p);
                while ((tempexpr = brexpr.readLine()) != null) {
                    if (tempexpr.startsWith("#")) continue;
                    expr = tempexpr.split("\t");
                    TSSMap.put(expr[3], expr[1]);
                    TSSset.add(expr[3]);
                }
                TSS = TSSset.toArray(new String[TSSset.size()]);
                Arrays.sort(TSS);
                int[][] count5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5k[i][j] = 0;
                    }
                }
                int[][] count10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count10k[i][j] = 0;
                    }
                }
                int[][] count0_1k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count0_1k[i][j] = 0;
                    }
                }
                int[][] count1_2k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count1_2k[i][j] = 0;
                    }
                }
                int[][] count2_3k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count2_3k[i][j] = 0;
                    }
                }
                int[][] count3_4k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count3_4k[i][j] = 0;
                    }
                }
                int[][] count4_5k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count4_5k[i][j] = 0;
                    }
                }
                int[][] count5_6k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count5_6k[i][j] = 0;
                    }
                }
                int[][] count6_7k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count6_7k[i][j] = 0;
                    }
                }
                int[][] count7_8k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count7_8k[i][j] = 0;
                    }
                }
                int[][] count8_9k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count8_9k[i][j] = 0;
                    }
                }
                int[][] count9_10k = new int[TSS.length][names.length];
                for (int i = 0; i < TSS.length; i++) {
                    for (int j = 0; j < names.length; j++) {
                        count9_10k[i][j] = 0;
                    }
                }
                System.out.println("complete create int array: " + p);
                while ((tempVCF = brVCF.readLine()) != null) {
                    if (tempVCF.startsWith("#")) {
                        VCF = tempVCF.split("\t");
                        for (int i = 0; i < 12; i++) {
                            bwS[i].write(VCF[0] + "\t");
                            for (int m = 9; m < names.length+9; m++) {
                                bwS[i].write(VCF[m] + "\t");
                            }
                            bwS[i].newLine();
                        }
                        continue;
                    }
                    VCF = tempVCF.split("\t");
                    int snpsite = parseInt(VCF[1]);
                    for (int i = 0; i < TSS.length; i++) {
                        int startsite = Integer.parseInt(TSSMap.get(TSS[i]));
                        int distance = (int) (startsite - snpsite);
                        if (distance >= 0 && distance < 5000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count5k[i][m - 9]++;
                                        count5k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count5k[i][m - 9]++;
                                        count5k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 0 && distance < 10000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count10k[i][m - 9]++;
                                        count10k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count10k[i][m - 9]++;
                                        count10k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 0 && distance < 1000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count0_1k[i][m - 9]++;
                                        count0_1k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count0_1k[i][m - 9]++;
                                        count0_1k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 1000 && distance < 2000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count1_2k[i][m - 9]++;
                                        count1_2k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count1_2k[i][m - 9]++;
                                        count1_2k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 2000 && distance < 3000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count2_3k[i][m - 9]++;
                                        count2_3k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count2_3k[i][m - 9]++;
                                        count2_3k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 3000 && distance < 4000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count3_4k[i][m - 9]++;
                                        count3_4k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count3_4k[i][m - 9]++;
                                        count3_4k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 4000 && distance < 5000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count4_5k[i][m - 9]++;
                                        count4_5k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count4_5k[i][m - 9]++;
                                        count4_5k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 5000 && distance < 6000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count5_6k[i][m - 9]++;
                                        count5_6k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count5_6k[i][m - 9]++;
                                        count5_6k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 6000 && distance < 7000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count6_7k[i][m - 9]++;
                                        count6_7k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count6_7k[i][m - 9]++;
                                        count6_7k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 7000 && distance < 8000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count7_8k[i][m - 9]++;
                                        count7_8k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count7_8k[i][m - 9]++;
                                        count7_8k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 8000 && distance < 9000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count8_9k[i][m - 9]++;
                                        count8_9k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count8_9k[i][m - 9]++;
                                        count8_9k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        if (distance >= 9000 && distance < 10000) {
                            int number0 = 0;
                            int number1 = 0;
                            for(int m = 9;m < VCF.length; m++){
                                if(VCF[m].equals("0/0")){
                                    number0 ++;
                                }
                                if(VCF[m].equals("1/1")){
                                    number1 ++;
                                }
                            }
                            if(number0 > number1){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("1/1")) {
                                        count9_10k[i][m - 9]++;
                                        count9_10k[i][m - 9]++;
                                    }
                                }
                            }
                            else if(number1 > number0){
                                for (int m = 9; m < VCF.length; m++) {
                                    if (VCF[m].equals("0/0")) {
                                        count9_10k[i][m - 9]++;
                                        count9_10k[i][m - 9]++;
                                    }
                                }
                            }
                        }
                        continue;
                    }
                }
                System.out.println("complete making int array : " + p);
                for (int i = 0; i < TSS.length; i++) {
                    bwS[0].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[0].write(count5k[i][j] + "\t");
                    }
                    bwS[0].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[11].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[11].write(count10k[i][j] + "\t");
                    }
                    bwS[11].newLine();
                }
                System.out.println("complete writing file : " + p + "_10k_all_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[1].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[1].write(count0_1k[i][j] + "\t");
                    }
                    bwS[1].newLine();

                }
                System.out.println("complete writing file : " + p + "_0k_1k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[2].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[2].write(count1_2k[i][j] + "\t");
                    }
                    bwS[2].newLine();
                }
                System.out.println("complete writing file : " + p + "_1k_2k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[3].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[3].write(count2_3k[i][j] + "\t");
                    }
                    bwS[3].newLine();
                }
                System.out.println("complete writing file : " + p + "_2k_3k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[4].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[4].write(count3_4k[i][j] + "\t");
                    }
                    bwS[4].newLine();
                }
                System.out.println("complete writing file : " + p + "_3k_4k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[5].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[5].write(count4_5k[i][j] + "\t");
                    }
                    bwS[5].newLine();
                }
                System.out.println("complete writing file : " + p + "_4k_5k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[6].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[6].write(count5_6k[i][j] + "\t");
                    }
                    bwS[6].newLine();
                }
                System.out.println("complete writing file : " + p + "_5k_6k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[7].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[7].write(count6_7k[i][j] + "\t");
                    }
                    bwS[7].newLine();
                }
                System.out.println("complete writing file : " + p + "_6k_7k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[8].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[8].write(count7_8k[i][j] + "\t");
                    }
                    bwS[8].newLine();
                }
                System.out.println("complete writing file : " + p + "_7k_8k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[9].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[9].write(count8_9k[i][j] + "\t");
                    }
                    bwS[9].newLine();
                }
                System.out.println("complete writing file : " + p + "_8k_9k_count.txt");
                for (int i = 0; i < TSS.length; i++) {
                    bwS[10].write(TSS[i] + "\t");
                    for (int j = 0; j < names.length; j++) {
                        bwS[10].write(count9_10k[i][j] + "\t");
                    }
                    bwS[10].newLine();
                }
                System.out.println("complete writing file : " + p + "_9k_10k_count.txt");
                System.out.println("complete writing file : " + p);
                for (int i = 0; i < 12; i++) {
                    bwS[i].flush();
                    bwS[i].close();
                }
                brexpr.close();
                brVCF.close();
            } catch (
                    Exception e) {
                e.printStackTrace();
            }
        });
    }
    
    public void changeName() {
        String infileS = "/data1/home/xiaohan/rareallele/fastQTL/chr36.vcf";
        String outputDir = "/data1/home/xiaohan/rareallele/fastQTL";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "Chr36.vcf").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < 1; i++) {
                        tems = temps[i].split("]");
                        bw.write("#" + tems[1] + "\t");
                    }
                    for (int i = 1; i < temps.length; i++) {
                        tems = temps[i].split("]");
                        bw.write(tems[1] + "\t");

                    }
                    bw.newLine();
                } else {
                    bw.write(temp);
                    bw.newLine();

                }

            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getGTvcf(String infileS, String outputDir) {
        BufferedReader br = IOUtils.getTextGzipReader(outputDir + infileS + ".new.vcf.gz");
        BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputDir, "genotypes" + infileS + ".vcf.gz").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        String tem = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                } else {
                    temps = temp.split("\t");
                    for (int j = 0; j < 8; j++) {
                        bw.write(temps[j] + "\t");
                    }
                    for (int i = 8; i < temps.length; i++) {
                        tems = temps[i].split(":");
                        if (tems[0].startsWith("0/2") || tems[0].startsWith("0/3")) {
                            tems[0] = "0/1";
                        }
                        if (tems[0].startsWith("1/2") || tems[0].startsWith("1/3") || tems[0].startsWith("2/2") || tems[0].startsWith("2/3") || tems[0].startsWith("3/3")) {
                            tems[0] = "1/1";
                        }
                        bw.write(tems[0] + "\t");
                    }
//                    for(int i = 8;i<9;i++){
//                        tems = temps[i].split(":");
//                        bw.write(tems[0]+"\t");
//                    }
//                    for(int m = 9;m<temps.length;m++) {
//                        bw.write(temps[m]+"\t");
//                    }
//                    for(int i = 9;i<temps.length;i++){
//                        tems = temps[i].split(":");
//                        //bw.write(tems[0]+tems[1]+"\t");
//                        bw.write(tems[0]+"\t");
//                    }
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getGTvcf2() {
        String inputDir = "/data2/xiaohan/sub7/";
        String outputDir = "/data2/xiaohan/GT/S7";
        HashSet<String> nameSet = new HashSet<String>();
        File[] fs = new File(inputDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf.gz");
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            String[] SampleName = fs[i].getName().split("\\.");
//            String Name = SampleName[0].replace(SampleName[0].substring(0, 3), "");
            String Name = fs[i].getName().split("\\.")[0].split("p")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                System.out.println("Start reading file snp" + p + "vcf.gz");
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDir, "snp" + p + ".vcf.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(new File(outputDir, "chr" + p + "GT.vcf.gz").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String[] tems = null;
//                S3
//                String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//                S7
                String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
                String[] names = name.split(",");
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##"))
                        continue;
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        for (int i = 0; i < 9; i++) {
                            bw.write(temps[i] + "\t");
                        }
                        for (int j = 0; j < names.length - 1; j++) {
                            bw.write(names[j] + "\t");
                        }
                        bw.write(names[names.length - 1]);
                        bw.newLine();
                        continue;
                    }
                    temps = temp.split("\t");
                    for (int i = 0; i < 2; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    //ID
                    bw.write("snp_" + temps[0] + "_" + temps[1] + "\t");
                    //REF ALT
                    for (int i = 3; i < 5; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    //QUAL FILTER
                    for (int i = 5; i < 7; i++) {
                        bw.write("." + "\t");
                    }
                    bw.write("INFO" + "\t" + "GT" + "\t");
                    for (int i = 9; i < temps.length; i++) {
                        tems = temps[i].split(":");
                        bw.write(tems[0] + "\t");
                    }
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void changeSampleName() {
        String inputDir = "/data2/xiaohan/sub3/";
        String outputDir = "/data2/xiaohan/GT/S3";
        int SampleCount = 2;
        BufferedReader[] br = new BufferedReader[SampleCount];
        BufferedWriter[] bw = new BufferedWriter[SampleCount];
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
//        String name = "B18-E002\tB18-E007\tB18-E008\tB18-E011\tB18-E014\tB18-E018\tB18-E023\tB18-E024\tB18-E029\tB18-E032\tB18-E035\tB18-E038\tB18-E043\tB18-E045\tB18-E046\tB18-E049\tB18-E051\tB18-E052\tB18-E062\tB18-E065\tB18-E070\tB18-E072\tB18-E074\tB18-E081\tB18-E082\tB18-E083\tB18-E087\tB18-E089\tB18-E097\tB18-E099\tB18-E115\tB18-E118\tB18-E124\tB18-E127\tB18-E134\tB18-E138\tB18-E139\tB18-E141\tB18-E152\tB18-E166\tB18-E170\tB18-E180\tB18-E184\tB18-E185\tB18-E188\tB18-E199\tB18-E203\tB18-E204\tB18-E205\tB18-E210\tB18-E214\tB18-E215\tB18-E218\tB18-E219\tB18-E227\tB18-E228\tB18-E233\tB18-E236\tB18-E237\tB18-E242\tB18-E244\tB18-E245\tB18-E251\tB18-E252\tB18-E256\tB18-E262\tB18-E265\tB18-E267\tB18-E270\tB18-E271\tB18-E273\tB18-E277\tB18-E280\tB18-E286\tB18-E288\tB18-E289\tB18-E290\tB18-E298\tB18-E299\tB18-E305\tB18-E306\tB18-E312\tB18-E316\tB18-E318\tB18-E320\tB18-E324\tB18-E330\tB18-E332\tB18-E335\tB18-E337\tB18-E346\tB18-E347\tB18-E348\tB18-E355\tB18-E356\tB18-E357"
//        String[] names = name.split("\t");
        try {
            for (int i = 2; i < SampleCount + 1; i++) {
                br[i - 1] = IOUtils.getTextReader(new File(inputDir, "snp" + i + ".vcf").getAbsolutePath());
            }
            for (int i = 2; i < SampleCount + 1; i++) {
                bw[i - 1] = IOUtils.getTextWriter(new File(outputDir, "genotypes" + i + ".vcf").getAbsolutePath());
            }
            for (int m = 2; m < SampleCount + 1; m++) {
                while ((temp = br[m - 1].readLine()) != null) {
                    if (temp.startsWith("##")) {
                        continue;
                    } else if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        for (int i = 0; i < 9; i++) {
                            bw[m - 1].write(temps[i] + "\t");
                        }
                        for (int j = 0; j < names.length - 1; j++) {
                            bw[m - 1].write(names[j] + "\t");
                        }
                        bw[m - 1].write(names[names.length - 1]);
                        bw[m - 1].newLine();
                        continue;
                    } else {
                        temps = temp.split("\t");
                        for (int i = 0; i < 3; i++) {
                            bw[m - 1].write(temps[i] + "\t");
                        }
//                        bw[m-1].write("snp_" + temps[0] + "_" + temps[1] + "\t");
                        for (int j = 3; j < 5; j++) {
                            bw[m - 1].write(temps[j] + "\t");
                        }
                        for (int i = 5; i < 8; i++) {
                            bw[m - 1].write("." + "\t");
                        }
                        for (int i = 8; i < temps.length; i++) {
//                            tems = temps[i].split(":");
//                            if(tems[0].equals("./.")){
//                                bw[m-1].write("0/0" + "\t");
//                            }
//                            else{
//                                bw[m-1].write(tems[0] + "\t");
//                            }
                            bw[m - 1].write(temps[i]);
                        }
                        bw[m - 1].newLine();
                        continue;
                    }
                }
            }
            for (int i = 1; i < SampleCount + 1; i++) {
                bw[i - 1].flush();
                bw[i - 1].close();
                br[i - 1].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void test() {
        for (int i = 0; i < 43; i++) {
            System.out.print(i + ".snp.vcf.gz" + "\t");
        }
    }

    public void countlines() {
        String infileS1 = "/data2/junxu/geneSNP/all.vcf.gz";
        String outfileDir = "/data2/xiaohan/test";
        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
        int count1 = 0;
        String temp1 = null;
        try {
            while ((temp1 = br1.readLine()) != null) {
                count1++;
                System.out.println(count1);
            }
            br1.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


//    public void countlines(){
//        String infileS1 = "/data2/junxu/geneSNP/0002.N349.snp.vcf.gz";
//        String outfileDir = "/data2/xiaohan/test";
//        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
//        BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir,"result1.txt").getAbsolutePath());
//        int count1 = 0;
//        String temp1 = null;
//        try {
//            if ((temp1 = br1.readLine()) != null) {
//                count1++;
//                bw.write(count1);
//                bw.newLine();
//            }
//            else {
//            bw.flush();bw.close();
//            br1.close();
//            }
//        }
//        catch (Exception e){
//                e.printStackTrace();
//            }
//    }

    public void checklines() {
        String infileS1 = "/data2/junxu/geneSNP/all.vcf.gz";
        String infileS2 = "/data1/home/xiaohan/rareallele/subvcf/colVCF.vcf";
        BufferedReader br1 = IOUtils.getTextGzipReader(infileS1);
        BufferedReader br2 = IOUtils.getTextReader(infileS2);
        BufferedWriter bw = IOUtils.getTextWriter(new File(infileS2, "ifTrue.txt").getAbsolutePath());
        int count1 = 0;
        int count2 = 0;
        String temp1 = null;
        String temp2 = null;
        try {
            while ((temp1 = br1.readLine()) != null) {
                count1++;
            }
            while ((temp2 = br2.readLine()) != null) {
                count2++;
            }
            if (count1 == count2) {
                System.out.println("True");
            } else {
                System.out.println("False");
            }
            bw.flush();
            bw.close();
            br1.close();
            br2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getsubVCF(String infileS, String inputDir, String outputDir) {
        String index = getVCFposition(infileS, inputDir);
//        String index ="12\n" +
//                "17\n" +
//                "18\n" +
//                "20\n" +
//                "21\n" +
//                "24\n" +
//                "26\n" +
//                "28\n" +
//                "33\n" +
//                "34\n" +
//                "38\n" +
//                "41\n" +
//                "44\n" +
//                "47\n" +
//                "52\n" +
//                "54\n" +
//                "55\n" +
//                "58\n" +
//                "59\n" +
//                "70\n" +
//                "73\n" +
//                "78\n" +
//                "80\n" +
//                "82\n" +
//                "88\n" +
//                "89\n" +
//                "94\n" +
//                "96\n" +
//                "104\n" +
//                "106\n" +
//                "119\n" +
//                "335\n" +
//                "126\n" +
//                "128\n" +
//                "134\n" +
//                "138\n" +
//                "139\n" +
//                "141\n" +
//                "149\n" +
//                "163\n" +
//                "166\n" +
//                "173\n" +
//                "177\n" +
//                "178\n" +
//                "181\n" +
//                "192\n" +
//                "195\n" +
//                "196\n" +
//                "197\n" +
//                "201\n" +
//                "205\n" +
//                "206\n" +
//                "209\n" +
//                "351\n" +
//                "328\n" +
//                "221\n" +
//                "224\n" +
//                "225\n" +
//                "229\n" +
//                "231\n" +
//                "342\n" +
//                "236\n" +
//                "237\n" +
//                "238\n" +
//                "241\n" +
//                "352\n" +
//                "249\n" +
//                "251\n" +
//                "254\n" +
//                "255\n" +
//                "257\n" +
//                "259\n" +
//                "9\n" +
//                "266\n" +
//                "329\n" +
//                "330\n" +
//                "268\n" +
//                "276\n" +
//                "277\n" +
//                "281\n" +
//                "282\n" +
//                "288\n" +
//                "292\n" +
//                "294\n" +
//                "296\n" +
//                "300\n" +
//                "357\n" +
//                "306\n" +
//                "309\n" +
//                "311\n" +
//                "317\n" +
//                "318\n" +
//                "331\n" +
//                "332\n" +
//                "322\n" +
//                "323\n";
        String[] indexes = null;
        indexes = index.split("\t");
        BufferedReader br = IOUtils.getTextGzipReader(inputDir + infileS + ".vcf");
        String temp = null;
        String[] temps = null;
        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, infileS + "-new.vcf").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    continue;
                }
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    bw.write(temps[0]);
                    for (int j = 1; j < 9; j++) {
                        bw.write("\t" + temps[j]);
                    }
                    for (int j = 0; j < names.length; j++) {
                        bw.write("\t" + names[j]);
                    }
                    bw.newLine();
                    continue;
                }
                if (!temp.startsWith("#")) {
                    temps = temp.split("\t");
                    bw.write(temps[0] + "\t" + temps[1]);
                    bw.write("\t" + "snp_" + temps[1] + "\t");
                    for (int i = 3; i < 9; i++) {
                        bw.write("\t" + temps[i]);
                    }
                    for (int i = 0; i < indexes.length; i++) {
                        bw.write("\t" + temps[parseInt(indexes[i])]);
                    }
                    bw.newLine();
                    continue;
                }
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public String getVCFposition(String infileS, String inputDir) {
        String infor = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/info.txt";
        BufferedReader info = IOUtils.getTextReader(infor);
        String sample = null;
        StringBuilder sb = new StringBuilder();
        try {
            while ((sample = info.readLine()) != null) {
                sb.append(sample + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        String SampleName = sb.toString();
//        String SampleName = "AT18488\n" +
//                "AT18493\n" +
//                "AT18494\n" +
//                "AT18496\n" +
//                "AT18497\n" +
//                "AT18500\n" +
//                "AT18502\n" +
//                "AT18504\n" +
//                "AT18509\n" +
//                "AT18510\n" +
//                "AT18514\n" +
//                "AT18517\n" +
//                "AT18520\n" +
//                "AT18523\n" +
//                "AT18528\n" +
//                "AT18530\n" +
//                "AT18531\n" +
//                "AT18534\n" +
//                "AT18535\n" +
//                "AT18546\n" +
//                "AT18549\n" +
//                "AT18554\n" +
//                "AT18556\n" +
//                "AT18558\n" +
//                "AT18564\n" +
//                "AT18565\n" +
//                "AT18570\n" +
//                "AT18572\n" +
//                "AT18580\n" +
//                "AT18582\n" +
//                "AT18597\n" +
//                "AT18984\n" +
//                "AT18606\n" +
//                "AT18608\n" +
//                "AT18615\n" +
//                "AT18619\n" +
//                "AT18620\n" +
//                "AT18622\n" +
//                "AT18632\n" +
//                "AT18646\n" +
//                "AT18650\n" +
//                "AT18659\n" +
//                "AT18663\n" +
//                "AT18664\n" +
//                "AT18667\n" +
//                "AT18678\n" +
//                "AT18681\n" +
//                "AT18682\n" +
//                "AT18683\n" +
//                "AT18688\n" +
//                "AT18692\n" +
//                "AT18693\n" +
//                "AT18696\n" +
//                "SYR-L1\n" +
//                "AT18837\n" +
//                "AT18710\n" +
//                "AT18713\n" +
//                "AT18714\n" +
//                "AT18719\n" +
//                "AT18721\n" +
//                "IRN-L2\n" +
//                "AT18727\n" +
//                "AT18728\n" +
//                "AT18729\n" +
//                "AT18732\n" +
//                "TJK-L1\n" +
//                "AT18741\n" +
//                "AT18743\n" +
//                "AT18746\n" +
//                "AT18747\n" +
//                "AT18749\n" +
//                "AT18752\n" +
//                "AFG-L1\n" +
//                "AT18761\n" +
//                "AT18838\n" +
//                "AT18839\n" +
//                "AT18765\n" +
//                "AT18773\n" +
//                "AT18774\n" +
//                "AT18779\n" +
//                "AT18780\n" +
//                "AT18786\n" +
//                "AT18790\n" +
//                "AT18792\n" +
//                "AT18794\n" +
//                "AT18798\n" +
//                "UZB-L1\n" +
//                "AT18805\n" +
//                "AT18808\n" +
//                "AT18810\n" +
//                "AT18819\n" +
//                "AT18820\n" +
//                "AT18842\n" +
//                "AT18843\n" +
//                "AT18828\n" +
//                "AT18829\n";
        String[] Sample = null;
        Sample = SampleName.split("\n");
        System.out.println(Sample[1]);
        Set<String> tempS = new HashSet<>();
        BufferedReader br = IOUtils.getTextReader(inputDir + infileS + ".vcf");
        String temp = null;
        String[] temps = null;
        String[] tempsOrigin = null;
//        StringBuilder sb = new StringBuilder();
        ArrayList<String> NameList = new ArrayList<>();
        HashMap<String, String> NameMap = new HashMap();
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    tempsOrigin = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].startsWith("AT")) {
                            NameMap.put(temps[i].substring(0, 7), temps[i]);
                            temps[i] = temps[i].substring(0, 7);
                        }
                    }
                    for (int j = 0; j < Sample.length; j++) {
                        for (int i = 0; i < temps.length; i++) {
                            if (temps[i].equals(Sample[j])) {
                                sb.append(i + "\t");
                                System.out.println(NameMap.get(temps[i]));
                                //System.out.println(i);
                                //System.out.println(tempsOrigin[i]+"\t"+Sample[j]);
                            }
                        }
                    }
                }
                if (!temp.startsWith("#")) {
                    break;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return sb.toString();
    }


    public void getTransNumber() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/S7_ZX-B18.txt";
        String infor = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/table/B18-AT.txt";
        BufferedReader br = IOUtils.getTextReader(infor);
        BufferedReader br1 = IOUtils.getTextReader(infileS);
        String temp1 = null;
        String temp = null;
        String[] temps = null;
        String[] temps1 = null;
        Set<String> Sample = new HashSet<String>();
        HashMap<String, String> transMap = new HashMap<>();
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                transMap.put(temps[0], temps[1]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        try {
            while ((temp1 = br1.readLine()) != null) {
                temps1 = temp1.split("\t");
                String SampleName = temps1[1].toString();
                if (!Sample.contains(SampleName)) {
                    Sample.add(SampleName);
                }
            }
            String[] SampleS = Sample.toArray(new String[Sample.size()]);
            Arrays.sort(SampleS);
            for (int i = 0; i < SampleS.length; i++) {
                System.out.println(SampleS[i] + "\t" + transMap.get(SampleS[i]) + "\t");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void getExampleVCF() {
        String inforS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/genotypetrans.txt";
        BufferedReader br1 = IOUtils.getTextReader(inforS);
        String temp1 = null;
        String[] temps1 = null;
        HashMap<String, String> NameMap = new HashMap<>();
        try {
            while ((temp1 = br1.readLine()) != null) {
                temps1 = temp1.split("\t");
                String Name = temps1[0];
                String transName = temps1[1];
                NameMap.put(Name, transName);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        String infileS = "/Users/yxh/Documents/RareAllele/003rawdata/header.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (NameMap.containsValue(temps[i])) {
                            System.out.println(i + "\t" + temps[i]);
                        }
                    }
                }
            }


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getup5kbSNP() throws IOException {
        //String SNPpositioninfile = "";
        String PhenotypeInfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/output/phenotypes.bed copy.txt";
        HashMap<String, Integer> GeneStartPointMap = new HashMap<String, Integer>();
        BufferedReader br = IOUtils.getTextReader(PhenotypeInfile);
        String temp = null;
        String[] temps = null;
        String geneName = null;
        String geneStartPoint = null;
        Set<String> geneSet = new HashSet<String>();
        while ((temp = br.readLine()) != null) {
            if (temp.startsWith("chr")) continue;
            temps = temp.split("\t");
            geneName = temps[3];
            geneSet.add(geneName);
            geneStartPoint = temps[1];
            GeneStartPointMap.put(geneName, Integer.valueOf(geneStartPoint));
            System.out.println(GeneStartPointMap.get(geneName));


        }


    }

    public void rankGenes() {

    }

    public static void main(String[] args) {
        new rareallele();
    }
}
