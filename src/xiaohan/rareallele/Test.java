package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Test {

    private String infile;

    public Test() throws IOException {
//        String[] subinDirS = {"chr1", "chr2", "chr3", "chr4", "chr5", "chr6"};
//        String[] suboutDirS = {"candidate", "count_correlation", "effect_correlation", "others"};
//        String inputDir1 = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S3/";
//        String inputDir2 = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/";
//        for (int i = 0; i < 6; i++) {
//            String inputDirS = inputDir1 + subinDirS[i] + "/";
//            String outputDirS = inputDirS;
//            for (int j = 0; j < suboutDirS.length; j++) {
//                new File(outputDirS, suboutDirS[j]).mkdir();
//            }
//            String outputDirS1 = outputDirS + suboutDirS[0];
//            String outputDirS2 = outputDirS + suboutDirS[1];
//            String outputDirS3 = outputDirS + suboutDirS[2];
//            String outputDirS4 = outputDirS + suboutDirS[3];
//            this.CalculateeGenes(inputDirS, outputDirS1);
//            this.CalculateBycount(outputDirS1, outputDirS2);
//            this.CalculateByeffect(outputDirS1, outputDirS3);
//            this.ExtractDistanceandEffect(outputDirS1, outputDirS4);
//        }
//        for (int i = 0; i < 6; i++) {
//            String inputDirS = inputDir2 + subinDirS[i] + "/";
//            String outputDirS = inputDirS;
//            for (int j = 0; j < suboutDirS.length; j++) {
//                new File(outputDirS, suboutDirS[j]).mkdir();
//            }
//            String outputDirS1 = outputDirS + suboutDirS[0];
//            String outputDirS2 = outputDirS + suboutDirS[1];
//            String outputDirS3 = outputDirS + suboutDirS[2];
//            String outputDirS4 = outputDirS + suboutDirS[3];
//            this.CalculateeGenes(inputDirS, outputDirS1);
//            this.CalculateBycount(outputDirS1, outputDirS2);
//            this.CalculateByeffect(outputDirS1, outputDirS3);
//            this.ExtractDistanceandEffect(outputDirS1, outputDirS4);
//        }
//        this.effectboxplot();
//        this.writecode();
//        this.findTaxon();
        //this.findSNPnumber();
        //this.findFalseSample();
//        this.callposition();
//        this.findDifference();
        //this.vcffiltering();
//        this.vcfmerge();
//        this.addinfo();
//        this.charDemo();
//         this.getTPM(f);
//        this.addinfobySample();
//        this.writeNumber();
//        this.writeSet();
//        this.CalculateBycount();
//        this.CalculateByeffect();
//        this.ExtractDistanceandEffect();
//        this.homologychr1();
//        this.candidate();
//        this.writecode1();
//        this.extractsub();
//        this.extractsub_5_10();
//        this.extractexpressionsub();
        this.extractExpressionHomo();
//        String infile = "name";
//        this.write(infile);
    }

    public void write(String infile){
        this.infile = infile;
        System.out.print(infile);
    }

    public void extractexpressionsub(){
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/";
        try {
            BufferedReader br = IOUtils.getTextReader(new File(inputDir, "countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwA = IOUtils.getTextWriter(new File(inputDir, "A_countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwB = IOUtils.getTextWriter(new File(inputDir, "B_countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwD = IOUtils.getTextWriter(new File(inputDir, "D_countResult7_DESeq2.txt").getAbsolutePath());
            String temp = null;
            String[] temps = null;
            String sub = "";
            temp = br.readLine();
            bwA.write(temp);
            bwB.write(temp);
            bwD.write(temp);
            bwA.newLine();
            bwB.newLine();
            bwD.newLine();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                System.out.println(temp);
                sub = temps[3].substring(8, 9);
//                        System.out.println(sub);
                if (sub.equals("A")) {
                    bwA.write(temp);
                    bwA.newLine();
                }
                if (sub.equals("B")) {
                    bwB.write(temp);
                    bwB.newLine();
                }
                if (sub.equals("D")) {
                    bwD.write(temp);
                    bwD.newLine();
                }
                continue;
            }
            br.close();
            bwA.flush();
            bwA.close();
            bwB.flush();
            bwB.close();
            bwD.flush();
            bwD.close();
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }

    public void extractExpressionHomo(){
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable";
        String inforDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data";
        String infor = "TheABD.txt";
        try {
            BufferedReader br1 = IOUtils.getTextReader(new File(inforDir, infor).getAbsolutePath());
            BufferedReader br = IOUtils.getTextReader(new File(inputDir, "countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwA = IOUtils.getTextWriter(new File(inputDir, "Ahomo_countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwB = IOUtils.getTextWriter(new File(inputDir, "Bhomo_countResult7_DESeq2.txt").getAbsolutePath());
            BufferedWriter bwD = IOUtils.getTextWriter(new File(inputDir, "Dhomo_countResult7_DESeq2.txt").getAbsolutePath());
            String temp = null;
            String[] temps = null;
            HashSet GeneSetA = new HashSet();
            HashSet GeneSetB = new HashSet();
            HashSet GeneSetD = new HashSet();
            HashMap<String, Integer> genesubAMap = new HashMap();
            HashMap<String, Integer> genesubBMap = new HashMap();
            HashMap<String, Integer> genesubDMap = new HashMap();
            while ((temp = br1.readLine()) != null) {
                temps = temp.split("\t");
                GeneSetA.add(temps[1]);
                GeneSetB.add(temps[2]);
                GeneSetD.add(temps[3]);
                genesubAMap.put(temps[1], Integer.parseInt(temps[0]));
                genesubBMap.put(temps[2], Integer.parseInt(temps[0]));
                genesubDMap.put(temps[3], Integer.parseInt(temps[0]));
                continue;
                }
            String[] rankA = new String[17108];
            String[] rankB = new String[17108];
            String[] rankD = new String[17108];
            temp = br.readLine();
            bwA.write(temp);
            bwB.write(temp);
            bwD.write(temp);
            bwA.newLine();
            bwB.newLine();
            bwD.newLine();
            String temp1 = null;
            String[] temps1 = null;
            while ((temp1 = br.readLine()) != null) {
                temps1 = temp1.split("\t");
                if (GeneSetA.contains(temps1[3])) {
                    String gene = temps1[3];
                    int indexA = genesubAMap.get(gene);
                    rankA[indexA - 1] = temp1;
                        }
                if (GeneSetB.contains(temps1[3])) {
                    String gene = temps1[3];
                    int indexB = genesubBMap.get(gene);
                    rankB[indexB - 1] = temp1;
                }
                if (GeneSetD.contains(temps1[3])) {
                    String gene = temps1[3];
                    int indexD = genesubDMap.get(gene);
                    rankD[indexD - 1] = temp1;
                }
            }
            for (int i = 0; i < rankA.length; i++) {
                bwA.write(rankA[i]);
                bwA.newLine();
            }
            for (int i = 0; i < rankB.length; i++) {
                bwB.write(rankB[i]);
                bwB.newLine();
            }
            for (int i = 0; i < rankD.length; i++) {
                bwD.write(rankD[i]);
                bwD.newLine();
            }
            br.close();
            br1.close();
            bwA.flush();
            bwA.close();
            bwB.flush();
            bwB.close();
            bwD.flush();
            bwD.close();
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }
    
    public void extractsub_5_10() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count";
        for (int i = 1; i < 3; i++) {
            int f = i * 5;
            try {
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, "all_" + f + "k_" + "all_count.txt").getAbsolutePath());
                BufferedWriter bwA = IOUtils.getTextWriter(new File(inputDir, "A_" + f + "k_" + "all_count.txt").getAbsolutePath());
                BufferedWriter bwB = IOUtils.getTextWriter(new File(inputDir, "B_" + f + "k_" + "all_count.txt").getAbsolutePath());
                BufferedWriter bwD = IOUtils.getTextWriter(new File(inputDir, "D_" + f + "k_" + "all_count.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String sub = "";
                temp = br.readLine();
                bwA.write(temp);
                bwB.write(temp);
                bwD.write(temp);
                bwA.newLine();
                bwB.newLine();
                bwD.newLine();
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    System.out.println(temp);
                    sub = temps[0].substring(8, 9);
//                        System.out.println(sub);
                    if (sub.equals("A")) {
                        bwA.write(temp);
                        bwA.newLine();
                    }
                    if (sub.equals("B")) {
                        bwB.write(temp);
                        bwB.newLine();
                    }
                    if (sub.equals("D")) {
                        bwD.write(temp);
                        bwD.newLine();
                    }
                    continue;
                }
                br.close();
                bwA.flush();
                bwA.close();
                bwB.flush();
                bwB.close();
                bwD.flush();
                bwD.close();
            } catch (
                    Exception e) {
                e.printStackTrace();
            }
        }
    }


    public void extractsub() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/snp_count";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("_")[1].split("k")[0];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                if (Integer.parseInt(f) < 10) {
                    int count = Integer.parseInt(f) + 1;
                    BufferedReader br = IOUtils.getTextReader(new File(inputDir, "all_" + f + "k_" + count + "k_count.txt").getAbsolutePath());
                    BufferedWriter bwA = IOUtils.getTextWriter(new File(inputDir, "A_" + f + "k_" + count + "k_count.txt").getAbsolutePath());
                    BufferedWriter bwB = IOUtils.getTextWriter(new File(inputDir, "B_" + f + "k_" + count + "k_count.txt").getAbsolutePath());
                    BufferedWriter bwD = IOUtils.getTextWriter(new File(inputDir, "D_" + f + "k_" + count + "k_count.txt").getAbsolutePath());
                    String temp = null;
                    String[] temps = null;
                    String sub = "";
                    temp = br.readLine();
                    bwA.write(temp);
                    bwB.write(temp);
                    bwD.write(temp);
                    bwA.newLine();
                    bwB.newLine();
                    bwD.newLine();
                    while ((temp = br.readLine()) != null) {
                        temps = temp.split("\t");
                        System.out.println(temp);
                        sub = temps[0].substring(8, 9);
//                        System.out.println(sub);
                        if (sub.equals("A")) {
                            bwA.write(temp);
                            bwA.newLine();
                        }
                        if (sub.equals("B")) {
                            bwB.write(temp);
                            bwB.newLine();
                        }
                        if (sub.equals("D")) {
                            bwD.write(temp);
                            bwD.newLine();
                        }
                        continue;
                    }
                    br.close();
                    bwA.flush();
                    bwA.close();
                    bwB.flush();
                    bwB.close();
                    bwD.flush();
                    bwD.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void writecode1() {
        for(int i = 0; i < 43 ; i++){
            int count = i + 1;
//            System.out.print("nohup vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/");
//            System.out.print(i);
//            System.out.print(".snp.vcf.gz --maf 0 --max-maf 0.05 --out ");
//            System.out.print(i + ".snp.maf005 --recode && tabix -p " + i + ".snp.maf005.recode.vcf.gz >log1.txt 2>&1 &");
//            System.out.println("nohup bcftools view -S S7SampleName.txt /data2/xiaohan/SNP/"+i+".snp.maf005.recode.vcf.gz -Ov > /data2/xiaohan/sub7/snp"+i+".vcf &");
//            System.out.println("nohup bcftools +dosage /data2/xiaohan/sub7/snp"+i+".vcf -- -t GT >/data2/xiaohan/DS/S7/col"+i+"DS.vcf &");
//            System.out.println("nohup bcftools view -S S3SampleName.txt /data2/xiaohan/SNP/"+i+".snp.maf005.mis01.recode.vcf.gz -Ov > /data2/xiaohan/sub3/snp"+i+".vcf && bgzip /data2/xiaohan/sub3/snp"+i+".vcf &");
//            System.out.println("nohup bcftools +dosage /data2/xiaohan/sub3/snp"+i+".vcf -- -t GT >/data2/xiaohan/DS/S3/col"+i+"DS.vcf &");
//              System.out.println("snp"+count+" <- snp"+count+"[,-1]\n" +
//                                "count"+count+"=c(rep(0,92)) \n" +
//                                "for(i in 1:35345){ \n" +
//                                "  Rk=rank(exp[i,],ties.method=\"first\") \n" +
//                                "  for(j in 1:92){ \n" +
//                                "    zz=as.integer(Rk[j]) \n" +
//                                "    if(snp"+count+"[i,j]!=0){ \n" +
//                                "      count"+count+"[zz]=count"+count+"[zz]+snp"+count+"[i,j]}\n" +
//                                "    j=j+1}\n" +
//                                "  i=i+1} \n" +
//                                "perc=1:92\n" +
//                                "count"+count);
//            int number = count - 14;
//            int number1 = number + 1;
//            System.out.println("###############" + number + "_" + number1 + "###################\n" +
//                    "x <- data[1,]\n" +
//                    "x <- x[,-1]\n" +
//                    "x <- c(x)\n" +
//                    "y <- data[" + count + ",]\n" +
//                    "y <- y[,-1]\n" +
//                    "y <- c(y)");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_am"+i+"_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq5M/output/subFastqs/TRU_am"+i+"_R1.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_pm"+i+"_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq5M/output/subFastqs/TRU_pm"+i+"_R1.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_am"+i+"_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq5M/output/subFastqs/TRU_am"+i+"_R2.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_pm"+i+"_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq5M/output/subFastqs/TRU_pm"+i+"_R2.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_am"+i+"_R1.fq.gz 10000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq10M/output/subFastqs/TRU_am"+i+"_R1.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_pm"+i+"_R1.fq.gz 10000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq10M/output/subFastqs/TRU_pm"+i+"_R1.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_am"+i+"_R2.fq.gz 10000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq10M/output/subFastqs/TRU_am"+i+"_R2.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq/output/subFastqs/TRU_pm"+i+"_R2.fq.gz 10000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/Truseq10M/output/subFastqs/TRU_pm"+i+"_R2.fq");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASURam"+i+"_200415_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASURam"+i+"_200415_R1.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASURam"+i+"_200415_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASURam"+i+"_200415_R2.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASRam"+i+"_200415_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASRam"+i+"_200415_R1.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASRam"+i+"_200415_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASRam"+i+"_200415_R2.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASUam"+i+"_200415_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASUam"+i+"_200415_R1.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASUam"+i+"_200415_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASUam"+i+"_200415_R2.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASam"+i+"_200415_R1.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASam"+i+"_200415_R1.fq.gz");
//            System.out.println("seqtk sample -s100 /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS/output/subFastqs/SiPASam"+i+"_200415_R2.fq.gz 5000000 > /data1/home/xiaohan/rareallele/SiPASpipeline/20200527/SiPAS5M/output/subFastqs/SiPASam"+i+"_200415_R2.fq.gz");
//              System.out.println("rownumber <- nrow(snp"+i+")\n" +
//                        "colnumber <- ncol(snp"+i+")\n" +
//                        "colnumber <- colnumber - 1\n" +
//                        "snp"+i+" <- snp"+i+"[,-1]\n" +
//                        "count"+i+"=c(rep(0,colnumber)) \n" +
//                        "for(i in 1:rownumber){ \n" +
//                        "  Rk=rank(exp[i,],ties.method=\"first\") \n" +
//                        "  for(j in 1:colnumber){ \n" +
//                        "    zz=as.integer(Rk[j]) \n" +
//                        "    count"+i+"[zz]=count"+i+"[zz]+snp"+i+"[i,j]\n" +
//                        "    j=j+1}\n" +
//                        "  i=i+1} \n" +
//                        "perc=1:colnumber\n" +
//                        "count"+i);
              System.out.println("vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/"+i+".snp.vcf.gz --maf 0.05 --out /data2/xiaohan/commonSNP/"+i+".snp.minmaf005.mis01 --recode ");
        }
    }

    public void CalculateByeffect() {
        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr6/candidate1";
        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr6/effect_correlation1/";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(inputDirS, "candidateGenes" + f + ".txt").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextReader(new File(inputDirS, "candidateGenes" + f + ".txt").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect_count.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect_boxplot.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                int[] count = new int[300];
                for (int i = 0; i < 300; i++) {
                    count[i] = 0;
                }
//                int[] count = new int[600];
//                for (int i = 0; i < 600; i++) {
//                        count[i] = 0;
//                }
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    Double effect = Double.parseDouble(temps[5]);
                    for (int i = 0; i < 100; i++) {
                        int startpoint = i * 10000;
                        int endpoint = (i + 1) * 10000;
                        if (distance > startpoint && distance < endpoint) {
                            if (effect >= 0.5 && effect < 1 || effect >= -1 && effect < -0.5) {
                                count[i]++;
                            }
                            if (effect >= 0.1 && effect < 0.5 || effect >= 0.5 && effect < -0.1) {
                                count[i + 100]++;
                            }
                            if (effect >= 0 && effect < 0.1 || effect >= -0.1 && effect < 0) {
                                count[i + 200]++;
                            }
                        }
                    }
//                    for (int i = 0; i < 100; i++) {
//                        int startpoint = -i * 10000;
//                        int endpoint = (-i + 1) * 10000;
//                        if (distance < startpoint && distance > endpoint) {
//                            if (effect >= 0.5 && effect < 1 || effect >= -1 && effect < -0.5) {
//                                count[i+300]++;
//                            }
//                            if (effect >= 0.1 && effect < 0.5 || effect >= 0.5 && effect < -0.1) {
//                                count[i+400]++;
//                            }
//                            if (effect >= 0 && effect < 0.1 || effect >= -0.1 && effect < 0) {
//                                count[i+000]++;
//                            }
//                        }
//                    }
                }
                double[] effectsize = new double[100];
                for (int i = 0; i < 100; i++) {
                    effectsize[i] = 0;
                }
                double[] effectsizeave = new double[100];
                for (int i = 0; i < 100; i++) {
                    effectsizeave[i] = 0;
                }
                int[] number = new int[100];
//                double[] effectsize = new double[200];
//                for (int i = 0; i < 200; i++) {
//                        effectsize[i] = 0;
//                    }
//                double[] effectsizeave = new double[200];
//                for (int i = 0; i < 200; i++) {
//                        effectsizeave[i] = 0;
//                }
//                int[] number = new int[200];
                int lines1 = 0;
                while ((temp1 = br1.readLine()) != null) {
                    temps1 = temp1.split("\\s");
                    int distance = Integer.parseInt(temps1[2]);
                    Double effect = Double.parseDouble(temps1[5]);
                    lines1++;
                    for (int i = 0; i < 100; i++) {
                        int startpoint = i * 10000;
                        int endpoint = (i + 1) * 10000;
                        if (distance > startpoint && distance < endpoint) {
                            effectsize[i] = effectsize[i] + Math.abs(effect);
                            number[i]++;
                        }
                    }
//                    for (int i = 0; i < 100; i++) {
//                        int startpoint = -i * 10000;
//                        int endpoint = (-i + 1) * 10000;
//                        if (distance < startpoint && distance > endpoint) {
//                            effectsize[i+100] = effectsize[i+1000] +Math.abs(effect);
//                            number[i+100]++;
//                        }
//                    }
                    for (int i = 0; i < 100; i++) {
                        effectsizeave[i] = effectsize[i] / number[i];
                    }
//                    for (int i = 0; i < 200; i++) {
//                        effectsizeave[i] = effectsize[i]/number[i];
//                    }
                }
                //不同影响效应值的count的TSS上游的分布
                bw.write("Type" + "\t" + "Distance" + "\t" + "Count");
                bw.newLine();
                for (int i = 0; i < 100; i++) {
                    int distance1 = i * 10000;
                    bw.write("Large" + "\t" + distance1 + "\t" + count[i] + "\n"
                            + "Medium" + "\t" + distance1 + "\t" + count[i + 100] + "\n"
                            + "Small" + "\t" + distance1 + "\t" + count[i + 200]);
                    bw.newLine();
                }
//                for (int i = 0; i < 100; i++) {
//                    int distance1 = -i * 10000 ;
//                    bw.write("Large" + "\t" + distance1 + "\t" + count[i+300] + "\n"
//                            + "Medium" + "\t" + distance1 + "\t" + count[i+400] + "\n"
//                            + "Small" + "\t" + distance1 + "\t" + count[i+500]);
//                    bw.newLine();
//                }
                //不同的影响效应值的平均值的TSS上游的分布
                bw1.write("Distance" + "\t" + "Effect");
                bw1.newLine();
                for (int i = 0; i < 100; i++) {
                    int distance1 = i * 10000;
                    double effectforanalysis = effectsizeave[i];
                    bw1.write(distance1 + "\t" + effectforanalysis);
                    bw1.newLine();
                }
//                for (int i = 0; i < 100; i++) {
//                    int distance1 = -i * 10000 ;
//                    double effectforanalysis = effectsizeave[i+100];
//                    bw1.write(distance1 + "\t"+ effectforanalysis);
//                    bw1.newLine();
//                }
                //所有影响效应值在TSS上游的分布
                bw.flush();
                bw1.flush();
                bw.close();
                bw1.close();
                br.close();
                br1.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void CalculateBycount() {
        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/candidate1";
        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/count_correlation1";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            System.out.println(fs[i]);
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(inputDirS, "candidateGenes" + f + ".txt").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_count.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                int[] count = new int[10000];
                for (int i = 0; i < 10000; i++) {
                    count[i] = 0;
                }
//                int[] count = new int[200];
//                for (int i = 0; i < 200; i++) {
//                    count[i] = 0;
//                }
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    for (int i = 0; i < 10000; i++) {
                        int startpoint = i * 100;
                        int endpoint = (i + 1) * 100;
                        int startpoint1 = -i * 100;
                        int endpoint1 = (-i - 1) * 100;
                        if (distance > startpoint && distance < endpoint) {
                            count[i]++;
                        }
//                        if (distance < startpoint1 && distance > endpoint1){
//                            count[i+100]++;
//                        }
                    }
                }
                bw.write("Distance" + "\t" + "Count");
                bw.newLine();
                for (int i = 0; i < 10000; i++) {
                    int distance1 = i * 100;
                    bw.write(distance1 + "\t" + count[i]);
                    bw.newLine();
                }
//                for (int i = 0 ; i < 100; i++){
//                    int distance1 = -i * 10000;
//                    bw.write(distance1 + "\t" + count[i+100]);
//                    bw.newLine();
//                }
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void ExtractDistanceandEffect() {
        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/candidate1/";
        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/others1/";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(inputDirS, "candidateGenes" + f + ".txt").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextReader(new File(inputDirS, "candidateGenes" + f + ".txt").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDirS, f + "upstream_distance_effect.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                int lines = 0;
                int line = 0;
                while ((temp1 = br1.readLine()) != null) {
                    lines++;
                }
                System.out.println(lines);
                int[] distancevalue = new int[lines];
                double[] effectvalue = new double[lines];
                bw1.write("Distance" + "\t" + "Effect");
                bw1.newLine();
                while ((temp = br.readLine()) != null) {
                    line++;
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    double effect = Double.parseDouble(temps[5]);
                    if (distance > 0) {
                        bw1.write(distance + "\t" + effect);
                        bw1.newLine();
                    }
                    distancevalue[line - 1] = distance;
                    effectvalue[line - 1] = effect;
                }
                bw.write("Distance" + "\t" + "Effect");
                bw.newLine();
                for (int i = 0; i < lines; i++) {
                    int value1 = distancevalue[i];
                    double value2 = effectvalue[i];
                    bw.write(value1 + "\t" + value2);
                    bw.newLine();
                }
                bw.flush();
                bw1.flush();
                bw.close();
                bw1.close();
                br.close();
                br1.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void candidate() {
        String infile = "/Users/yxh/Documents/RareAllele/006information/chr1-1.txt";
        String output = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/candidate1";
        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/candidate/";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedReader brS = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
            System.out.println(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
            String temp = null;
            String temp1 = null;
            String[] temps = null;
            String[] temps1 = null;
            String geneName = null;
            BufferedWriter bw = IOUtils.getTextWriter(new File(output, "candidateGenes" + f + ".txt").getAbsolutePath());
            HashSet<String> geneNameList = new HashSet<>();
            try {
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    geneName = temps[0];
                    if (!geneNameList.contains(geneName)) {
                        geneNameList.add(geneName);
                    }
                    continue;
                }
                while ((temp1 = brS.readLine()) != null) {
                    temps1 = temp1.split("\\s");
                    if (geneNameList.contains(temps1[0])) {
                        System.out.println(temps1[0]);
                        bw.write(temp1);
                        bw.newLine();
                    }
                    continue;
                }
                br.close();
                brS.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }

        });

    }


    public void homologychr1() {
        String innfile = "/Users/yxh/Documents/RareAllele/006information/TheABD.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/006information/";
        BufferedReader br = IOUtils.getTextReader(innfile);
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "chr1homeology.txt").getAbsolutePath());
        String temp = null;
        String[] temps = null;
        String tems1 = null;
        String tems2 = null;
        String tems3 = null;

        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                tems1 = temps[0].substring(7, 8);
                tems2 = temps[1].substring(7, 8);
                tems3 = temps[2].substring(7, 8);
                if (tems1.equals("1") && tems2.equals("1") && tems3.equals("1")) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void writeSet() {
        String infor = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression/S7expression1.bed";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/set/set7";
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "setFile1").getAbsolutePath());
        BufferedReader br = IOUtils.getTextReader(infor);
        String temp = null;
        String[] temps = null;
        ArrayList<Integer> startList = new ArrayList<>();
        HashSet<Integer> Start = new HashSet<>();
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                temps = temp.split("\t");
                Start.add(Integer.parseInt(temps[1]));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        for (Integer s : Start) {
            startList.add(s);
        }
        System.out.println(Start.size());
        Collections.sort(startList);
        int setnumber = startList.get(startList.size() - 1);
        System.out.println(setnumber);
        int sets = setnumber / 10000 + 1;
        int start = 0;
        int end = 0;
        int start1 = 10000;
        try {
            for (int j = 0; j < sets + 1; j++) {
                int set = j + 1;
                bw.write("set" + set + " ");
                bw.write(j + ":" + start + "-" + start1);
                start = start1;
                start1 = start1 + 10000;
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void effectboxplot() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/chr1/others/1distance_effect.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                int distance = Integer.parseInt(temps[1]);
                for (int i = 0; i < 10; i++) {
                    int range1 = i * 100000;
                    int range2 = (i + 1) * 100000;
                    if (distance > range1 && distance < range2) {
                        System.out.println(range1 + "\t" + range2 + "\t" + temps[2]);
                    }
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void CalculateeGenes(String inputDirS, String outputDirS) {
//        String inputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6";
//        String outputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/candidate";
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


    public void ExtractDistanceandEffect(String inputDirS, String outputDirS) {
//        String inputDirS = "/Users/yxh/Documents/RareAllele/004test/eGene/nominals1/candidate/";
//        String outputDirS = "/Users/yxh/Documents/RareAllele/004test/eGene/nominals1/effect_correlation/";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDirS, f + "upstream_distance_effect.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                int lines = 0;
                int line = 0;
                while ((temp1 = br1.readLine()) != null) {
                    lines++;
                }
                System.out.println(lines);
                int[] distancevalue = new int[lines];
                double[] effectvalue = new double[lines];
                bw1.write("Distance" + "\t" + "Effect");
                bw1.newLine();
                while ((temp = br.readLine()) != null) {
                    line++;
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    double effect = Double.parseDouble(temps[5]);
                    if (distance > 0) {
                        bw1.write(distance + "\t" + effect);
                        bw1.newLine();
                    }
                    distancevalue[line - 1] = distance;
                    effectvalue[line - 1] = effect;
                }
                bw.write("Distance" + "\t" + "Effect");
                bw.newLine();
                for (int i = 0; i < lines; i++) {
                    int value1 = distancevalue[i];
                    double value2 = effectvalue[i];
                    bw.write(value1 + "\t" + value2);
                    bw.newLine();
                }
                bw.flush();
                bw1.flush();
                bw.close();
                bw1.close();
                br.close();
                br1.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void writecode() {
        for (int i = 0; i < 1000; i++) {
            int startpoint = i * 1000;
            int endpoint = (i + 1) * 1000;
            System.out.println("if(distancde > " + startpoint + " && distance < " + endpoint + "){");
            System.out.println("count[" + i + "]" + "++;");
            System.out.println("}else ");
        }
    }

    public void CalculateBycount(String inputDirS, String outputDirS) {
//        String inputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/candidate";
//        String outputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/count_correlation";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            System.out.println(fs[i]);
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_count.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                int[] count = new int[100];
                for (int i = 0; i < 100; i++) {
                    count[i] = 0;
                }
//                int[] count = new int[200];
//                for (int i = 0; i < 200; i++) {
//                    count[i] = 0;
//                }
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    for (int i = 0; i < 100; i++) {
                        int startpoint = i * 10000;
                        int endpoint = (i + 1) * 10000;
                        int startpoint1 = -i * 10000;
                        int endpoint1 = (-i - 1) * 10000;
                        if (distance > startpoint && distance < endpoint) {
                            count[i]++;
                        }
//                        if (distance < startpoint1 && distance > endpoint1){
//                            count[i+100]++;
//                        }
                    }
                }
                bw.write("Distance" + "\t" + "Count");
                bw.newLine();
                for (int i = 0; i < 100; i++) {
                    int distance1 = i * 10000;
                    bw.write(distance1 + "\t" + count[i]);
                    bw.newLine();
                }
//                for (int i = 0 ; i < 100; i++){
//                    int distance1 = -i * 10000;
//                    bw.write(distance1 + "\t" + count[i+100]);
//                    bw.newLine();
//                }
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void CalculateByeffect(String inputDirS, String outputDirS) {
//        String inputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/candidate";
//        String outputDirS = "/data1/home/xiaohan/rareallele/fastQTL/nominal/S7/chr6/effect_correlation/";
        File[] fs = new File(inputDirS).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            String Name = fs[i].getName().split("\\.")[0].split("s")[1];
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextGzipReader(new File(inputDirS, "candidateGenes" + f + ".txt.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect_count.txt").getAbsolutePath());
                BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDirS, f + "distance_effect_boxplot.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                int[] count = new int[300];
                for (int i = 0; i < 300; i++) {
                    count[i] = 0;
                }
//                int[] count = new int[600];
//                for (int i = 0; i < 600; i++) {
//                        count[i] = 0;
//                }
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\\s");
                    int distance = Integer.parseInt(temps[2]);
                    Double effect = Double.parseDouble(temps[5]);
                    for (int i = 0; i < 100; i++) {
                        int startpoint = i * 10000;
                        int endpoint = (i + 1) * 10000;
                        if (distance > startpoint && distance < endpoint) {
                            if (effect >= 0.5 && effect < 1 || effect >= -1 && effect < -0.5) {
                                count[i]++;
                            }
                            if (effect >= 0.1 && effect < 0.5 || effect >= 0.5 && effect < -0.1) {
                                count[i + 100]++;
                            }
                            if (effect >= 0 && effect < 0.1 || effect >= -0.1 && effect < 0) {
                                count[i + 200]++;
                            }
                        }
                    }
//                    for (int i = 0; i < 100; i++) {
//                        int startpoint = -i * 10000;
//                        int endpoint = (-i + 1) * 10000;
//                        if (distance < startpoint && distance > endpoint) {
//                            if (effect >= 0.5 && effect < 1 || effect >= -1 && effect < -0.5) {
//                                count[i+300]++;
//                            }
//                            if (effect >= 0.1 && effect < 0.5 || effect >= 0.5 && effect < -0.1) {
//                                count[i+400]++;
//                            }
//                            if (effect >= 0 && effect < 0.1 || effect >= -0.1 && effect < 0) {
//                                count[i+000]++;
//                            }
//                        }
//                    }
                }
                double[] effectsize = new double[100];
                for (int i = 0; i < 100; i++) {
                    effectsize[i] = 0;
                }
                double[] effectsizeave = new double[100];
                for (int i = 0; i < 100; i++) {
                    effectsizeave[i] = 0;
                }
                int[] number = new int[100];
//                double[] effectsize = new double[200];
//                for (int i = 0; i < 200; i++) {
//                        effectsize[i] = 0;
//                    }
//                double[] effectsizeave = new double[200];
//                for (int i = 0; i < 200; i++) {
//                        effectsizeave[i] = 0;
//                }
//                int[] number = new int[200];
                int lines1 = 0;
                while ((temp1 = br1.readLine()) != null) {
                    temps1 = temp1.split("\\s");
                    int distance = Integer.parseInt(temps1[2]);
                    Double effect = Double.parseDouble(temps1[5]);
                    lines1++;
                    for (int i = 0; i < 100; i++) {
                        int startpoint = i * 10000;
                        int endpoint = (i + 1) * 10000;
                        if (distance > startpoint && distance < endpoint) {
                            effectsize[i] = effectsize[i] + Math.abs(effect);
                            number[i]++;
                        }
                    }
//                    for (int i = 0; i < 100; i++) {
//                        int startpoint = -i * 10000;
//                        int endpoint = (-i + 1) * 10000;
//                        if (distance < startpoint && distance > endpoint) {
//                            effectsize[i+100] = effectsize[i+1000] +Math.abs(effect);
//                            number[i+100]++;
//                        }
//                    }
                    for (int i = 0; i < 100; i++) {
                        effectsizeave[i] = effectsize[i] / number[i];
                    }
//                    for (int i = 0; i < 200; i++) {
//                        effectsizeave[i] = effectsize[i]/number[i];
//                    }
                }
                //不同影响效应值的count的TSS上游的分布
                bw.write("Type" + "\t" + "Distance" + "\t" + "Count");
                bw.newLine();
                for (int i = 0; i < 100; i++) {
                    int distance1 = i * 10000;
                    bw.write("Large" + "\t" + distance1 + "\t" + count[i] + "\n"
                            + "Medium" + "\t" + distance1 + "\t" + count[i + 100] + "\n"
                            + "Small" + "\t" + distance1 + "\t" + count[i + 200]);
                    bw.newLine();
                }
//                for (int i = 0; i < 100; i++) {
//                    int distance1 = -i * 10000 ;
//                    bw.write("Large" + "\t" + distance1 + "\t" + count[i+300] + "\n"
//                            + "Medium" + "\t" + distance1 + "\t" + count[i+400] + "\n"
//                            + "Small" + "\t" + distance1 + "\t" + count[i+500]);
//                    bw.newLine();
//                }
                //不同的影响效应值的平均值的TSS上游的分布
                bw1.write("Distance" + "\t" + "Effect");
                bw1.newLine();
                for (int i = 0; i < 100; i++) {
                    int distance1 = i * 10000;
                    double effectforanalysis = effectsizeave[i];
                    bw1.write(distance1 + "\t" + effectforanalysis);
                    bw1.newLine();
                }
//                for (int i = 0; i < 100; i++) {
//                    int distance1 = -i * 10000 ;
//                    double effectforanalysis = effectsizeave[i+100];
//                    bw1.write(distance1 + "\t"+ effectforanalysis);
//                    bw1.newLine();
//                }
                //所有影响效应值在TSS上游的分布
                bw.flush();
                bw1.flush();
                bw.close();
                bw1.close();
                br.close();
                br1.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void writeNumber() {
        String inputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/rvtest/data/set/set7";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".bed");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            System.out.println(fs[i].getName().split("\\.")[0]);
            String Name = fs[i].getName().split("\\.")[0].split("n")[1];
            nameSet.add(Name);
//            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "setFile" + f).getAbsolutePath());
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, "S7expression" + f + ".bed").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                ArrayList<Integer> startList = new ArrayList<>();
                HashMap<String, Integer> StartMap = new HashMap<>();
                HashSet<String> geneName = new HashSet<>();
                HashSet<Integer> Start = new HashSet<>();
                HashSet<Integer> SetStart = new HashSet<>();
                ArrayList<Integer> SetList = new ArrayList<>();
                ArrayList<String> geneNameList = new ArrayList<>();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    Start.add(Integer.parseInt(temps[1]));
                    geneName.add(temps[3]);
                    StartMap.put(temps[3], Integer.parseInt(temps[1]));
                }
                for (String s : geneName) {
                    geneNameList.add(s);
                }
                Collections.sort(geneNameList);
                int start = 0;
                int end = 0;
                for (int j = 0; j < geneNameList.size(); j++) {
                    String set = geneNameList.get(j);
                    System.out.println(set);
                    System.out.println(StartMap.get(set));
                    start = StartMap.get(set) - 1000000;
                    if (start < 0) {
                        start = 0;
                    }
                    end = StartMap.get(set);
                    int count = 0;
                    while (start < end) {
                        int latter = start + 10000;
                        count++;
                        bw.write(set + "_" + count + " " + f + ":" + start + "-" + latter);
                        bw.newLine();
                        start = start + 10000;
                    }
                }
                for (Integer s : SetStart) {
                    SetList.add(s);
                }
                Collections.sort(SetList);
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public String getTPM(String f) {
        StringBuilder sb = new StringBuilder();
        sb.append("x=read.table(\"" + f + "\",sep = \"\\t\")" + "\n");
        sb.append("countToTpm <- function(counts, effLen)" + "\n").append("{" + "\n");
        sb.append("\t" + "rate <- log(counts) - log(effLen)" + "\n");
        sb.append("\t" + "denom <- log(sum(exp(rate)))" + "\n");
        sb.append("\t" + "exp(rate - denom + log(1e6))" + "\n");
        sb.append("}" + "\n");
        sb.append("countTable<-data.frame(count=x[,6:101],length=x$End-x$Start)" + "\n");
        sb.append("for(i in 1:95){" + "\n");
        sb.append("\t" + "countTable[,i]<-with(countTable,countToTpm(countTable[,i],length))" + "\n");
        sb.append("}" + "\n");
        sb.append("write.table(countTable,file=\"" + f.replace(".txt", "TPM.txt") + "\",append = F,quote = F,sep = \"\\t\",eol = \"\\n\",row.names = F,col.names = T)");
        String statement = sb.toString();
        System.out.println(statement);
        return statement;
    }

    public void charDemo() {
        char a = 'A';
        char b = (char) (a + 1);
        int c = a + b;
        System.out.println(c);
        System.out.println(b);
        System.out.println(a + b);
        System.out.println("a + b is " + a + b);
    }

    public void addinfo() {
        String DSdir = "/data2/xiaohan/DS/S7";
        String genoDir = "/data2/xiaohan/sub7";
        String outputDir = "/data2/xiaohan/addinfo/S7";
//        String DSdir = "/data2/xiaohan/DS/S3";
//        String genoDir = "/data2/xiaohan/sub3";
//        String outputDir = "/data2/xiaohan/addinfo/S3";
//        String DSdir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/DS/S7";
//        String genoDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/sub7";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/addinfo/S7";
        HashSet<String> nameSet = new HashSet<String>();
        //S3
//          String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        //S7
        String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
        File[] fs = new File(DSdir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf");
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
//            String[] SampleName = fs[i].getName().split("\\.");
//            String Name = SampleName[0].replace(SampleName[0].substring(0, 3), "");
            String Name = fs[i].getName().split("\\.")[0].split("l")[1].split("D")[0];
            nameSet.add(Name);
            System.out.print(Name);
        }
        nameSet.stream().forEach((String p) -> {
            try {
                String[] temps = null;
                String temp = null;
                String[] tems = null;
                String[] temps1 = null;
                String temp1 = null;
                BufferedReader br = IOUtils.getTextReader(new File(DSdir, "col" + p + "DS.vcf").getAbsolutePath());
                BufferedReader br1 = IOUtils.getTextReader(new File(genoDir, "snp" + p + ".vcf").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "genotypes" + p + ".vcf").getAbsolutePath());
                for (int m = 0; m < 45; m++) {
                    temp1 = br1.readLine();
                }
                while ((temp = br.readLine()) != null) {
                    temp1 = br1.readLine();
                    temps1 = temp1.split("\t");
                    if (temp.startsWith("#")) {
                        temps = temp.split("\t");
                        bw.write("#");
                        //1.Chrom+2.position
                        for (int i = 0; i < 2; i++) {
                            tems = temps[i].split("]");
                            bw.write(tems[1] + "\t");
                        }
                        //3.ID
                        bw.write("ID" + "\t");
                        //4.REF+5.ALT
                        for (int i = 2; i < 4; i++) {
                            tems = temps[i].split("]");
                            bw.write(tems[1] + "\t");
                        }
                        //6.QUAL+7.FILTER+8.INFO+9.FORMAT
                        bw.write(temps1[5] + "\t" + temps1[6] + "\t" + temps1[7] + "\t" + temps1[8] + "\t");
                        //96Sample
                        for (int i = 0; i < 95; i++) {
                            bw.write(names[i] + "\t");
                        }
                        bw.write(names[95]);
                        //                    for (int i = 4; i < 179; i++) {
                        //                        bw.write(temps[i] + "\t");
                        //                    }
                        bw.newLine();
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
                        bw.write(temps1[5] + "\t" + temps1[6] + "\t" + "INFO" + "\t" + "DS" + "\t");
                        for (int i = 4; i < 99; i++) {
                            bw.write(temps[i] + "\t");
                        }
                        bw.write(temps[99]);
                        bw.newLine();
                    }
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void addinfobySample() {
        String DSdir = "/data2/xiaohan/DS";
        String genoDir = "/data2/xiaohan/sub3";
        String outputDir = "/data2/xiaohan/addinfo/S3";
//        String DSdir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/DS/";
//        String genoDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/sub";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/test/addinfo/";
        HashSet<String> nameSet = new HashSet<String>();
        //3
        String name = "B18-E002,B18-E007,B18-E008,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        //7
//        String name = "B18-E002,B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E035,B18-E038,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E251,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E348,B18-E355,B18-E356,B18-E357";
        String[] names = name.split(",");
        String[] temps = null;
        String temp = null;
        String[] tems = null;
        String[] temps1 = null;
        String temp1 = null;
        for (int Samplecount = 4; Samplecount < 5; Samplecount++) {
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "genotypes" + Samplecount + ".vcf").getAbsolutePath());
            BufferedReader br = IOUtils.getTextReader(new File(DSdir, "col" + Samplecount + "DS.vcf").getAbsolutePath());
            BufferedReader br1 = IOUtils.getTextReader(new File(genoDir, "snp" + Samplecount + ".vcf").getAbsolutePath());
            try {
                for (int m = 0; m < 45; m++) {
                    temp1 = br1.readLine();
                }
                while ((temp = br.readLine()) != null) {
                    temp1 = br1.readLine();
                    temps1 = temp1.split("\t");
                    if (temp.startsWith("#")) {
                        temps = temp.split("\t");
                        bw.write("#");
                        //1.Chrom+2.position
                        for (int i = 0; i < 2; i++) {
                            tems = temps[i].split("]");
                            bw.write(tems[1] + "\t");
                        }
                        //3.ID
                        bw.write("ID" + "\t");
                        //4.REF+5.ALT
                        for (int i = 2; i < 4; i++) {
                            tems = temps[i].split("]");
                            bw.write(tems[1] + "\t");
                        }
                        //6.QUAL+7.FILTER+8.INFO+9.FORMAT
                        bw.write(temps1[5] + "\t" + temps1[6] + "\t" + temps1[7] + "\t" + temps1[8] + "\t");
                        //96Sample
                        for (int i = 0; i < 95; i++) {
                            bw.write(names[i] + "\t");
                        }
                        bw.write(names[95]);
                        bw.newLine();
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
                        bw.write(temps1[5] + "\t" + temps1[6] + "\t" + "INFO" + "\t" + "DS" + "\t");
                        for (int i = 4; i < 99; i++) {
                            bw.write(temps[i] + "\t");
                        }
                        bw.write(temps[99]);
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
    }


    public void vcfmerge() {
        StringBuilder sb = new StringBuilder();
        sb.append("nohup vcf-concat ");
        for (int i = 0; i < 45; i++) {
            sb.append("/data3/wgs/vcf/GATK/vmap3/1.SNP/" + i + ".snp.vcf.gz ");
        }
        sb.append("| bgzip -c > all.non-filtering.vcf.gz &");
        System.out.println(sb.toString());
    }

    public void vcffiltering() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < 45; i++) {
//            sb.append("nohup vcftools --gzvcf /data3/wgs/vcf/GATK/vmap3/1.SNP/");
//            sb.append(i+".snp.vcf.gz --maf 0 --max-maf 0.05 --out ");
//            sb.append("/data2/xiaohan/SNP/");
//            sb.append(i+".snp.maf005 --recode "+" || ");
            sb.append("bgzip " + i + ".snp.maf005.recode.vcf" + "\n");
//            sb.append("log"+i+".txt 2>&1 &"+"\n");
        }
        String words = sb.toString();
        System.out.println(words);
    }

    public static void mkshell(File file) {
        File[] fl = file.listFiles();
        for (File f : fl) {
            if (f.isDirectory()) {
                mkshell(f);
            }
            if (f.isFile()) {
                System.out.println(f);
            }
        }
    }

    public void callposition() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/origin/wheat_v1.1_Lulab.gtf";
        String geneNameS = null;
        int gfIndex = 0;
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
            for (int i = 0; i < geneNames.length; i++) {
                String geneS = geneNames[i];
                System.out.println(geneChrMap.get(geneS) + "\t" + geneMinMap.get(geneS) + "\t" + geneMaxMap.get(geneS) + "\t" + geneS);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
            System.out.println(geneNameS + "\t" + gfIndex);
        }
    }

    public void findDifference() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/expressionboxcox.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedReader br1 = IOUtils.getTextReader(infileS);
        try {
            String temp = br.readLine();
            String temp1 = br1.readLine();
            temp = br.readLine();
            temp1 = br1.readLine();
            while ((temp = br.readLine()) != null) {
                String[] temps = temp.split("\t");
                String[] temps1 = temp1.split("\t");
                String chr = temps[0];
                String chr1 = temps1[0];
                String largerNumber = temps[1];
                String smallNumber = temps1[1];
                int larger = Integer.parseInt(largerNumber);
                int small = Integer.parseInt(smallNumber);
                if (chr.equals(chr1) && larger < small) {
                    System.out.println(temps1[0] + "\t" + temps1[1]);
                }
                temp1 = br1.readLine();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void findFalseSample() {
        String infileS = "/Users/yxh/Documents/eQTL/003experiment/浓度汇总/RNAconcentration.txt";
        BufferedReader br = IOUtils.getTextReader(infileS);
        String[] temps = null;
        String temp = null;
        HashSet<String> NameSet = new HashSet<String>();
        HashMap<String, Double> A280Map = new HashMap<>();
        HashMap<String, Double> A230Map = new HashMap<>();
        HashMap<String, Double> ConcentrationMap = new HashMap<>();
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                temps = temp.split("\t");
                String Name = temps[0];
                Double A280 = Double.valueOf(temps[1]);
                Double A230 = Double.valueOf(temps[2]);
                Double Concentration = Double.valueOf(temps[3]);
                NameSet.add(Name);
                A280Map.put(Name, A280);
                A230Map.put(Name, A230);
                ConcentrationMap.put(Name, Concentration);
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
                if ((A280 > 2.3 || A280 < 1.6) || (A230 > 2.5 || A230 < 1.6) || (Concentration < 30)) {
                    System.out.println(Name + "\t" + A280 + "\t" + A230 + "\t" + Concentration);
                }

            }

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void findTaxon() {

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

    public static void main(String[] args) throws IOException {
        new Test();
//        String path = "/data2/junxu/SiPASData/data_P101SC18112845-01-F004-B4-21/1.rawdata";
//        File file = new File(path);
//        mkshell(file);
    }
}

