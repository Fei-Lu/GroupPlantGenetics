package xiaohan.rareallele;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.math.BigDecimal;
import java.util.*;

public class RareAlleleCount {

    String Samplename = "B18-E007,B18-E008,B18-E011,B18-E014,B18-E018,B18-E023,B18-E024,B18-E029,B18-E032,B18-E038,B18-E043,B18-E046,B18-E049,B18-E051,B18-E052,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E083,B18-E087,B18-E089,B18-E097,B18-E099,B18-E118,B18-E124,B18-E127,B18-E138,B18-E139,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E203,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E227,B18-E228,B18-E233,B18-E237,B18-E242,B18-E245,B18-E252,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E320,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E347,B18-E355,B18-E356,B18-E357";
    String SNPfileDir = "/data2/xiaohan/genotype_root/genotype_rootMaf005";//根据不同MAF值过滤得到的VCF文件存放位置
    String TSSpositionfileDir = "/data1/home/xiaohan/rareallele/fastQTL/expression/S7";//存储不同基因的位置区域的文件存放位置chr + start + end + gene
    String rareAlleleCountDir = "/data1/home/xiaohan/rareallele/rankcorrelation/root/rareAlleleCount/Gerp200k";//不同的上游rare allele数目统计的文件存放位置
    String inforDir = "/data2/xiaohan/GerpOrigin/chr";

    String donor02GeneNamelistfile = "/data1/home/xiaohan/rareallele/rankcorrelation/infor/donor02GeneName.txt";//表达个体数目大于百分之二十的基因名字列表文件
    String homoGeneNamelistfile = "/data1/home/xiaohan/rareallele/rankcorrelation/infor/TheABD.txt";//ABD同源基因的基因列表文件

    String[] subDir = {"chr1-42", "all", "donor02", "sub", "subhomo"};
    //    String[] FileNames = {"0k_1000k"};
//    String[] FileNames = {"0k_2k", "2k_4k", "4k_6k", "6k_8k", "8k_10k", "0k_10k", "10k_12k", "12k_14k", "14k_16k", "16k_18k", "18k_20k", "10k_20k"};
//    String[] FileNames = {"0k_100k", "100k_200k", "200k_300k", "300k_400k", "400k_500k", "500k_600k", "600k_700k", "700k_800k", "800k_900k", "900k_1000k"};
//    String[] FileNames = {"0k_100k"};
    String[] FileNames = {"100k_200k", "200k_300k", "300k_400k", "400k_500k", "500k_600k", "600k_700k", "700k_800k", "800k_900k", "900k_1000k"};
//  String[] FileNames = {"200k_300k", "300k_400k"};
//  String[] FileNames = {"400k_500k", "500k_600k"};
//  String[] FileNames = {"600k_700k", "700k_800k"};
//  String[] FileNames = {"800k_900k", "900k_1000k"};
//  
//    String[] FileNames = {"0k_100k"};

    String mafDir = "/data2/xiaohan/VCF_information/frqMaf005";
    String annotationfile = "/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3";

    public RareAlleleCount() {
//        this.parseParameters(arg);
        this.mkFileDir();
        this.getupstreamSNPcount();//根据上游的不同位置call出rare allele count
//        this.getcodingSNPcount();
        this.mergeUprareCount();
        this.getdonor02File();
        this.getsubdonor02File();
    }

//    private void parseParameters(String infileS) {
//        List<String> pLineList = new ArrayList<>();
//        try {
//            BufferedReader br = IOUtils.getTextReader(infileS);
//            String temp = null;
//            boolean ifOut = false;
//            if (!(temp = br.readLine()).equals("RareAllele-Count")) {
//                ifOut = true;
//            }
//            if (!(temp = br.readLine()).equals("Author: Xiaohan Yang")) {
//                ifOut = true;
//            }
//            if (!(temp = br.readLine()).equals("Email: yangxh0120@gmail.com")) {
//                ifOut = true;
//            }
//            if (ifOut) {
//                System.out.println("Thanks for using RareAllele-Count.");
//                System.out.println("Please keep the authorship in the parameter file. Program stops.");
//                System.exit(0);
//            }
//            while ((temp = br.readLine()) != null) {
//                if (temp.startsWith("#")) {
//                    continue;
//                }
//                if (temp.isEmpty()) {
//                    continue;
//                }
//                pLineList.add(temp);
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(1);
//        }
////        this.Samplename=pLineList.get(0);
////        this.SNPfileDir=pLineList.get(1);
////        this.TSSpositionfileDir=pLineList.get(2);
////        this.rareAlleleCountDir=pLineList.get(3);
////        this.FileName=pLineList.get(4);
////        FileNames = FileName.split(",");
////        this.donor02GeneNamelistfile=pLineList.get(5);
//        this.homoGeneNamelistfile = pLineList.get(6);
//    }
    public void mkFileDir() {
        for (int i = 0; i < this.subDir.length; i++) {
            new File(this.rareAlleleCountDir, subDir[i]).mkdir();
        }
    }

    public void getsubdonor02File() {
        String inputDir = this.rareAlleleCountDir + "/" + subDir[2];
//        String inputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/root/expressionTable/donor02";
        String infor = this.homoGeneNamelistfile;
        String outputDirsub = this.rareAlleleCountDir + "/" + subDir[3];
        String outputDirhomo = this.rareAlleleCountDir + "/" + subDir[4];
//        String outputDirsub = "/data1/home/xiaohan/rareallele/rankcorrelation/root/expressionTable/donor02/sub";
//        String outputDirhomo = "/data1/home/xiaohan/rareallele/rankcorrelation/root/expressionTable/donor02/subhomo";

        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
//        fs = IOUtils.listFilesStartsWith(fs, "donor02");
//        String prefix = "DE";
        String prefix = "donor";
        fs = IOUtils.listFilesStartsWith(fs, prefix);
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName();
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String info = null;
                String[] infos = null;
                String temp = null;
                String[] temps = null;
                BufferedReader brinfo = IOUtils.getTextReader(infor);
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, f).getAbsolutePath());
//                BufferedWriter bwA = IOUtils.getTextWriter(new File(outputDirsub, f.replace("donor02", "A_donor02")).getAbsolutePath());
//                BufferedWriter bwB = IOUtils.getTextWriter(new File(outputDirsub, f.replace("donor02", "B_donor02")).getAbsolutePath());
//                BufferedWriter bwD = IOUtils.getTextWriter(new File(outputDirsub, f.replace("donor02", "D_donor02")).getAbsolutePath());
//                BufferedWriter bwAhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace("donor02", "Ahomo_donor02")).getAbsolutePath());
//                BufferedWriter bwBhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace("donor02", "Bhomo_donor02")).getAbsolutePath());
//                BufferedWriter bwDhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace("donor02", "Dhomo_donor02")).getAbsolutePath());
                BufferedWriter bwA = IOUtils.getTextWriter(new File(outputDirsub, f.replace(prefix, "A_" + prefix)).getAbsolutePath());
                BufferedWriter bwB = IOUtils.getTextWriter(new File(outputDirsub, f.replace(prefix, "B_" + prefix)).getAbsolutePath());
                BufferedWriter bwD = IOUtils.getTextWriter(new File(outputDirsub, f.replace(prefix, "D_" + prefix)).getAbsolutePath());
                BufferedWriter bwAhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace(prefix, "Ahomo_" + prefix)).getAbsolutePath());
                BufferedWriter bwBhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace(prefix, "Bhomo_" + prefix)).getAbsolutePath());
                BufferedWriter bwDhomo = IOUtils.getTextWriter(new File(outputDirhomo, f.replace(prefix, "Dhomo_" + prefix)).getAbsolutePath());
                HashSet<String> AsubGeneSet = new HashSet();
                HashSet<String> BsubGeneSet = new HashSet();
                HashSet<String> DsubGeneSet = new HashSet();
                while ((info = brinfo.readLine()) != null) {
                    infos = info.split("\t");
                    AsubGeneSet.add(infos[1]);
                    BsubGeneSet.add(infos[2]);
                    DsubGeneSet.add(infos[3]);
                }
                String title = br.readLine();
                bwA.write(title);
                bwB.write(title);
                bwD.write(title);
                bwA.newLine();
                bwB.newLine();
                bwD.newLine();
                bwAhomo.write(title);
                bwBhomo.write(title);
                bwDhomo.write(title);
                bwAhomo.newLine();
                bwBhomo.newLine();
                bwDhomo.newLine();
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    String Name = temps[0].substring(8, 9);
                    String Name1 = temps[0];
                    if (Name.equals("A")) {
                        bwA.write(temp);
                        bwA.newLine();
                    }
                    if (Name.equals("B")) {
                        bwB.write(temp);
                        bwB.newLine();
                    }
                    if (Name.equals("D")) {
                        bwD.write(temp);
                        bwD.newLine();
                    }
                    if (AsubGeneSet.contains(Name1)) {
                        bwAhomo.write(temp);
                        bwAhomo.newLine();
                    }
                    if (BsubGeneSet.contains(Name1)) {
                        bwBhomo.write(temp);
                        bwBhomo.newLine();
                    }
                    if (DsubGeneSet.contains(Name1)) {
                        bwDhomo.write(temp);
                        bwDhomo.newLine();
                    }
                }
                br.close();
                brinfo.close();
                bwA.flush();
                bwAhomo.flush();
                bwA.close();
                bwAhomo.close();
                bwB.flush();
                bwBhomo.flush();
                bwB.close();
                bwBhomo.close();
                bwD.flush();
                bwDhomo.flush();
                bwD.close();
                bwDhomo.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getdonor02File() {
        String inputDir = this.rareAlleleCountDir + "/" + subDir[1];

//        String infile = "/data1/home/xiaohan/rareallele/rankcorrelation/root/expressionTable/DEnorm7_87chr1-42.txt";
        String infor = this.donor02GeneNamelistfile;
        String outputDir = this.rareAlleleCountDir + "/" + subDir[2];
//        String outputDir = "/data1/home/xiaohan/rareallele/rankcorrelation/root/expressionTable/donor02";
        File[] fs = new File(inputDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesStartsWith(fs, "all");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName();
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String info = null;
                String temp = null;
                String[] temps = null;
                BufferedReader brinfo = IOUtils.getTextReader(infor);
                BufferedReader br = IOUtils.getTextReader(new File(inputDir, f).getAbsolutePath());
//            BufferedReader br = IOUtils.getTextReader(infile);
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, f.replace("all", "donor02")).getAbsolutePath());
//            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir,"DEnorm7_87chr1-42_donor02.txt").getAbsolutePath());
                HashSet<String> subGeneSet = new HashSet();
                String[] geneName = null;
                HashMap<String, String> tempGene = new HashMap<>();
                int countline = 0;
                while ((info = brinfo.readLine()) != null) {
                    subGeneSet.add(info);
                    countline++;
                }
                geneName = subGeneSet.toArray(new String[subGeneSet.size()]);
                Arrays.sort(geneName);
                String title = br.readLine();
                bw.write(title);
                bw.newLine();
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    tempGene.put(temps[0], temp);
                }
                for (int i = 0; i < geneName.length; i++) {
//                    System.out.println(geneName[i]);
                    bw.write(tempGene.get(geneName[i]));
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void mergeUprareCount() {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < FileNames.length; i++) {
            sb.append("cat chr*_" + FileNames[i] + "_RACount.txt | sort | uniq > " + this.rareAlleleCountDir + "/" + subDir[1] + "/all_" + FileNames[i] + "_RACount.txt \n");
            String command = sb.toString();
            try {
                File dir = new File(new File(this.rareAlleleCountDir, subDir[0]).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
//        StringBuilder sb1 = new StringBuilder();
//        sb1.append("cat chr*" + "_codingRACount.txt | sort | uniq > " + this.rareAlleleCountDir + "/" + subDir[1] + "/all" + "_codingRACount.txt \n");
//        String command1 = sb1.toString();
//        try {
//            File dir = new File(new File(this.rareAlleleCountDir, subDir[0]).getAbsolutePath());
//            String[] cmdarry = {"/bin/bash", "-c", command1};
//            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//            p.waitFor();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
    }

    public void getupstreamSNPcount() {
        String VCFfileDir = this.SNPfileDir;
        String outputDir = this.rareAlleleCountDir + "/" + subDir[0];
        String positionfileDir = this.TSSpositionfileDir;
//        S3
//        String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//        S7
        //建立对应的name位点表
        String[] Samplenames = Samplename.split(",");
        //对于多个区间"i"，多个chr"j"进行rare allele count计数
        System.out.println("This program is going to calculate minor allele count in every gene in different distance to TSS");
        System.out.println("The calculation form is : 1.Deciding the minor allele: 0 or 1");
        System.out.println("                          2.calculate the MAF and discard SNPs which MAF > 0.05");
        System.out.println("                          3.0/0 or 1/1 will be counting as 0 allele count or 2 allele count decided by which one is minor allele");
        System.out.println("                          4.Sum of every snp upstream genes");
        System.out.println("                          5.Output file form : gene B18-E007    B18-E009");
        System.out.println("                                               TraeCS4D355900 0 12");
        HashSet<String> nameSet = new HashSet();
        File[] fs = new File(VCFfileDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf.gz");
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName().split("\\.")[1].replace("chr", "");
//            if (21 <= Integer.parseInt(Name) && Integer.parseInt(Name) <= 42) {
                nameSet.add(Name);
                System.out.println(Name);
//            }
        }
        nameSet.parallelStream().forEach((String j) -> {
            for (int fn = 0; fn < FileNames.length; fn++) {
                try {
                    /*
                    如果需要Gerp文件或者是Maf文件
                     */
                    String base = null;
                    double maf1 = 0.0;
                    double maf2 = 0.0;
                    double maf = 0.0;
                    int chrNumber = Integer.parseInt(String.valueOf(j));
                    BufferedReader brinfo = IOUtils.getTextGzipReader(new File(this.inforDir, "chr" + j + ".bed.gz").getAbsolutePath());
                    //创建scoreMap
                    HashMap<Integer, Double> scoreMap = new HashMap();
                    Set scoreSet = new HashSet();
                    String info = null;
                    String[] infos = null;
                    while ((info = brinfo.readLine()) != null) {
                        infos = info.split("\t");
                        double score = Double.parseDouble(infos[4]);
                        score = new BigDecimal(score).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue();
                        scoreMap.put(Integer.parseInt(infos[2]), score);
                        scoreSet.add(Integer.parseInt(infos[2]));
                    }
                    System.out.println("Finished creating scoreMap………………………………");
                    //创建MAFMap
                    BufferedReader brmaf = IOUtils.getTextReader(new File(this.mafDir, "chr" + j + ".frq").getAbsolutePath());
                    System.out.println("This is reading maffile : " + new File(this.mafDir, "chr" + j + ".frq").getAbsolutePath());
                    HashMap<Integer, Double> mafMap = new HashMap();
                    HashMap<Integer, String> mafbaseMap = new HashMap();
                    String maftemp = null;
                    String[] maftemps = null;
                    while ((maftemp = brmaf.readLine()) != null) {
                        if (maftemp.startsWith("CHR")) {
                            continue;
                        }
                        maftemps = maftemp.split("\t");
                        maf1 = Double.parseDouble(maftemps[maftemps.length - 2].split(":")[1]);
                        maf2 = Double.parseDouble(maftemps[maftemps.length - 1].split(":")[1]);
                        if (maf1 <= maf2) {
//                            maf = new BigDecimal(maf1).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue();
                            base = maftemps[maftemps.length - 2].split(":")[0];
                        }
                        if (maf1 >= maf2) {
//                            maf = new BigDecimal(maf2).setScale(4, BigDecimal.ROUND_HALF_UP).doubleValue();
                            base = maftemps[maftemps.length - 1].split(":")[0];
                        }
//                        mafMap.put(Integer.parseInt(maftemps[1]), maf);
                        mafbaseMap.put(Integer.parseInt(maftemps[1]), base);
                    }
                    System.out.println("Finished creating mafMap………………………………");
                    //0.读入文件和创建输出文件
                    String ksite = FileNames[fn].split("_")[0];
                    System.out.println("This is running :" + "chr" + j + " file at " + ksite + " site");
                    BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "87B18.chr" + j + ".maf005.vcf.gz").getAbsolutePath());
                    BufferedReader brexpr = IOUtils.getTextGzipReader(new File(positionfileDir, "S7expression" + j + ".bed.gz").getAbsolutePath());
                    BufferedWriter bwS = IOUtils.getTextWriter(new File(outputDir, "chr" + j + "_" + FileNames[fn] + "_RACount.txt").getAbsolutePath());

                    //1.根据输出的矩阵SNPfile来对不同基因进行计数。
                    bwS.write("gene\t");
                    for (int m = 0; m < Samplenames.length; m++) {
                        bwS.write(Samplenames[m] + "\t");
                    }
                    bwS.newLine();
                    //2.创建一个基因list
                    GeneFeature gf = new GeneFeature(this.annotationfile);
                    HashSet<String> geneName = new HashSet<>();
                    for (int i = 0; i < gf.genes.length; i++) {
                        String name = gf.genes[i].geneName;
                        int chr = gf.genes[i].geneRange.chr;
                        if (chr == chrNumber) {
                            geneName.add(name);
                        }
                    }
                    String[] geneNames = geneName.toArray(new String[geneName.size()]);
                    Arrays.sort(geneNames);
                    //3.创建计数矩阵并将其初始化为0
                    int[][] count = new int[geneNames.length][Samplenames.length];
                    for (int m = 0; m < geneNames.length; m++) {
                        for (int n = 0; n < Samplenames.length; n++) {
                            count[m][n] = 0;
                        }
                    }
                    System.out.println("Finished initializing count matrix");
                    //4.对VCF文件读入每一行进行处理输出MAF<0.05的SNPfile矩阵
                    String[] VCF = null;
                    String tempVCF = null;
                    //4.1根据FileName定下所筛选的TSS上游的位置范围
                    int DiscontrolS = Integer.parseInt(FileNames[fn].split("_")[0].replace("k", "")) * 1000;
                    int DiscontrolE = Integer.parseInt(FileNames[fn].split("_")[1].replace("k", "")) * 1000;
                    System.out.println(DiscontrolE + "\t" + DiscontrolS);
                    HashMap<String, Integer> VCFMap = new HashMap();
                    int countline = 0;
                    while ((tempVCF = brVCF.readLine()) != null) {
                        countline++;
                        if (countline % 5000 == 0) {
                            System.out.print(countline + "\n");
                        }
                        if (tempVCF.startsWith("##")) {
                            continue;
                        }
                        if (tempVCF.startsWith("#C")) {
                            VCF = tempVCF.split("\t");
                            for (int m = 9; m < VCF.length; m++) {
                                VCFMap.put(VCF[m], m);
                            }
                            continue;
                        } else {
                            VCF = tempVCF.split("\t");
                            int snpsite = Integer.parseInt(VCF[1]);
                            int TSSrow = Integer.MAX_VALUE;
                            for (int l = 0; l < geneNames.length; l++) {
                                int strand = gf.getGeneStrand(gf.getGeneIndex(geneNames[l]));
                                int startsite = gf.getGeneStart(gf.getGeneIndex(geneNames[l]));
                                int endsite = gf.getGeneEnd(gf.getGeneIndex(geneNames[l]));
                                int distance = Integer.MAX_VALUE;
                                if (strand == 1) {
                                    distance = startsite - snpsite;
                                } else {
                                    distance = snpsite - endsite;
                                }
                                if(distance > DiscontrolE)continue;
                                if (distance >= DiscontrolS && distance < DiscontrolE) {
                                    String[] VCFforGT = new String[VCF.length - 9];
                                    for (int m = 0; m < VCFforGT.length; m++) {
                                        VCFforGT[m] = VCF[m + 9].split(":")[0];
                                    }
                                    double score = 0.000;
                                    if (scoreSet.contains(snpsite)) {
                                        score = scoreMap.get(snpsite);
                                        if (score > 1 && VCF[3].equals(mafbaseMap.get(snpsite))) {
//                                    if (VCF[3].equals(mafbaseMap.get(snpsite))) {
                                        for (int m = 0; m < Samplenames.length; m++) {
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/0")) {
                                                count[l][m] += 2;
                                            }
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("1/1")) {
                                                count[l][m] += 0;
                                            }
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/1")) {
                                                count[l][m] += 1;
                                            }
                                        }
                                        } else if (score > 1 && VCF[4].equals(mafbaseMap.get(snpsite))) {
//                                    } else if (VCF[4].equals(mafbaseMap.get(snpsite))) {
                                        for (int m = 0; m < Samplenames.length; m++) {
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/0")) {
                                                count[l][m] += 0;
                                            }
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("1/1")) {
                                                count[l][m] += 2;
                                            }
                                            if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/1")) {
                                                count[l][m] += 1;
                                            }
                                        }
                                        }
                                    }
                                } else if(distance < 0){
                                    break;
                                }
                            }
                        }
                    }
                    for (int m = 0; m < geneNames.length; m++) {
                        bwS.write(geneNames[m] + "\t");
                        for (int n = 0; n < Samplenames.length; n++) {
                            bwS.write(count[m][n] + "\t");
                        }
                        bwS.newLine();
                    }
                    System.out.println("complete writing file : " + "chr" + chrNumber + "at " + ksite + " site");
                    brexpr.close();
                    brVCF.close();
                    bwS.flush();
                    bwS.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            }
        });
    }

    public void getcodingSNPcount() {
        String VCFfileDir = this.SNPfileDir;
        String outputDir = this.rareAlleleCountDir + "/" + subDir[0];
        String positionfileDir = this.TSSpositionfileDir;
//        S3
//        String name = "B18-E007,B18-E010,B18-E011,B18-E014,B18-E016,B18-E018,B18-E029,B18-E032,B18-E035,B18-E043,B18-E045,B18-E046,B18-E049,B18-E051,B18-E062,B18-E065,B18-E070,B18-E072,B18-E074,B18-E081,B18-E082,B18-E089,B18-E097,B18-E099,B18-E115,B18-E118,B18-E124,B18-E127,B18-E134,B18-E138,B18-E141,B18-E152,B18-E166,B18-E170,B18-E180,B18-E184,B18-E185,B18-E188,B18-E199,B18-E204,B18-E205,B18-E210,B18-E214,B18-E215,B18-E218,B18-E219,B18-E228,B18-E233,B18-E236,B18-E237,B18-E242,B18-E244,B18-E245,B18-E252,B18-E253,B18-E256,B18-E262,B18-E265,B18-E267,B18-E270,B18-E271,B18-E273,B18-E277,B18-E280,B18-E286,B18-E288,B18-E289,B18-E290,B18-E298,B18-E299,B18-E305,B18-E306,B18-E312,B18-E316,B18-E318,B18-E324,B18-E330,B18-E332,B18-E335,B18-E337,B18-E346,B18-E348,B18-E355,B18-E356,B18-E357";
//        S7
        //建立对应的name位点表
        String[] Samplenames = Samplename.split(",");
//
        //对于多个区间"i"，多个chr"j"进行rare allele count计数
        System.out.println("This program is going to calculate minor allele count in every gene in different distance to TSS");
        System.out.println("The calculation form is : 1.Deciding the minor allele: 0 or 1");
        System.out.println("                          2.calculate the MAF and discard SNPs which MAF > 0.05");
        System.out.println("                          3.0/0 or 1/1 will be counting as 0 allele count or 2 allele count decided by which one is minor allele");
        System.out.println("                          4.Sum of every snp upstream genes");
        System.out.println("                          5.Output file form : gene B18-E007    B18-E009");
        System.out.println("                                               TraeCS4D355900 0 12");
        HashSet<String> nameSet = new HashSet();
        File[] fs = new File(VCFfileDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "vcf.gz");
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName().split("\\.")[1].replace("chr", "");
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach((String j) -> {
            try {
                //0.读入文件和创建输出文件
                BufferedReader brVCF = IOUtils.getTextGzipReader(new File(VCFfileDir, "87B18.chr" + j + ".maf005.vcf.gz").getAbsolutePath());
                BufferedReader brexpr = IOUtils.getTextGzipReader(new File(positionfileDir, "S7expression" + j + ".bed.gz").getAbsolutePath());
                BufferedWriter bwS = IOUtils.getTextWriter(new File(outputDir, "chr" + j + "_codingRACount.txt").getAbsolutePath());
                //1.根据输出的矩阵SNPfile来对不同基因进行计数。
                bwS.write("gene\t");
                for (int m = 0; m < Samplenames.length; m++) {
                    bwS.write(Samplenames[m] + "\t");
                }
                bwS.newLine();
                //2.创建一个关于基因和转录起始位点与终止位点的对应关系表
                String tempexpr = null;
                String[] expr = null;
                HashMap<String, String[]> TSSMap = new HashMap<>();
                Set<String> TSSset = new HashSet<>();
                while ((tempexpr = brexpr.readLine()) != null) {
                    if (tempexpr.startsWith("#")) {
                        continue;
                    }
                    List<String> tList = pgl.infra.utils.PStringUtils.fastSplit(tempexpr);
                    expr = tList.toArray(new String[tList.size()]);
                    String[] site = new String[2];
                    site[0] = expr[1];
                    site[1] = expr[2];
                    TSSMap.put(expr[3], site);
                    TSSset.add(expr[3]);
                    //TSSMap --> gene -- start
                    //TSSset --> gene
                }
                System.out.println("Finish creating gene site map");
                String[] TSS = TSSset.toArray(new String[TSSset.size()]);
                Arrays.sort(TSS);
                //3.创建计数矩阵并将其初始化为0
                double[][] count = new double[TSS.length][Samplenames.length];
                for (int m = 0; m < TSS.length; m++) {
                    for (int n = 0; n < Samplenames.length; n++) {
                        count[m][n] = 0;
                    }
                }
                System.out.println("Finished initializing count matrix");
                //4.对VCF文件读入每一行进行处理输出MAF<0.05的SNPfile矩阵
                String[] VCF = null;
                String tempVCF = null;
                //4.1根据FileName定下所筛选的TSS上游的位置范围
                HashMap<String, Integer> VCFMap = new HashMap();
                int countline = 0;
                while ((tempVCF = brVCF.readLine()) != null) {
                    countline++;
                    if (countline % 5000 == 0) {
                        System.out.print(countline + "\n");
                    }
                    if (tempVCF.startsWith("##")) {
                        continue;
                    }
                    if (tempVCF.startsWith("#C")) {
                        VCF = tempVCF.split("\t");
                        for (int m = 9; m < VCF.length; m++) {
                            VCFMap.put(VCF[m], m);
                        }
                        continue;
                    } else {
                        VCF = tempVCF.split("\t");
                        int snpsite = Integer.parseInt(VCF[1]);
                        int TSSrow = Integer.MAX_VALUE;
                        for (int l = 0; l < TSS.length; l++) {
                            int startsite = Integer.parseInt(TSSMap.get(TSS[l])[0]);
                            int stopsite = Integer.parseInt(TSSMap.get(TSS[l])[1]);
                            if (snpsite >= startsite && snpsite <= stopsite) {
                                TSSrow = l;
                                String[] VCFforGT = new String[VCF.length - 9];
                                for (int m = 0; m < VCFforGT.length; m++) {
                                    VCFforGT[m] = VCF[m + 9].split(":")[0];
                                }
                                int site0 = 0;
                                int site1 = 0;
                                for (int m = 0; m < VCFforGT.length; m++) {//读出0和1的个数，按照小值计算MAF。
                                    if (VCFforGT[m].equals("0/0")) {
                                        site0++;
                                        site0++;
                                    }
                                    if (VCFforGT[m].equals("0/1")) {
                                        site0++;
                                        site1++;
                                    }
                                    if (VCFforGT[m].equals("1/1")) {
                                        site1++;
                                        site1++;
                                    }
                                }
//                                     double MAF = 0.0000;
                                int minorAllele = Integer.MAX_VALUE;
                                int sitesum = site0 + site1;
                                if (site0 > site1) {
                                    minorAllele = 1;
                                }
                                if (site1 > site0) {
                                    minorAllele = 0;
                                }
                                /*
                                    获取对应的Gerp与Maf
                                 */
//                                    double score = 0.000;
//                                    if (scoreSet.contains(snpsite)) {
//                                        score = scoreMap.get(snpsite);
//                                        if (score > 0) {
//                                            score = score;
//                                        }
//                                        if (score <= 0) {
//                                            score = 0.00;
//                                        }
//                                    }
//                                    double MAF = mafMap.get(snpsite);
                                if (minorAllele == 0) {
                                    for (int m = 0; m < Samplenames.length; m++) {
                                        if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/0")) {
//                                                count[TSSrow][m] += 2 * score * (1 / MAF) / 10;
                                            count[TSSrow][m] += 2;
                                        }
                                        if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("1/1")) {
                                            count[TSSrow][m] += 0;
                                        }
                                    }
                                }
                                if (minorAllele == 1) {
                                    for (int m = 0; m < Samplenames.length; m++) {
                                        if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("0/0")) {
                                            count[TSSrow][m] += 0;
                                        }
                                        if (VCFforGT[VCFMap.get(Samplenames[m]) - 9].equals("1/1")) {
//                                                count[TSSrow][m] += 2 * score * (1 / MAF) / 10;
                                            count[TSSrow][m] += 2;
                                        }
                                    }
                                }
                            } else {
                                continue;
                            }
                        }
                    }
                }
                for (int m = 0; m < TSS.length; m++) {
                    bwS.write(TSS[m] + "\t");
                    for (int n = 0; n < Samplenames.length; n++) {
                        bwS.write(count[m][n] + "\t");
                    }
                    bwS.newLine();
                }
                System.out.println("complete writing file : " + "chr" + j);
                brexpr.close();
                brVCF.close();
                bwS.flush();
                bwS.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public static void main(String[] args) {

        new RareAlleleCount();
    }
}
