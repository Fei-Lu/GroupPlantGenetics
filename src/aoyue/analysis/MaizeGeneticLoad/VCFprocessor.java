/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.table.RowTable;

/**
 *
 * @author Aoyue
 */
public class VCFprocessor {

    public VCFprocessor() {
        //this.getBasedSNP();
        //this.subsetBasedSNP();
        //this.mergeSubsetVCF();
        this.randomVCF();

    }
    
    public void randomVCF(){
        String infileS = "/Users/Aoyue/project/maizeGeneticLoad/maf_fei/popgen/group/random_chr10_hmp321.vcf.gz";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/013_structure/chr10_hmp321.subset.vcf.gz";
        //String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/013_structure/chr10_hmp321.subset.small.vcf.gz";
        try{
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > 0.3) {
                            //if(r > 0.001){
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
            }
            br.close();bw.flush();bw.close();
            System.out.println(cnt + "  snps are subseted");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * Goal: 将抽样的10条染色体合并成一条完整的VCF； Step1:先写表头；Step2:在按行读。
     */
    public void mergeSubsetVCF() {
        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/010_subsetVCF";
        String outfileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/011_mergesubsetVCF/hmp321_agpv4_merge.filt.subset.vcf.gz";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".filt.subset.vcf.gz");
        try {
            String headerfileS = new File(infileDirS, "hmp321_agpv4_chr001.filt.subset.vcf.gz").getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(headerfileS);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
            }
            br.close();
            int cnt = 0;
            for (int i = 0; i < fs.length; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i + 1);
                String infileS = new File(infileDirS, "hmp321_agpv4_chr" + chr + ".filt.subset.vcf.gz").getAbsolutePath();
                br = IOUtils.getTextGzipReader(infileS);
                int cntchr = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                    cnt++;
                    cntchr++;
                }
                System.out.println(chr + "  " + cntchr + "   sites.");
            }
            bw.flush();
            bw.close();
            System.out.println("Total SNPs is   " + " " + cnt + ".");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * this is the entrance of maize GLoadEntrance 001 1319 sites. 002 1079
         * sites. 003 1035 sites. 004 1070 sites. 005 992 sites. 006 790 sites.
         * 007 740 sites. 008 829 sites. 009 695 sites. 010 687 sites. Total
         * SNPs is 9236.
         */
    }

    /**
     * 根据生成的BasedSNP集合，进行抽样，合计抽取 10000 sites，来分析PCA和STRUCTURE.
     * 1.按照染色体长度来抽取每条染色体的sites数； 2.同时过滤含有2个alt的基因型。
     */
    public void subsetBasedSNP() {
        /**
         * *** HPC *****
         */
        String infileDirS = "/data1/home/aoyue/maizeGeneticLoad/001_variantSummary/009_vcf_filtration";
        String outfileDirS = "/data1/home/aoyue/maizeGeneticLoad/001_variantSummary/010_subsetVCF";

        /**
         * *** local test *****
         */
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/009_test_vcf_filtration";
//        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/010_test_subsetVCF";
//        String mergefileS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/011_test_mergeVCF/";
        int SNPsites = 10000; //待会儿要修改呀！！！！！
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".filt.vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        double ratio1 = 0.145770384;
        double ratio2 = 0.116050825;
        double ratio3 = 0.111885092;
        double ratio4 = 0.117262562;
        double ratio5 = 0.106299287;
        double ratio6 = 0.082623568;
        double ratio7 = 0.086587021;
        double ratio8 = 0.085989346;
        double ratio9 = 0.075851916;
        double ratio10 = 0.071679999;
        double snp1 = 9398581;
        double snp2 = 7227800;
        double snp3 = 6762066;
        double snp4 = 7299366;
        double snp5 = 6157059;
        double snp6 = 4728300;
        double snp7 = 4804064;
        double snp8 = 4891714;
        double snp9 = 4102494;
        double snp10 = 4304177;
        //1.确定玉米每条染色体长度占总长度的比例，然后进行均匀抽样。
        /**
         * Chromosome	Length(V4)	CentromereS	CentromereE	Ratio	SNPs 1	307041717
         * 136770000	137120000	0.145770384	9398581 2	244442276	95510000	97490000
         * 0.116050825	7227800 3	235667834	85780000	86930000	0.111885092	6762066
         * 4	246994605	109070000	110500000	0.117262562	7299366 5	223902240
         * 104540000	106820000	0.106299287	6157059 6	174033170	52300000	53110000
         * 0.082623568	4728300 7	182381542	56380000	56680000	0.086587021	4804064
         * 8	181122637	50530000	52070000	0.085989346	4891714 9	159769782
         * 53750000	57760000	0.075851916	4102494 10	150982314	51390000	52780000
         * 0.071679999	4304177 2106338117 1 59,675,621
         */

        fsList.parallelStream().forEach(p -> {
            try {
                String chr = p.getName().split("_chr")[1].split(".filt.vcf")[0];
                System.out.println("The chromosome we are processing is " + chr + ".");
                String outfileS = new File(outfileDirS, "hmp321_agpv4_chr" + chr + ".filt.subset.vcf.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextGzipReader(p.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                if (chr.equals("001")) { //need remodify
                    double dd = SNPsites * ratio1; //need remodify
                    double ratio = dd / snp1; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("002")) {
                    double dd = SNPsites * ratio2; //need remodify
                    double ratio = dd / snp2; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("003")) {
                    double dd = SNPsites * ratio3; //need remodify
                    double ratio = dd / snp3; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("004")) {
                    double dd = SNPsites * ratio4; //need remodify
                    double ratio = dd / snp4; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("005")) {
                    double dd = SNPsites * ratio5; //need remodify
                    double ratio = dd / snp5; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("006")) {
                    double dd = SNPsites * ratio6; //need remodify
                    double ratio = dd / snp6; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("007")) {
                    double dd = SNPsites * ratio7; //need remodify
                    double ratio = dd / snp7; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("008")) {
                    double dd = SNPsites * ratio8; //need remodify
                    double ratio = dd / snp8; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("009")) {
                    double dd = SNPsites * ratio9; //need remodify
                    double ratio = dd / snp9; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");

                } else if (chr.equals("010")) {
                    double dd = SNPsites * ratio10; //need remodify
                    double ratio = dd / snp10; //need remodify
                    int cnt = 0;
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = PStringUtils.fastSplit(temp);
                            if (l.get(4).contains(",")) {
                                continue;
                            }
                            if (l.get(4).startsWith("<")) {
                                continue;
                            }
                            // 第3列是alt的信息，若有2个等位基因，则去除这一行
                            double r = Math.random();
                            if (r > ratio) {
                                continue; //返回带正号的 double 值，该值大于等于 0.0 且小于 1.0。返回值是一个伪随机选择的数，在该范围内（近似）均匀分布
                            }
                            bw.write(temp);
                            bw.newLine();
                            cnt++;
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("**************** chr " + chr + " subset is finished****************\n" + chr + " " + cnt + "   sites was effectually keeped ");
                    System.out.println(chr + " " + dd + "   sites was planed to keeped ");
                }
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 通过技术过滤得到的位点库：/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter
     * 将原始hmp321_agpv4_chr*.vcf.gz 过滤，然后下一步进行其他群体遗传学分析。 先把vcf读进去，再进行多线程操作。
     */
    public void getBasedSNP() {
        /**
         * *** HPC *****
         */
        String infileDirS = "/data2/aoyue/maizeData/hmp321_agp4/";
        String dbfileDirS = "/data1/home/aoyue/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter/";
        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/009_test_vcf_filtration";

        /**
         * *** local test *****
         */
//        String infileDirS = "/Users/Aoyue/project/maizeGeneticLoad/oriData/hmp321_agpv4/";
//        String dbfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/001_hmp321Info_filter/";
//        String outfileDirS = "/Users/Aoyue/project/maizeGeneticLoad/001_variantSummary/009_vcf_filtration/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(p -> {
            try {
                String chr = p.getName().split("_chr")[1].split(".vcf")[0];
                chr = PStringUtils.getNDigitNumber(3, Integer.parseInt(chr));
                System.out.println("The chromosome we are processing is " + chr + ".");

                String dbfileS = new File(dbfileDirS, "hmp321Info_filter_chr" + chr + "_AGPv4_AnnoDB.txt").getAbsolutePath();
                String outfileS = new File(outfileDirS, "hmp321_agpv4_chr" + chr + ".filt.vcf.gz").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(dbfileS);
                ArrayList<Integer> posList = new ArrayList<>();  //将pos变成整型类型
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Chr")) {
                        continue;
                    }
                    temp = temp.substring(0, 40);
                    int pos = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                    posList.add(pos);
                }
                br.close();
                System.out.println("****************" + chr + " database is finished****************\n" + "chr " + chr + "  has pos count " + posList.size());

                Collections.sort(posList);
                br = IOUtils.getTextGzipReader(p.getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                int cnt = 0;
                while ((temp = br.readLine()) != null) {
                    StringBuilder sb = new StringBuilder();
                    if (temp.startsWith("#")) {
                        bw.write(temp);
                        bw.newLine();
                    } else {
                        int pos = Integer.parseInt(PStringUtils.fastSplit(temp).get(1));
                        int index = Collections.binarySearch(posList, pos);
                        if (index >= 0) {
                            cnt++;
                            sb.append(temp);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println("****************" + chr + " filtration is finished****************\n" + chr + cnt + "   sites was keeped ");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        //nohup java -jar getBasedSNP.jar > log_getBasedSNP.txt & 
    }
}
