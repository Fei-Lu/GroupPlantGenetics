/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import pgl.infra.table.RowTable;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import java.awt.Point;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Collections;
import java.util.List;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
//import static org.apache.spark.ml.tree.DecisionTreeModel$class.depth;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Wheat120bamProcessor {

    public Wheat120bamProcessor() {
        //this.mkFastCallTaxaBamMap();
        //this.mkFastCallTaxaBamMapchr1_42();
        //this.mkParameterchr1_42();
        this.mkJavaCmdchr1_42();
        //this.mergewheatchrVCF();
        //this.countVCF();
        //this.subsetVCF();
        /**
         * ***只含有SD Depth*******
         */
        //this.statVcfCoverage();
        //this.subsetCoveVSSD();
        //this.calEllipse();
        /**
         * ***含有PV SD Depth*******
         */
        //this.statVcfPValue();
        //this.subsetCovevsSDvsPV();
        //this.calEllipse2();
        //this.calEllipse_onechr();
        this.filterVCF();

    }

    /**
     * 目的：根据已做的测试，在chr1A chr1B chr1D染色体中对选定的椭圆进行过滤。 
     * 1.计算chr1A的平均覆盖度和SD，输出到 CHROM POS	AverageDepth	SD 文件中；多线程运行！
     * 2.根据散点图中和椭圆的2a位置判断，将椭圆内需要保留的位点放入一个List
     * posinEllipseList库中; 单线程，考虑到每个椭圆的大小位置不一样
     * 3.将chr1A文件读进去，根究posinEllipseList库进行过滤筛选；
     */
    public void filterVCF() {
        //this.filterVCF_step1();
        this.filterVCF_step2();

    }

    /**
     * 计算点到椭圆2a的距离，并进行posinEllipseList库的建立 将VCF读入文件，进行筛选！
     */
    public void filterVCF_step2() {
        // local test
//        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/007_filterVCF_step1/";
//        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/008_filterVCF_step2/";
//        String vcfFileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/000_subSet/";
        //HPC
        String infileDirS = "/data2/aoyue/output_wheatVMapII/001_header/";
        String outfileDirS = "/data2/aoyue/output_wheatVMapII/002_ellipseFilter/";
        String vcfFileDirS = "/data2/sharedData/wheat_Jiao/vcf/vcf_chr1Ato7D/";

        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".header.txt.gz");
        for (int i = 0; i < fs.length; i++) {
            if (i == 0) {
                long startTime = System.nanoTime();
                double a = 2; //**************remodify
                double b = 0.58; //**************remodify
                double h = 9.8; //**************remodify
                double k = 3.6; //**************remodify
                double t = 9; //**************remodify
                //local test
//                String infileS = new File(infileDirS, "chr1A_subset.header.txt.gz").getAbsolutePath();
//                String outfileS = new File(outfileDirS, "chr1A.filt.vcf.gz").getAbsolutePath();
//                String vcfFileS = new File(vcfFileDirS, "chr1A_subset.vcf.gz").getAbsolutePath();
                //HPC run
                String infileS = new File(infileDirS,"chr1A.header.txt.gz").getAbsolutePath(); //**************remodify
                String outfileS = new File(outfileDirS,"chr1A.filt.vcf.gz").getAbsolutePath(); //**************remodify
                String vcfFileS = new File(vcfFileDirS,"chr1A.vcf.gz").getAbsolutePath(); //**************remodify
                double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
                System.out.println("The length c is  " + c);
                double sint = Math.sin((Math.PI * t / 180));
                double cost = Math.cos((Math.PI * t / 180));
                double ccost = c * cost;
                double csint = c * sint;
                double f1x = h - (c * cost);
                double f1y = k - (c * sint);
                double f2x = h + (c * cost);
                double f2y = k + (c * sint);
                System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
                System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
                int cntinellipse = 0;
                List<String> posinEllipseList = new ArrayList<>();
                try {
                    BufferedReader br = IOUtils.getTextGzipReader(infileS);
                    //CHROM	POS	REF	ALT	AverageDepth	SD	PValue
                    String temp = br.readLine();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        l = PStringUtils.fastSplit(temp, "\t");
                        String pos = l.get(1);
                        double x = Double.parseDouble(l.get(4));
                        double y = Double.parseDouble(l.get(5));
                        double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                        double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                        double distance = distance1 + distance2;
                        //System.out.println("两点间的距离是:" + distance);
                        double a2 = 2 * a;
                        if (distance <= a2) {
                            posinEllipseList.add(pos);
                            cntinellipse++;
                        }
                    }
                    System.out.println(cntinellipse + " in ellipse on chr1A"); //**************remodify
                    br.close();
                    System.out.println("****************chr1A database is finished****************"); //**************remodify
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

                try {
                    System.out.println("************Starting chr1A filter****************"); //**************remodify
                    BufferedReader br = IOUtils.getTextGzipReader(vcfFileS);
                    BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                    String temp = null;
                    int cnt = 0;
                    String[] posinEllipse = posinEllipseList.toArray(new String[posinEllipseList.size()]);
                    Arrays.sort(posinEllipse);
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        StringBuilder sb = new StringBuilder();
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            l = PStringUtils.fastSplit(temp, "\t");
                            String pos = l.get(1);
                            int index = Arrays.binarySearch(posinEllipse, pos);
                            if (index >= 0) {
                                sb.append(temp);
                                bw.write(sb.toString());
                                bw.newLine();
                                cnt++;
                            }
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("****************Filter chr1A is finished****************"); //**************remodify
                    System.out.println(cnt + " in ellipse on chr1A"); //**************remodify
                    long endTime = System.nanoTime();
                    float excTime = (float) (endTime - startTime) / 1000000000;
                    System.out.println("chr1A Execution time: " + String.format("%.2f", excTime) + "s");

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            } else if (i == 1) {
                long startTime = System.nanoTime();
                double a = 2.8;
                double b = 1;
                double h = 9.8;
                double k = 3.7;
                double t = 9;
                //local test
//                String infileS = new File(infileDirS, "chr1A_subset.header.txt.gz").getAbsolutePath();
//                String outfileS = new File(outfileDirS, "chr1A.filt.vcf.gz").getAbsolutePath();
//                String vcfFileS = new File(vcfFileDirS, "chr1A_subset.vcf.gz").getAbsolutePath();
                //HPC run
                String infileS = new File(infileDirS, "chr1B.header.txt.gz").getAbsolutePath(); //**************remodify
                String outfileS = new File(outfileDirS, "chr1B.filt.vcf.gz").getAbsolutePath(); //**************remodify
                String vcfFileS = new File(vcfFileDirS, "chr1B.vcf.gz").getAbsolutePath(); //**************remodify
                double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
                System.out.println("The length c is  " + c);
                double sint = Math.sin((Math.PI * t / 180));
                double cost = Math.cos((Math.PI * t / 180));
                double ccost = c * cost;
                double csint = c * sint;
                double f1x = h - (c * cost);
                double f1y = k - (c * sint);
                double f2x = h + (c * cost);
                double f2y = k + (c * sint);
                System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
                System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
                int cntinellipse = 0;
                List<String> posinEllipseList = new ArrayList<>();
                try {
                    BufferedReader br = IOUtils.getTextGzipReader(infileS);
                    //CHROM	POS	REF	ALT	AverageDepth	SD	PValue
                    String temp = br.readLine();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        l = PStringUtils.fastSplit(temp, "\t");
                        String pos = l.get(1);
                        double x = Double.parseDouble(l.get(4));
                        double y = Double.parseDouble(l.get(5));
                        double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                        double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                        double distance = distance1 + distance2;
                        //System.out.println("两点间的距离是:" + distance);
                        double a2 = 2 * a;
                        if (distance <= a2) {
                            posinEllipseList.add(pos);
                            cntinellipse++;
                        }
                    }
                    System.out.println(cntinellipse + " in ellipse on chr1B"); //**************remodify
                    br.close();
                    System.out.println("****************chr1B database is finished****************"); //**************remodify
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

                try {
                    System.out.println("************Starting chr1B filter****************"); //**************remodify
                    BufferedReader br = IOUtils.getTextGzipReader(vcfFileS);
                    BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                    String temp = null;
                    int cnt = 0;
                    String[] posinEllipse = posinEllipseList.toArray(new String[posinEllipseList.size()]);
                    Arrays.sort(posinEllipse);
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        StringBuilder sb = new StringBuilder();
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            l = PStringUtils.fastSplit(temp, "\t");
                            String pos = l.get(1);
                            int index = Arrays.binarySearch(posinEllipse, pos);
                            if (index >= 0) {
                                sb.append(temp);
                                bw.write(sb.toString());
                                bw.newLine();
                                cnt++;
                            }
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("****************Filter chr1B is finished****************"); //**************remodify
                    System.out.println(cnt + " in ellipse on chr1B"); //**************remodify
                    
                    long endTime = System.nanoTime();
                    float excTime = (float) (endTime - startTime) / 1000000000;
                    System.out.println("chr1B Execution time: " + String.format("%.2f", excTime) + "s");

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

            } else if (i == 2) {
                long startTime = System.nanoTime();
                double a = 2.6;
                double b = 1.2;
                double h = 9.8;
                double k = 3.6;
                double t = 9;
                //local test
//                String infileS = new File(infileDirS, "chr1D_subset.header.txt.gz").getAbsolutePath();
//                String outfileS = new File(outfileDirS, "chr1D.filt.vcf.gz").getAbsolutePath();
//                String vcfFileS = new File(vcfFileDirS, "chr1D_subset.vcf.gz").getAbsolutePath();
                //HPC run
                String infileS = new File(infileDirS, "chr1D.header.txt.gz").getAbsolutePath(); //**************remodify
                String outfileS = new File(outfileDirS, "chr1D.filt.vcf.gz").getAbsolutePath(); //**************remodify
                String vcfFileS = new File(vcfFileDirS, "chr1D.vcf.gz").getAbsolutePath(); //**************remodify
                double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
                System.out.println("The length c is  " + c);
                double sint = Math.sin((Math.PI * t / 180));
                double cost = Math.cos((Math.PI * t / 180));
                double ccost = c * cost;
                double csint = c * sint;
                double f1x = h - (c * cost);
                double f1y = k - (c * sint);
                double f2x = h + (c * cost);
                double f2y = k + (c * sint);
                System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
                System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
                int cntinellipse = 0;
                List<String> posinEllipseList = new ArrayList<>();
                try {
                    BufferedReader br = IOUtils.getTextGzipReader(infileS);
                    //CHROM	POS	REF	ALT	AverageDepth	SD	PValue
                    String temp = br.readLine();
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        l = PStringUtils.fastSplit(temp, "\t");
                        String pos = l.get(1);
                        double x = Double.parseDouble(l.get(4));
                        double y = Double.parseDouble(l.get(5));
                        double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                        double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                        double distance = distance1 + distance2;
                        //System.out.println("两点间的距离是:" + distance);
                        double a2 = 2 * a;
                        if (distance <= a2) {
                            posinEllipseList.add(pos);
                            cntinellipse++;
                        }
                    }
                    System.out.println(cntinellipse + " in ellipse on chr1D"); //**************remodify
                    br.close();
                    System.out.println("****************chr1D database is finished****************"); //**************remodify
                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }

                try {
                    System.out.println("************Starting chr1D filter****************"); //**************remodify
                    BufferedReader br = IOUtils.getTextGzipReader(vcfFileS);
                    BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
                    String temp = null;
                    int cnt = 0;
                    String[] posinEllipse = posinEllipseList.toArray(new String[posinEllipseList.size()]);
                    Arrays.sort(posinEllipse);
                    while ((temp = br.readLine()) != null) {
                        List<String> l = new ArrayList<>();
                        StringBuilder sb = new StringBuilder();
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            l = PStringUtils.fastSplit(temp, "\t");
                            String pos = l.get(1);
                            int index = Arrays.binarySearch(posinEllipse, pos);
                            if (index >= 0) {
                                sb.append(temp);
                                bw.write(sb.toString());
                                bw.newLine();
                                cnt++;
                            }
                        }
                    }
                    br.close();
                    bw.flush();
                    bw.close();
                    System.out.println("****************Filter chr1D is finished****************"); //**************remodify
                    System.out.println(cnt + " in ellipse on chr1D"); //**************remodify
                    long endTime = System.nanoTime();
                    float excTime = (float) (endTime - startTime) / 1000000000;
                    System.out.println("chr1D Execution time: " + String.format("%.2f", excTime) + "s");

                } catch (Exception e) {
                    e.printStackTrace();
                    System.exit(1);
                }
            }
        }
    }

    /**
     * mk header 并过滤有2个alt的位点
     */
    public void filterVCF_step1() {
//        String infileDirS = "/data2/sharedData/wheat_Jiao/vcf/vcf_chr1Ato7D/";
//        String outfileDirS = "/data2/aoyue/output_wheatVMapII/001_header/";
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/000_subSet";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/007_filterVCF_step1/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            //开始处理文件，计时
            long startTime = System.nanoTime();
            Calendar cal = Calendar.getInstance();
            //cal.set(Calendar.HOUR_OF_DAY, 23);
            int hour = cal.get(Calendar.HOUR_OF_DAY);
            int minute = cal.get(Calendar.MINUTE);
            int second = cal.get(Calendar.SECOND);
            System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + "Now starting to " + f.getName() + ";");

            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".header.txt.gz")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                //CHROM POS REF ALT AverageDepth SD PValue 
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + l.get(3) + "\t" + l.get(4) + "\t"
                        + "AverageDepth\tSD\tPValue");
                bw.newLine();
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();
                    cnt++; // 对snp开始计数
                    if (cnt % 1000000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    if (l.get(3).contains(",")) {
                        continue;
                    }
                    String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                    String ifPVzero = null;
                    String[] taxa = new String[l.size() - 9];
                    for (int i = 0; i < taxa.length; i++) {
                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)).append("\t")
                            .append(relativeMean).append("\t").append(sd)
                            .append("\t").append(pvalue);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                //文件处理完毕，计时
                hour = cal.get(Calendar.HOUR_OF_DAY);
                minute = cal.get(Calendar.MINUTE);
                second = cal.get(Calendar.SECOND);
                System.out.println(String.format("%d:%d:%d\t", hour, minute, second) + f.getName() + " is finished!!!");
                long endTime = System.nanoTime();
                float excTime = (float) (endTime - startTime) / 1000000000;
                System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 目的：进行单条染色体的测试；由于每个染色体的抽样椭圆分布可能不能，要一一进行测试。 根据计算好的pv sd
     * depth文件来过滤椭圆,在文件键入2列 CHROM	POS	AverageDepth	SD	Distance	Group	PValue
     * IfPVzero	AFG-L1
     */
    public void calEllipse_onechr() {
        // 进行chr1B测试
//        String infileS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/005_pvCalSample/chr1B_subset.depth.pv_5000sites.txt";
//        String outfileS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/006_addellipse2/chr1B_subset.depth.pv_5000sites.addellipse2.txt";
        // 进行chr1D测试
        String infileS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/005_pvCalSample/chr1D_subset.depth.pv_5000sites.txt";
        String outfileS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/006_addellipse2/chr1D_subset.depth.pv_5000sites.addellipse2.txt";

//        double a = 2.8;
//        double b = 1;
//        double h = 9.8;
//        double k = 3.7;
//        double t = 9;
        double a = 2.6;
        double b = 1.2;
        double h = 9.8;
        double k = 3.6;
        double t = 9;
        //测试用double a = 5 ; double b = 3; double h = 10; double k = 3.7 ; double t = 30 ; 
        double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
        System.out.println("The length c is  " + c);
        double sint = Math.sin((Math.PI * t / 180));
        double cost = Math.cos((Math.PI * t / 180));
        double ccost = c * cost;
        double csint = c * sint;
//        System.out.println("c " + c);
//        System.out.println("sint  " + sint);
//        System.out.println("cost  " + cost);
//        System.out.println("csint  " + csint);
//        System.out.println("ccost  " + ccost);
        double f1x = h - (c * cost);
        double f1y = k - (c * sint);
        double f2x = h + (c * cost);
        double f2y = k + (c * sint);
        System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
        System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
        //设置两个焦点的坐标p1 p2
        Point p1 = new Point();
        p1.setLocation(f1x, f1y);
        Point p2 = new Point();
        p2.setLocation(f2x, f2y);

        System.out.println("p1 " + p1.getX() + "    " + p1.getY());
        System.out.println("p2 " + p2.getX() + "    " + p2.getY());
        //定位坐标
//            System.out.println("p1的x坐标为"+p1.getX());
//            System.out.println("p1的y坐标为"+p1.getY());
//            System.out.println("p2的x坐标为"+p2.getX());
//            System.out.println("p2的y坐标为"+p2.getY());
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
        int cntinellipse = 0;
        try {
            //先处理表头
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            l = PStringUtils.fastSplit(temp);
            bw.write(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + "Distance" + "\t" + "Group");
            String[] taxa = new String[l.size() - 4];
            for (int i = 0; i < taxa.length; i++) {
                bw.write("\t" + l.get(i + 4));
            }
            bw.newLine();
            /*
            椭圆内分组为0，椭圆外分组为1.
             */
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp, "\t");
                double x = Double.parseDouble(l.get(2));
                double y = Double.parseDouble(l.get(3));
                //Point p = new Point();
                //p.setLocation(x, y);
                // 计算两点间距离公式 (此方法有误)
                //double distance1 = Math.sqrt(Math.abs((p1.getX() - p.getX())* (p1.getX() - p.getX())+(p1.getY() - p.getY())* (p1.getY() - p.getY())));
                //double distance2 = Math.sqrt(Math.abs((p2.getX() - p.getX())* (p2.getX() - p.getX())+(p2.getY() - p.getY())* (p2.getY() - p.getY())));
                double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                double distance = distance1 + distance2;
                System.out.println("两点间的距离是:" + distance);

                //开始写出文件
                StringBuilder sb = new StringBuilder();
                sb.append(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + String.format("%.4f", distance) + "\t");
                double a2 = 2 * a;
                if (distance < a2) {
                    sb.append(String.valueOf(1));
                    cntinellipse++;
                } else {
                    sb.append(String.valueOf(0));
                }
                for (int i = 0; i < taxa.length; i++) {
                    sb.append("\t").append(l.get(i + 4));
                }
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println(cntinellipse);
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 根据计算好的pv sd depth文件来过滤椭圆,在文件键入2列 CHROM	POS	AverageDepth	SD	Distance	Group
     * PValue	IfPVzero	AFG-L1
     */
    public void calEllipse2() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/005_pvCalSample";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/006_addellipse2";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fsList = Arrays.asList(fs);

        /**
         * Mode: 进行椭圆的计算和判断 椭圆是平面上到两个固定点的距离之和为常数 2a 的点之轨迹。 再进行计算判断；
         * 已知椭圆圆心是C(h,k),半长轴是a,半短轴是b，旋转角度是t，半焦距是c，有a2-b2=c2, 1.先求半焦距c的长度;
         * 2.求焦点1的坐标F1X F1Y h和焦点2的坐标F2X F2Y 3.定义一个新的 点坐标对象；
         *
         * Mode: 进行椭圆的计算和判断 椭圆是平面上到两个固定点的距离之和为常数 2a 的点之轨迹。 再进行计算判断；
         * 已知椭圆圆心是C(h,k),半长轴是a,半短轴是b，旋转角度是t，半焦距是c，有a2-b2=c2,
         * 1.先求半焦距c的长度;2.求焦点1的坐标F1X F1Y h和焦点2的坐标F2X F2Y 3.定义一个新的点坐标对象；
         */
        double a = 2;
        double b = 0.58;
        double h = 9.8;
        double k = 3.6;
        double t = 9;
        //测试用double a = 5 ; double b = 3; double h = 10; double k = 3.7 ; double t = 30 ; 
        double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
        System.out.println("The length c is  " + c);
        double sint = Math.sin((Math.PI * t / 180));
        double cost = Math.cos((Math.PI * t / 180));
        double ccost = c * cost;
        double csint = c * sint;
        System.out.println("c " + c);
        System.out.println("sint  " + sint);
        System.out.println("cost  " + cost);
        System.out.println("csint  " + csint);
        System.out.println("ccost  " + ccost);
        double f1x = h - (c * cost);
        double f1y = k - (c * sint);
        double f2x = h + (c * cost);
        double f2y = k + (c * sint);
        System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
        System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
        //设置两个焦点的坐标p1 p2
        Point p1 = new Point();
        p1.setLocation(f1x, f1y);
        Point p2 = new Point();
        p2.setLocation(f2x, f2y);

        System.out.println("p1 " + p1.getX() + "    " + p1.getY());
        System.out.println("p2 " + p2.getX() + "    " + p2.getY());
        //定位坐标
//            System.out.println("p1的x坐标为"+p1.getX());
//            System.out.println("p1的y坐标为"+p1.getY());
//            System.out.println("p2的x坐标为"+p2.getX());
//            System.out.println("p2的y坐标为"+p2.getY());
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", ".addellipse2.txt")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            int cntinellipse = 0;

            try {
                //先处理表头
                String temp = br.readLine();
                List<String> l = new ArrayList<>();
                l = PStringUtils.fastSplit(temp);
                bw.write(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + "Distance" + "\t" + "Group");
                String[] taxa = new String[l.size() - 4];
                for (int i = 0; i < taxa.length; i++) {
                    bw.write("\t" + l.get(i + 4));
                }
                bw.newLine();
                /*
                椭圆内分组为0，椭圆外分组为1.
                 */
                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp, "\t");
                    double x = Double.parseDouble(l.get(2));
                    double y = Double.parseDouble(l.get(3));
                    //Point p = new Point();
                    //p.setLocation(x, y);
                    // 计算两点间距离公式
                    //double distance1 = Math.sqrt(Math.abs((p1.getX() - p.getX())* (p1.getX() - p.getX())+(p1.getY() - p.getY())* (p1.getY() - p.getY())));
                    //double distance2 = Math.sqrt(Math.abs((p2.getX() - p.getX())* (p2.getX() - p.getX())+(p2.getY() - p.getY())* (p2.getY() - p.getY())));
                    double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                    double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                    double distance = distance1 + distance2;
                    System.out.println("两点间的距离是:" + distance);

                    //开始写出文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + String.format("%.4f", distance) + "\t");
                    double a2 = 2 * a;
                    if (distance < a2) {
                        sb.append(String.valueOf(1));
                        cntinellipse++;
                    } else {
                        sb.append(String.valueOf(0));
                    }
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(l.get(i + 4));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(cntinellipse);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 对已经生成的2万行结果(包含Depth Variance PV)进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCovevsSDvsPV() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/004_pvCal/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/005_pvCalSample/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".depth.pv.txt");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".depth.pv.txt", ".depth.pv_5000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }

    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计和Pvalue获取，制成一个表格，chr pos averageDepth SD PValue
     * IfPVzero taxa1Depth ..... 注意表格不要以#开头，否则被注释，看不到
     * 代码，即在之前depth和SD的基础上又加了2列，第一列是PValue 第二列是PValue是否是0的判断,如果是，那么值为1，如果不是那么值为0
     */
    public void statVcfPValue() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/000_subSet";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/004_pvCal/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".depth.pv.txt")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD\tPValue\tIfPVzero");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    TDoubleArrayList PValueList = new TDoubleArrayList();

                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    String pvalue = l.get(7).split("PV=")[1].split(";")[0];
                    String ifPVzero = null;
                    if (pvalue.equals("0.0")) {
                        ifPVzero = "1";
                    } else {
                        ifPVzero = "0";
                    }
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    //计算完毕，接下来开始写入文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(relativeMean).append("\t").append(sd)
                            .append("\t").append(pvalue).append("\t").append(ifPVzero);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }

    public void calEllipse() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/002_depthCalSample/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/003_addellipse/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        List<File> fsList = Arrays.asList(fs);

        /**
         * Mode: 进行椭圆的计算和判断 椭圆是平面上到两个固定点的距离之和为常数 2a 的点之轨迹。 再进行计算判断；
         * 已知椭圆圆心是C(h,k),半长轴是a,半短轴是b，旋转角度是t，半焦距是c，有a2-b2=c2,
         * 1.先求半焦距c的长度;2.求焦点1的坐标F1X F1Y h和焦点2的坐标F2X F2Y 3.定义一个新的 点坐标对象；
         */
        double a = 2;
        double b = 0.58;
        double h = 10;
        double k = 3.7;
        double t = 8;
        //测试用double a = 5 ; double b = 3; double h = 10; double k = 3.7 ; double t = 30 ; 
        double c = Math.sqrt(Math.pow(a, 2) - Math.pow(b, 2));
        System.out.println("The length c is  " + c);
        double sint = Math.sin((Math.PI * t / 180));
        double cost = Math.cos((Math.PI * t / 180));
        double ccost = c * cost;
        double csint = c * sint;
        System.out.println("c " + c);
        System.out.println("sint  " + sint);
        System.out.println("cost  " + cost);
        System.out.println("csint  " + csint);
        System.out.println("ccost  " + ccost);
        double f1x = h - (c * cost);
        double f1y = k - (c * sint);
        double f2x = h + (c * cost);
        double f2y = k + (c * sint);
        System.out.println("The coordinate f1x is  " + f1x + "\n" + "The coordinate f1y is  " + f1y);
        System.out.println("The coordinate f2x is  " + f2x + "\n" + "The coordinate f2y is  " + f2y);
        //设置两个焦点的坐标p1 p2
        Point p1 = new Point();
        p1.setLocation(f1x, f1y);
        Point p2 = new Point();
        p2.setLocation(f2x, f2y);

        System.out.println("p1 " + p1.getX() + "    " + p1.getY());
        System.out.println("p2 " + p2.getX() + "    " + p2.getY());
        //定位坐标
//            System.out.println("p1的x坐标为"+p1.getX());
//            System.out.println("p1的y坐标为"+p1.getY());
//            System.out.println("p2的x坐标为"+p2.getX());
//            System.out.println("p2的y坐标为"+p2.getY());

        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".txt", ".addellipse.txt")).getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            int cntinellipse = 0;

            try {
                //先处理表头
                String temp = br.readLine();
                List<String> l = new ArrayList<>();
                l = PStringUtils.fastSplit(temp);
                bw.write(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + "Distance" + "\t" + "Group");
                String[] taxa = new String[l.size() - 4];
                for (int i = 0; i < taxa.length; i++) {
                    bw.write("\t" + l.get(i + 4));
                }
                bw.newLine();
                /*
                椭圆内分组为0，椭圆外分组为1.
                 */

                while ((temp = br.readLine()) != null) {
                    l = PStringUtils.fastSplit(temp, "\t");
                    double x = Double.parseDouble(l.get(2));
                    double y = Double.parseDouble(l.get(3));
                    //Point p = new Point();
                    //p.setLocation(x, y);
                    // 计算两点间距离公式
                    //double distance1 = Math.sqrt(Math.abs((p1.getX() - p.getX())* (p1.getX() - p.getX())+(p1.getY() - p.getY())* (p1.getY() - p.getY())));
                    //double distance2 = Math.sqrt(Math.abs((p2.getX() - p.getX())* (p2.getX() - p.getX())+(p2.getY() - p.getY())* (p2.getY() - p.getY())));
                    double distance1 = Math.sqrt(Math.pow((x - f1x), 2) + Math.pow((y - f1y), 2));
                    double distance2 = Math.sqrt(Math.pow((x - f2x), 2) + Math.pow((y - f2y), 2));
                    double distance = distance1 + distance2;
                    System.out.println("两点间的距离是:" + distance);

                    //开始写出文件
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0) + "\t" + l.get(1) + "\t" + l.get(2) + "\t" + l.get(3) + "\t" + String.format("%.4f", distance) + "\t");
                    double a2 = 2 * a;
                    if (distance < a2) {
                        sb.append(String.valueOf(1));
                        cntinellipse++;
                    } else {
                        sb.append(String.valueOf(0));
                    }
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(l.get(i + 4));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
                System.out.println(cntinellipse);
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
    }

    /**
     * 对已经生成的2万行结果进行随机抽样，抽出5000行进行接下来的随机分析
     */
    public void subsetCoveVSSD() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/001_depthCal/";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/002_depthCalSample/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".depth.txt");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".depth.txt", ".depth_5000sites.txt")).getAbsolutePath();
            RowTable<String> t = new RowTable<>(f.getAbsolutePath());
            boolean[] ifOut = new boolean[t.getRowNumber()];
            int totallines = t.getRowNumber();
            double ratio = (double) 5000 / totallines; //注意一定要在5000千加上 强制类型转换，不然不能得出小数
            for (int i = 0; i < t.getRowNumber(); i++) {
                double r = Math.random();
                if (r > ratio) {

                } else {
                    ifOut[i] = true;
                }
            }
            t.writeTextTable(outfileS, IOFileFormat.Text, ifOut);
        });
    }

    /**
     * 对抽样的vcf文件进行每个位点每个taxa的深度统计，制成一个表格，chr pos averageDepth SD taxa1Depth
     * ..... 注意表格不要以#开头，否则被注释，看不到
     */
    public void statVcfCoverage() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/000_subSet";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/001_depthCal/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        List<File> fsList = Arrays.asList(fs);
        fsList.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(f.getAbsolutePath());
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf", ".depth.txt")).getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            int cnt = 0;
            try {
                while ((temp = br.readLine()).startsWith("##")) {
                }
                List<String> l = PStringUtils.fastSplit(temp, "\t");
                bw.write(l.get(0).replaceFirst("#", "") + "\t" + l.get(1) + "\t" + "AverageDepth\tSD");
                String[] taxa = new String[l.size() - 9];
                for (int i = 0; i < taxa.length; i++) {
                    taxa[i] = l.get(i + 9);
                    bw.write("\t" + taxa[i]);
                }
                bw.newLine();
                //建立一个整型list类型的数组，每个元素是一个list,一共有 taxa.length个list
                cnt = 0;
                while ((temp = br.readLine()) != null) {
                    TDoubleArrayList depthList = new TDoubleArrayList();
                    cnt++; // 对snp开始计数
                    if (cnt % 1000 == 0) {
                        System.out.println(String.valueOf(cnt) + " lines");
                    }
                    l = PStringUtils.fastSplit(temp, "\t");
                    for (int i = 0; i < taxa.length; i++) {

                        String genoS = l.get(i + 9);
                        if (genoS.startsWith(".")) {
                            depthList.add(0);
                            continue;
                        }
                        List<String> ll = PStringUtils.fastSplit(genoS, ":");
                        List<String> lll = PStringUtils.fastSplit(ll.get(1), ",");
                        int depth = Integer.valueOf(lll.get(0)) + Integer.valueOf(lll.get(1));
                        depthList.add(depth);
                    }
                    double[] dep = depthList.toArray();
                    DescriptiveStatistics d = new DescriptiveStatistics(dep);
                    double relativeMean = d.getMean();
                    double sd = d.getStandardDeviation();
                    StringBuilder sb = new StringBuilder();
                    sb.append(l.get(0)).append("\t").append(l.get(1)).append("\t").append(relativeMean).append("\t").append(sd);
                    for (int i = 0; i < taxa.length; i++) {
                        sb.append("\t").append(String.format("%.0f", depthList.get(i)));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
            System.out.println(f.getName() + " is calculated well done");
        });
    }

    /**
     * 对chr1A 1B 1D进行随机抽样，抽取2万个sites进行测试，过滤掉有2个alt的位点。 30651417 0.001786084
     * 20000	11197684 chr1A.vcf.gz 0.00136431	20000	14659430 chr1B.vcf.gz
     * 0.004177731	20000	4787288 chr1D.vcf.gz
     */
    public void subsetVCF() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/originData";
        String outfileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/001_filterVCF/000_subSet";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");

        try {
            long startTime = System.nanoTime();
            //long endTime = System.nanoTime();
            //float excTime = (float) (endTime - startTime) / 1000000000;
            //System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
            String chr = "1A";
            String infileS = new File(infileDirS, "chr1A.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, "chr1A_subset.vcf").getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            double ratioA = 0.001786084;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                double r = Math.random();
                if (r > ratioA) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(4).contains(",")) {
                    continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished chr1A sampling");
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            long startTime = System.nanoTime();
            String chr = "1B";
            String infileS = new File(infileDirS, "chr1B.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, "chr1B_subset.vcf").getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            double ratioB = 0.00136431;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                double r = Math.random();
                if (r > ratioB) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(4).contains(",")) {
                    continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished chr1B sampling");
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        try {
            long startTime = System.nanoTime();
            String chr = "1D";
            String infileS = new File(infileDirS, "chr1D.vcf.gz").getAbsolutePath();
            String outfileS = new File(outfileDirS, "chr1D_subset.vcf").getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            double ratioB = 0.004177731;
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                }
                double r = Math.random();
                if (r > ratioB) {
                    continue;
                }
                List<String> l = PStringUtils.fastSplit(temp);
                if (l.get(4).contains(",")) {
                    continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("Finished chr1D sampling");
            long endTime = System.nanoTime();
            float excTime = (float) (endTime - startTime) / 1000000000;
            System.out.println("Execution time: " + String.format("%.2f", excTime) + "s");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

        /**
         * Here is the entrance of wheatGload ! Finished chr1A sampling
         * Execution time: 65.63s Finished chr1B sampling Execution time: 81.63s
         * Finished chr1D sampling Execution time: 27.61s
         */
    }

    /**
     * 对合并后的vcf进行变异数数统计
     */
    public void countVCF() {
        String infileDirS = "/Users/Aoyue/project/wheatVMapII/005_preTest/originData";
        File[] fs = new File(infileDirS).listFiles();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                System.out.println(fs[i].getName() + " is hidden");
                fs[i].delete();
            }
        } //注意有隐藏文件，需要进行删除后重新列出文件目录。
        fs = new File(infileDirS).listFiles();
        int sum = 0;
        for (int i = 0; i < fs.length; i++) {
            int cnt = 0;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(fs[i].getAbsolutePath());
                String temp = null;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    cnt++;
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println(String.valueOf(cnt) + "\t" + fs[i].getName());
            sum += cnt;
        }
        System.out.println(String.valueOf(sum));

        /**
         * 14659430	chr1B.vcf.gz 7015	chr1D.vcf.gz.tbi 4787288	chr1D.vcf.gz
         * 11197684	chr1A.vcf.gz 30651417
         */
    }

    /**
     * 对chr 1 2 3 4 5 6 进行合并，组成1A 1B 1D染色体 1	0	471304005	chr1A	0	471304005 2	0
     * 122798051	chr1A	471304005	594102056 3	0	438720154	chr1B	0	438720154 4	0
     * 251131716	chr1B	438720154	689851870 5	0	452179604	chr1D	0	452179604 6	0
     * 43273582	chr1D	452179604	495453186 举例：chr1A 0-10(bed) 实际1based 1-10 chr1A
     * 10-15(bed) 实际1based 11-15
     */
    public void mergewheatchrVCF() {
        String infileDirS = "/data2/aoyue/fastCall_project/vcf";
        String outfileDirS = "/data2/aoyue/fastCall_project/mergedChrABDvcf";

        //String infileDirS = "/Users/Aoyue/Documents/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/000_sampleVCF/";
        //String outfileDirS = "/Users/Aoyue/Documents/fastCall/004_fastV2_JiaoDataParameters/003_testvcf/001_testMergedVCF/";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".VCF.txt");
        int A1 = 471304005;
        int B1 = 438720154;
        int D1 = 452179604;
        try {
            for (int i = 0; i < fs.length; i++) {
                if (i == 0) {
                    System.out.println("Now we are starting to merge 1A chomosome");
                    String chrABD = "1A";
                    String chr = PStringUtils.getNDigitNumber(3, i + 1);
                    String chrSecond = PStringUtils.getNDigitNumber(3, i + 2);
//                    String infileS = new File(infileDirS,"chr"+chr+"_5000.VCF.txt").getAbsolutePath();
//                    String infile2S = new File(infileDirS,"chr"+chrSecond+"_5000.VCF.txt").getAbsolutePath();
//                    String outfileS = new File(outfileDirS,"chr"+chrABD +"_test.vcf").getAbsolutePath();

                    String infileS = new File(infileDirS, "chr" + chr + ".VCF.txt").getAbsolutePath();
                    String infile2S = new File(infileDirS, "chr" + chrSecond + ".VCF.txt").getAbsolutePath();
                    String outfileS = new File(outfileDirS, "chr" + chrABD + ".vcf").getAbsolutePath();

                    BufferedReader br = IOUtils.getTextReader(infileS);
                    BufferedReader br2 = IOUtils.getTextReader(infile2S);
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = new ArrayList();
                            l = PStringUtils.fastSplit(temp);
                            bw.write(chrABD);
                            StringBuilder sb = new StringBuilder();
                            for (int j = 1; j < l.size(); j++) {
                                sb.append("\t");
                                sb.append(l.get(j));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                    br.close();
                    /*
                    开始进行chr1和chr2染色体的合并
                     */
                    while ((temp = br2.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            continue;
                        }
                        List<String> l = new ArrayList();
                        l = PStringUtils.fastSplit(temp);
                        int pos = Integer.parseInt(l.get(1)) + A1;
                        bw.write(chrABD);
                        bw.write("\t");
                        bw.write(String.valueOf(pos));
                        StringBuilder sb = new StringBuilder();
                        for (int j = 2; j < l.size(); j++) {
                            sb.append("\t");
                            sb.append(l.get(j));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    br2.close();
                    bw.flush();
                    bw.close();
                    System.out.println("Merging 1A chomosome is finished");
                }
                if (i == 2) {
                    System.out.println("Now we are starting to merge 1B chomosome"); //need revise
                    String chrABD = "1B"; //need revise
                    String chr = PStringUtils.getNDigitNumber(3, i + 1);
                    String chrSecond = PStringUtils.getNDigitNumber(3, i + 2);
                    String infileS = new File(infileDirS, "chr" + chr + ".VCF.txt").getAbsolutePath();
                    String infile2S = new File(infileDirS, "chr" + chrSecond + ".VCF.txt").getAbsolutePath();
                    String outfileS = new File(outfileDirS, "chr" + chrABD + ".vcf").getAbsolutePath();

                    BufferedReader br = IOUtils.getTextReader(infileS);
                    BufferedReader br2 = IOUtils.getTextReader(infile2S);
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = new ArrayList();
                            l = PStringUtils.fastSplit(temp);
                            bw.write(chrABD);
                            StringBuilder sb = new StringBuilder();
                            for (int j = 1; j < l.size(); j++) {
                                sb.append("\t");
                                sb.append(l.get(j));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                    br.close();
                    /**
                     * 开始进行chr3和chr4染色体的合并 //need revise
                     */
                    while ((temp = br2.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            continue;
                        }
                        List<String> l = new ArrayList();
                        l = PStringUtils.fastSplit(temp);
                        int pos = Integer.parseInt(l.get(1)) + B1; //need revise
                        bw.write(chrABD);
                        bw.write("\t");
                        bw.write(String.valueOf(pos));
                        StringBuilder sb = new StringBuilder();
                        for (int j = 2; j < l.size(); j++) {
                            sb.append("\t");
                            sb.append(l.get(j));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    br2.close();
                    bw.flush();
                    bw.close();
                    System.out.println("Merging 1B chomosome is finished"); //need revise
                }
                if (i == 4) {
                    System.out.println("Now we are starting to merge 1D chomosome"); //need revise
                    String chrABD = "1D"; //need revise
                    String chr = PStringUtils.getNDigitNumber(3, i + 1);
                    String chrSecond = PStringUtils.getNDigitNumber(3, i + 2);
                    String infileS = new File(infileDirS, "chr" + chr + ".VCF.txt").getAbsolutePath();
                    String infile2S = new File(infileDirS, "chr" + chrSecond + ".VCF.txt").getAbsolutePath();
                    String outfileS = new File(outfileDirS, "chr" + chrABD + ".vcf").getAbsolutePath();

                    BufferedReader br = IOUtils.getTextReader(infileS);
                    BufferedReader br2 = IOUtils.getTextReader(infile2S);
                    BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                    String temp = null;
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            bw.write(temp);
                            bw.newLine();
                        } else {
                            List<String> l = new ArrayList();
                            l = PStringUtils.fastSplit(temp);
                            bw.write(chrABD);
                            StringBuilder sb = new StringBuilder();
                            for (int j = 1; j < l.size(); j++) {
                                sb.append("\t");
                                sb.append(l.get(j));
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                    br.close();
                    /*
                    开始进行chr5和chr6染色体的合并 //need revise
                     */
                    while ((temp = br2.readLine()) != null) {
                        if (temp.startsWith("#")) {
                            continue;
                        }
                        List<String> l = new ArrayList();
                        l = PStringUtils.fastSplit(temp);
                        int pos = Integer.parseInt(l.get(1)) + D1; //need revise
                        bw.write(chrABD);
                        bw.write("\t");
                        bw.write(String.valueOf(pos));
                        StringBuilder sb = new StringBuilder();
                        for (int j = 2; j < l.size(); j++) {
                            sb.append("\t");
                            sb.append(l.get(j));
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    br2.close();
                    bw.flush();
                    bw.close();
                    System.out.println("Merging 1D chomosome is finished"); //need revise
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立42条染色体的java 运行脚本。
     */
    public void mkJavaCmdchr1_42() {
        String outfileS = "/Users/Aoyue/Documents/fastCall/004_fastV2_JiaoDataParameters/002_javaCMD/fastCall_chr1_42.sh";
        //nohup java -Xmx200g -jar fastCall_V2.jar parameters_001_FastCallV2_Jiao.txt > log_001_fastcall.txt &
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < 45; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                bw.write("java -Xmx200g -jar fastCall_V2.jar parameters_");
                bw.write(chr);
                bw.write("_FastCallV2_Jiao.txt > log_");
                bw.write(chr);
                bw.write("_fastcall.txt");
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }

    }

    /**
     * 本方法的目的是：建立42条染色体的parameters文件。
     */
    public void mkParameterchr1_42() {
        String outfileDirS = "/Users/Aoyue/Documents/fastCall/004_fastV2_JiaoDataParameters/001_parameters";
        try {
            for (int i = 0; i < 45; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "parameters_" + chr + "_FastCallV2_Jiao.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                StringBuilder sb = new StringBuilder();
                sb.append("Author: Fei Lu\n");
                sb.append("Email: flu@genetics.ac.cn; dr.lufei@gmail.com\n");
                sb.append("Homepage: https://plantgeneticslab.weebly.com/\n").append("\n");
                sb.append("#This WGS SNP discovery pipeline are designed for both heterozygous and inbred species, especially the depth is high (e.g. >10X each taxon).\n");
                sb.append("#To run and pipeline, the machine should have Java 8 and samtools installed. The lib directory should stay with FastCall.jar in the same folder. Command (e.g.): java -Xmx200g -jar FastCall.jar parameter.txt > log.txt &\n");
                sb.append("#To specify reference, bam direcotory, chromosome, and output direcotory, please edit the the parameters below. Also, please keep the order of parameters.\n").append("\n").append("\n");
                sb.append("#Parameter1:\tReference genome file with an index file (.fai). The reference should be in FastA format. Chromosomes are labled as 1-based numbers (1,2,3,4,5...).\n");
                sb.append("/data1/publicData/wheat/reference/v1.0/ABD/abd_iwgscV1.fa").append("\n").append("\n");
                sb.append("#Parameter2:\tBam directory, where your bam files are\n");
                sb.append("/data2/sharedData/wheat_Jiao/splitBamfile/").append(chr).append("/").append("\n").append("\n");
                sb.append("#Parameter3:\tTaxa bam information file, including the info about what bams are included for each taxon\n");
                sb.append("/data2/aoyue/fastCall_project/taxaBamMap/taxaBamMap_").append(chr).append("_Jiao.txt").append("\n").append("\n");
                sb.append("#Parameter4:\tChromosome on which genotyping will be performed (e.g. 1)\n");
                sb.append(i).append("\n").append("\n");
                sb.append("#Parameter5:\tVCF output directory\n");
                sb.append("/data2/aoyue/fastCall_project/vcf/").append("\n");
                bw.write(sb.toString());
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立42条染色体的taxa-bam map图。
     */
    public void mkFastCallTaxaBamMapchr1_42() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/Jiao/002_script/wheat119SM_Jiao.t.txt";
        String outfileDirS = "/Users/Aoyue/Documents/fastCall/004_fastV2_JiaoDataParameters/000_taxaBamMap/";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        String[] name = namelist.toArray(new String[namelist.size()]);
        Arrays.sort(name);
        try {
            for (int i = 0; i < 45; i++) {
                String chr = PStringUtils.getNDigitNumber(3, i);
                String outfileS = new File(outfileDirS, "taxaBamMap_" + chr + "_Jiao.txt").getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(outfileS);
                for (int j = 0; j < name.length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(name[j]).append("\t").append(name[j]).append(".chr").append(chr).append(".bam");
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }

        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }

    /**
     * 本方法的目的是：建立一张taxa-bam map,1号染色体，先做一遍测试，成功后则建立42条染色体的map图
     */
    public void mkFastCallTaxaBamMap() {
        String infileS = "/Users/Aoyue/project/wheatVMapII/Jiao/002_script/wheat119SM_Jiao.t.txt";
        String outfileS = "/Users/Aoyue/Documents/taxaBamMap_001_Jiao.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> namelist = t.getColumn(0);
        String[] name = namelist.toArray(new String[namelist.size()]);
        Arrays.sort(name);
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < name.length; i++) {
                StringBuilder sb = new StringBuilder();
                sb.append(name[i]).append("\t").append(name[i]).append(".chr001.bam");
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
}
