package xiaohan.GBS;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.koloboke.collect.map.IntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import daxing.common.LibraryOfGRT;
import pgl.AppUtils;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import xiaohan.utils.RowTable;
import xiaohan.utils.SNPmappingInGene;
import xiaohan.utils.IOUtils;

import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.atomic.LongAdder;

import static cern.jet.math.Arithmetic.factorial;

public class heterozygosity {
    String barcodeFileS = null;
    String libraryFastqMapFileS = null;
    String cutter1 = null;
    String cutter2 = null;
    String kinship = "/data1/home/xiaohan/nam/KinshipMap.txt";
    String workingDirS = null;
    String[] subDir = {"subFastqs", "bam", "sortedbam", "IBS", "heter", "sum"};

    String genotypeDir = null;
    String BamDir = null;
    String ref = null;
    String plate = null;
    String bamsuffix = ".sorted.bam";
    String threads = null;
    String HapscannerDir = null;

    String bwa = null;
    String bwalib = null;
    String samtools = null;

    Multimap<String, String> kinshipMap = ArrayListMultimap.create();
    HashSet<String> HBSet = new HashSet<>();

    String parameterDir = null;
    String taxaRefBAMDir = null;
    String posDir = null;
    String posAlleleDir = null;

    String taxaRefBamFileS = null;
    //The posAllele file (with header), the format is Chr\tPos\tRef\tAlt (from VCF format). The positions come from haplotype library.
    String posAlleleFileS = null;
    //The pos files (without header), the format is Chr\tPos. The positions come from haplotype library.
    String posFileS = null;
    //The chromosome which will be genotyped
    int chr = -1;
    //The path of samtools
    String samtoolsPath = null;
    //The directory of output
    String outputDirS = null;
    int regionStart = Integer.MIN_VALUE;
    int regionEnd = Integer.MIN_VALUE;

    int nThreads = -1;

    HashMap<String, List<String>> taxaBamsMap = new HashMap<>();

    HashMap<String, String> taxaRefMap = new HashMap<>();

    String[] subDirS = {"mpileup", "indiVCF", "VCF"};

    IntDoubleMap factorialMap = null;
    int maxFactorial = 150;
    //combined: sequencing error and alignment error
    double combinedErrorRate = 0.05;

    HashMap<Integer, Integer> chrMap = new HashMap<>();

    public heterozygosity(String[] args) {
        this.parseparameters(args[0]);
        this.mkdir();
//        this.parsefastqs();
//        this.bwa(args);
//        this.sortbam(args);
//        this.parameter();
//        this.taxaRefBAM();
//        this.callsnp();
        this.getkinship();
//        this.getheterozygosity();
//        this.mergeFile();
        this.getIBS();
//        this.getsum();
//        this.test1(args);
//        this.test(args);

    }

    public void test1(String... args) {
        String infile1 = args[0];
        String infile2 = args[1];
        String outfile = args[2];
        GenotypeGrid g1 = new GenotypeGrid(infile1, GenoIOFormat.VCF_GZ);
        GenotypeGrid g2 = new GenotypeGrid(infile2, GenoIOFormat.VCF_GZ);
        GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
        SumTaxaDivergence std = new SumTaxaDivergence(g);
        std.writeDxyMatrix(outfile, IOFileFormat.Text);
        g.getIBSDistanceMatrix();
    }

    public void test(String[] args) {
//        LibraryOfGRT li = new LibraryOfGRT(args[0],args[1],args[2],args[3]);
//        li.splitBarcode(args[4]);
//        int a = chrUtils.getChrIndextoChrABDIndex(17);
//        for (int i = 1; i < 22; i++) {
//            System.out.println(chrUtils.getNamechrABD(i));
//        }
//        System.out.print(a);
        String infile = kinship;
        BufferedReader br = IOUtils.getTextReader(infile);
        try {
            String temp = null;
            String[] temps = null;
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temp.startsWith("F2")) continue;
                HBSet.add(temps[0]);
                kinshipMap.put(temps[0], temps[1]);
                kinshipMap.put(temps[0], temps[2]);
            }
            br.close();

            RowTable<String> rt = new RowTable<>(args[2]);
//        System.out.println(rt.getRowNumber());
//        System.out.println(rt.getColumnNumber());
//        System.out.println(rt.getHeader().toArray().length);

            for (int i = 0; i < rt.getRowNumber(); i++) {
                for (int j = 0; j < rt.getColumnNumber(); j++) {
                    if (rt.getCell(i, j).startsWith("N")) {
                        rt.setCell(i, j, "1");
                    }
                }
            }

            List<String> HBset = new ArrayList<>();
            List<String> Pset = new ArrayList<>();
            String[] header = rt.getHeader().toArray(new String[0]);
            for (int i = 0; i < header.length; i++) {
                if (header[i].startsWith("HB")) {
                    HBset.add(header[i]);
                } else if (!header[i].startsWith("Dxy")) {
                    Pset.add(header[i]);
                }
            }


            System.out.println(HBset.size());
            System.out.println(Pset.size());

            BufferedWriter bw = IOUtils.getTextWriter(args[3]);
            for (int i = 0; i < HBset.size(); i++) {
                String HB = HBset.get(i);
                System.out.println(HB);
                Collection<String> P = kinshipMap.get(HB);
//            String[] pa = rt.getColumn(rt.getColumnIndex("Dxy")).toArray(new String[0]);
//            System.out.println(pa[0]);
                double[] IBS = rt.getColumnAsDoubleArray(rt.getColumnIndex(HB));
                double[] IBSsub = Arrays.copyOfRange(IBS, HBset.size(), header.length - 1);
                double[] minMin = SNPmappingInGene.minMin(IBSsub);
//            System.out.println(HBset.size());
//            System.out.println(Pset.size());
                System.out.println(IBSsub.length);
//            System.out.println(header.length);
                String p1 = header[HBset.size() + (int) minMin[2] + 1];
                String p2 = header[HBset.size() + (int) minMin[3] + 1];
                String P1 = P.toArray(new String[0])[0];
                String P2 = P.toArray(new String[0])[1];
                int P1index = Pset.indexOf(P1);
                int P2index = Pset.indexOf(P2);
                int[] rank = getIncreaseRanksArray(IBSsub);

//                System.out.println(rank[P1index]);
//                System.out.println(rank[P2index]);
//                System.out.println(rank[Pset.indexOf("Z19-H28")]);
//                System.out.println(rank[P1index]);
//                System.out.println(rank[P2index]);
                System.out.println(HB + " of " + p1 + " & " + p2);
                if (P.contains(p1) && P.contains(p2)) {
                    bw.write(HB + "\tcorrect\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\t" + rank[P1index] + "\t" + rank[P2index] + "\n");
                    continue;
                } else {
                    bw.write(HB + "\twrong\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\t" + rank[P1index] + "\t" + rank[P2index] + "\n");
                    continue;
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static int[] getRanksArray(int[] array) {
        int[] result = new int[array.length];

        for (int i = 0; i < array.length; i++) {
            int count = 0;
            for (int j = 0; j < array.length; j++) {
                if (array[j] > array[i]) {
                    count++;
                }
            }
            result[i] = count + 1;
        }
        return result;
    }

    public static int[] getDecreaseRanksArray(double[] array) {
        int[] result = new int[array.length];

        for (int i = 0; i < array.length; i++) {
            int count = 0;
            for (int j = 0; j < array.length; j++) {
                if (array[j] > array[i]) {
                    count++;
                }
            }
            result[i] = count + 1;
        }
        return result;
    }

    public static int[] getIncreaseRanksArray(double[] array) {
        int[] result = new int[array.length];

        for (int i = 0; i < array.length; i++) {
            int count = 0;
            for (int j = 0; j < array.length; j++) {
                if (array[j] < array[i]) {
                    count++;
                }
            }
            result[i] = count + 1;
        }
        return result;
    }

    public void parseparameters(String infileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(infileS);
        List<String> pLineList = d.getFirstElement();
        workingDirS = pLineList.get(0);
        barcodeFileS = pLineList.get(1);
        libraryFastqMapFileS = pLineList.get(2);
        kinship = pLineList.get(3);
        cutter1 = pLineList.get(4);
        cutter2 = pLineList.get(5);
        plate = pLineList.get(6);

        genotypeDir = pLineList.get(7);
        ref = pLineList.get(8);
        threads = pLineList.get(9);
        HapscannerDir = pLineList.get(10);

        bwa = pLineList.get(11);
        bwalib = pLineList.get(12);
        samtools = pLineList.get(13);
    }

    public void mkdir() {
        for (int i = 0; i < subDir.length; i++) {
            File f = new File(workingDirS, subDir[i]);
            if (!f.exists()) {
                f.mkdir();
            }
        }
    }

    private void getheterozygosity() {
        String inputDir = new File(HapscannerDir, "output/" + plate + "/VCF").getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        for (int i = 0; i < 21; i++) {
            chrSet.add(i);
        }
        StringBuilder sb = new StringBuilder();
        chrSet.parallelStream().forEach(m -> {
            try {
                int[][] count = new int[HBSet.size()][3];
                for (int i = 0; i < HBSet.size(); i++) {
                    for (int j = 0; j < 3; j++) {
                        count[i][j] = 0;
                    }
                }
                for (int k = 1; k <= 2; k++) {
                    int chr = m * 2 + k;
                    int dep1 = 0;
                    int dep2 = 0;
                    int chrABDindex = xiaohan.utils.chrUtils.getChrIndextoChrABDIndex(chr);
                    System.out.println("This is reading " + chr + " which is in " + chrABDindex);
                    String vcf = new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf.gz").getAbsolutePath();
                    String genovcf = new File(genotypeDir, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf.gz").getAbsolutePath();
                    BufferedReader br = IOUtils.getTextGzipReader(vcf);
                    BufferedReader br1 = IOUtils.getTextGzipReader(genovcf);
                    String temp = null;
                    String temp1 = null;
                    temp = br.readLine();
                    temp1 = br1.readLine();
                    while ((!temp.startsWith("#C"))) {
                        if (temp.startsWith("#C")) break;
                        temp = br.readLine();
                    }
                    while ((!temp1.startsWith("#C"))) {
                        if (temp1.startsWith("#C")) break;
                        temp1 = br1.readLine();
                    }
                    HashMap<Integer, String> HBMap = new HashMap<>();
                    HashMap<String, Integer> PMap = new HashMap<>();
                    String[] temps = temp.split("\t");
                    String[] temps1 = temp1.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        HBMap.put(i, temps[i]);
                    }
                    for (int i = 9; i < temps1.length; i++) {
                        PMap.put(temps1[i], i);
                    }
                    String[] HBnames = Arrays.copyOfRange(temps, 9, temps.length);
                    while ((temp = br.readLine()) != null) {
                        temp1 = br1.readLine();
                        temps = temp.split("\t");
                        temps1 = temp1.split("\t");
                        for (int i = 9; i < temps.length; i++) {
                            int het = -1;
//                        System.out.println(temps[i]);
                            String[] P = kinshipMap.get(HBMap.get(i)).toArray(new String[0]);
                            if (temps1[PMap.get(P[0])].startsWith("./.") || temps1[PMap.get(P[1])].startsWith("./.")) {
                                continue;
                            }
                            if (temps1[PMap.get(P[0])].startsWith("0/1") || temps1[PMap.get(P[1])].startsWith("0/1")) {
                                continue;
                            }
                            if (temps1[PMap.get(P[0])].split(":")[0].equals(temps1[PMap.get(P[1])].split(":")[0])) {
                                continue;
                            } else if (!temps1[PMap.get(P[0])].split(":")[0].equals(temps1[PMap.get(P[1])].split(":")[0])) {
                                if (temps[i].startsWith("./.")) {
                                    count[i - 9][2] += 1;
                                    het = -1;
                                } else {
                                    dep1 = Integer.parseInt(temps[i].split(":")[1].split(",")[0]);
                                    dep2 = Integer.parseInt(temps[i].split(":")[1].split(",")[1]);
                                    if ((dep1 + dep2) >= 2) {
                                        if (temps[i].startsWith("0/1")) {
                                            count[i - 9][0] += 1;
                                            het = 1;
                                        } else {
                                            count[i - 9][1] += 1;
                                            het = 0;
                                        }
                                        continue;
                                    } else continue;
                                }
                            }
                        }
                    }
                    br.close();
                    br1.close();
                }
                String chrABDName = xiaohan.utils.chrUtils.getNamechrABD(m + 1);
                BufferedWriter bw = IOUtils.getTextWriter(new File(this.workingDirS, subDir[4] + "/" + chrABDName + ".txt").getAbsolutePath());
                bw.write("Taxa\theterSite\tHomoSite\tMissingSite\n");
                String[] names = HBSet.toArray(new String[0]);
                for (int j = 0; j < HBSet.size(); j++) {
                    bw.write(names[j] + "\t" + count[j][0] + "\t" + count[j][1] + "\t" + count[j][2]);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    private void getsum() {
        try {
            for (int m = 0; m < 21; m++) {
                int chr1 = m * 2 + 1;
                int chr2 = m * 2 + 2;
                String chrABDName = xiaohan.utils.chrUtils.getNamechrABD(m + 1);
                String ibsOutfileS = new File(workingDirS, subDir[3] + "/check_" + chrABDName + ".txt").getAbsolutePath();
                RowTable<String> rt = new RowTable<>(ibsOutfileS);
                for (int i = 0; i < rt.getRowNumber(); i++) {
                    for (int j = 0; j < rt.getColumnNumber(); j++) {
                        if (rt.getCell(i, j).startsWith("N")) {
                            rt.setCell(i, j, "1");
                        }
                    }
                }

                List<String> HBset = new ArrayList<>();
                List<String> Pset = new ArrayList<>();
                String[] header = rt.getHeader().toArray(new String[0]);
                for (int i = 0; i < header.length; i++) {
                    if (header[i].startsWith("HB")) {
                        HBset.add(header[i]);
                    } else if (!header[i].startsWith("Dxy")) {
                        Pset.add(header[i]);
                    }
                }

                BufferedWriter bw = IOUtils.getTextWriter(new File(workingDirS, subDir[subDir.length - 1] + "/sum" + chrABDName + ".txt").getAbsolutePath());
                for (int i = 0; i < HBset.size(); i++) {
                    String HB = HBset.get(i);
                    Collection<String> P = kinshipMap.get(HB);
                    double[] IBS = rt.getColumnAsDoubleArray(rt.getColumnIndex(HB));
                    double[] IBSsub = Arrays.copyOfRange(IBS, HBset.size(), header.length - 1);
                    double[] minMin = SNPmappingInGene.minMin(IBSsub);
                    String p1 = header[HBset.size() + (int) minMin[2] + 1];
                    String p2 = header[HBset.size() + (int) minMin[3] + 1];
                    String P1 = P.toArray(new String[0])[0];
                    String P2 = P.toArray(new String[0])[1];
                    int P1index = Pset.indexOf(P1);
                    int P2index = Pset.indexOf(P2);
                    int[] rank = getIncreaseRanksArray(IBSsub);
                    System.out.println(rank[P1index]);
                    System.out.println(rank[P2index]);
                    System.out.println(HB + " of " + p1 + " & " + p2);
                    if (P.contains(p1) && P.contains(p2)) {
                        bw.write(HB + "\tcorrect\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\t" + rank[P1index] + "\t" + rank[P2index] + "\n");
                        continue;
                    } else {
                        bw.write(HB + "\twrong\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\t" + rank[P1index] + "\t" + rank[P2index] + "\n");
                        continue;
                    }
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void mergeFile(){
        String inputDir = new File(HapscannerDir, "output/" + plate + "/VCF").getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        try {
            for (int m = 0; m < 21; m++) {
                int chr1 = m * 2 + 1;
                int chr2 = m * 2 + 2;
                String chrABDName = xiaohan.utils.chrUtils.getNamechrABD(m + 1);
                sb.setLength(0);
                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr1) + ".vcf").getAbsolutePath() + "\n");
                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr2) + ".vcf").getAbsolutePath() + "\n");
                sb.append("vcf-concat " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr1) + ".vcf.gz").getAbsolutePath());
                sb.append(" " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr2) + ".vcf.gz").getAbsolutePath());
                sb.append(" > " + new File(inputDir,  chrABDName+".vcf").getAbsolutePath() + "\n");
                sb.append("bgzip " + new File(inputDir, chrABDName+".vcf").getAbsolutePath() + "\n");
//                sb.append("rm " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr1) + ".vcf.gz").getAbsolutePath() + "\n");
//                sb.append("rm " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr2) + ".vcf.gz").getAbsolutePath() + "\n");
                System.out.println(sb.toString());
                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", sb.toString()};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private void getIBS() {
        String inputDir = new File(HapscannerDir, "output/" + plate + "/VCF").getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        try {
            for (int m = 0; m < 21; m++) {
                int chr1 = m * 2 + 1;
                int chr2 = m * 2 + 2;
                String chrABDName = xiaohan.utils.chrUtils.getNamechrABD(m + 1);
//                sb.setLength(0);
//                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr1) + ".vcf").getAbsolutePath() + "\n");
//                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr2) + ".vcf").getAbsolutePath() + "\n");
//                System.out.println(sb.toString());
//                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
//                String[] cmdarry = {"/bin/bash", "-c", sb.toString()};
//                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                p.waitFor();

                String vcf1 = new File(inputDir, chrABDName + ".vcf.gz").getAbsolutePath();
                String vcf2 = new File(genotypeDir, chrABDName + ".vcf.gz").getAbsolutePath();
//
                GenotypeGrid g1 = new GenotypeGrid(vcf1, GenoIOFormat.VCF_GZ);
                GenotypeGrid g2 = new GenotypeGrid(vcf2, GenoIOFormat.VCF_GZ);

                String ibsOutfileS = new File(workingDirS, subDir[3] + "/check_" + chrABDName + ".txt").getAbsolutePath();
//
                GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
                SumTaxaDivergence std = new SumTaxaDivergence(g);
                std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
                g.getIBSDistanceMatrix();
                RowTable<String> rt = new RowTable<>(ibsOutfileS);
                for (int i = 0; i < rt.getRowNumber(); i++) {
                    for (int j = 0; j < rt.getColumnNumber(); j++) {
                        if (rt.getCell(i, j).startsWith("N")) {
                            rt.setCell(i, j, "1");
                        }
                    }
                }

                List<String> HBset = new ArrayList<>();
                List<String> Pset = new ArrayList<>();
                String[] header = rt.getHeader().toArray(new String[0]);
                for (int i = 0; i < header.length; i++) {
                    if (header[i].startsWith("HB")) {
                        HBset.add(header[i]);
                    } else if (!header[i].startsWith("Dxy")) {
                        Pset.add(header[i]);
                    }
                }

                BufferedWriter bw = IOUtils.getTextWriter(new File(workingDirS, subDir[subDir.length - 1] + "/sum" + chrABDName + ".txt").getAbsolutePath());
                for (int i = 0; i < HBset.size(); i++) {
                    String HB = HBset.get(i);
                    String HBnickname = HB.substring(0,5);
                    Collection<String> P = kinshipMap.get(HBnickname);
                    double[] IBS = rt.getColumnAsDoubleArray(rt.getColumnIndex(HB));
                    System.out.println(IBS.length);
                    System.out.println(HBset.size());
                    double[] IBSsub = Arrays.copyOfRange(IBS, HBset.size(), header.length - 1);
                    System.out.println(IBSsub.length);
                    double[] minMin = SNPmappingInGene.minMin(IBSsub);
                    System.out.println(minMin.length);
                    System.out.println(HBset.size());
                    System.out.println(Pset.size());
                    System.out.println(IBSsub.length);
                    System.out.println(header.length);
                    String p1 = header[HBset.size() + (int) minMin[2] + 1];
                    String p2 = header[HBset.size() + (int) minMin[3] + 1];
                    String P1 = P.toArray(new String[0])[0];
                    String P2 = P.toArray(new String[0])[1];
                    int P1index = Pset.indexOf(P1);
                    int P2index = Pset.indexOf(P2);
                    System.out.println(HB + " of " + p1 + " & " + p2);
                    if (P.contains(p1) && P.contains(p2)) {
                        bw.write(HB + "\tcorrect\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\n");
                        continue;
                    } else {
                        bw.write(HB + "\twrong\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\t" + P1 + "\t" + P2 + "\t" + IBSsub[P1index] + "\t" + IBSsub[P2index] + "\n");
                        continue;
                    }
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void getkinship() {
        String infile = kinship;
        BufferedReader br = IOUtils.getTextReader(infile);
        try {
            String temp = null;
            String[] temps = null;
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                HBSet.add(temps[0]);
                kinshipMap.put(temps[0], temps[1]);
                kinshipMap.put(temps[0], temps[2]);
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void callsnp() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is calling snp *************************************************************************");
        File[] fs = new File(this.parameterDir).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, plate + "_");
        Arrays.sort(fs);
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i].getAbsolutePath());
            this.parseParameters(fs[i].getAbsolutePath());

            this.mkDir();
            this.scanIndiVCFByThreadPool();
            this.mkFinalVCF();
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Calling snp takes " + (endTime - startTime) + "ms");
    }

    public void poswithAllele() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing pos file ****************************************************");
        HashSet<String> nameSet = new HashSet<>();
        File[] fs = new File(genotypeDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        for (int i = 0; i < fs.length; i++) {
            System.out.println(fs[i].getAbsolutePath());
            nameSet.add(fs[i].getName());
        }
        nameSet.parallelStream().forEach(p -> {
            System.out.println("Reading file :" + p);
            int chr = Integer.parseInt(p.split("\\.")[0]);
            String infile = new File(genotypeDir, p).getAbsolutePath();
            String outfile = new File(posDir, "pos_chr" + chr + ".txt").getAbsolutePath();
            String outfile1 = new File(posAlleleDir, "posAllele_chr" + chr + ".txt").getAbsolutePath();
            BufferedReader br = IOUtils.getTextGzipReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
            try {
                String temp = null;
                String[] temps = null;
                bw1.write("Chr\tPos\tRef\tAlt(maximum 2 alternative alleles, which is seperated by \",\", e.g. A,C)\n");
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    bw.write(temps[0] + "\t" + temps[1] + "\n");
                    bw1.write(temps[0] + "\t" + temps[1] + "\t" + temps[3] + "\t" + temps[4] + "\n");
                }
                br.close();
                bw.flush();
                bw.close();
                bw1.flush();
                bw1.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Writing pos and posAllele files takes " + (endTime - startTime) + "ms");
    }

    public void parameter() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing parameters files ***********************************************************");
        parameterDir = new File(HapscannerDir, "parameter").getAbsolutePath();
        File paradir = new File(parameterDir);
        if (!paradir.exists()) paradir.mkdir();

        File posdir = new File(new File(HapscannerDir, "pos").getAbsolutePath());
        File posAlleledir = new File(new File(HapscannerDir, "posAllele").getAbsolutePath());
        if (!posdir.exists() || !posAlleledir.exists()) {
            posdir.mkdir();
            posAlleledir.mkdir();
            this.poswithAllele();
        }

        posDir = new File(HapscannerDir, "pos").getAbsolutePath();
        posAlleleDir = new File(HapscannerDir, "posAllele").getAbsolutePath();

        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(parameterDir, plate + "_parameter_chr" + chr + ".txt").getAbsolutePath());
                bw.write("@App:\tHapScanner\n" +
                        "@Author:\tFei Lu\n" +
                        "@Email:\tflu@genetics.ac.cn; dr.lufei@gmail.com\n" +
                        "@Homepage:\thttps://plantgeneticslab.weebly.com/\n" +
                        "\n" +
                        "#HapScanner is used to perform genotyping of diploid species from whole genome sequenceing data, based on an existing genetic variation library.\n" +
                        "#To run and pipeline, the machine should have both Java 8 and samtools installed. The lib directory should stay with TIGER.jar in the same folder.\n" +
                        "#Command line example. java -Xmx100g -jar TIGER.jar -a HapScanner -p parameter_hapscanner.txt > log.txt &\n" +
                        "#To specify options, please edit the the parameters below. Also, please keep the order of parameters.\n" +
                        "\n" +
                        "#Parameter 1: The taxaRefBam file containing information of taxon and its corresponding refernece genome and bam files. The bam file should have .bai file in the same folder\n" +
                        "#If one taxon has n bam files, please list them in n rows.\n");
                bw.write(new File(HapscannerDir, "taxaRefBAM/" + plate + "_taxaRefBAM_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library. \n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n");
                bw.write(new File(HapscannerDir, "posAllele/posAllele_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n");
                bw.write(new File(HapscannerDir, "pos/pos_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 4: The chromosome which will be scanned.\n" +
                        chr + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        new File(samtools).getAbsolutePath() + "\n" +
                        "\n" +
                        "#Parameter 7: Number of threads\n" +
                        threads + "\n" +
                        "\n" +
                        "#Parameter 8: The directory of output\n");
                bw.write(HapscannerDir + "/output/" + plate + "\n");
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis();
        System.out.println("******* Writing parameter files takes " + (endTime - startTime) + "ms");
    }

    public void taxaRefBAM() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing taxaRefBam files *************************************************************************");
        ArrayList<String> files = new ArrayList<>();
        BamDir = new File(workingDirS, subDir[2]).getAbsolutePath();
        taxaRefBAMDir = new File(HapscannerDir, "taxaRefBAM").getAbsolutePath();
        File taxadir = new File(taxaRefBAMDir);
        if (!taxadir.exists()) taxadir.mkdir();
        File[] fs1 = new File(BamDir).listFiles();
        fs1 = IOUtils.listFilesEndsWith(fs1, bamsuffix);
        for (int j = 0; j < fs1.length; j++) {
            System.out.println(fs1[j].getAbsolutePath());
            files.add(fs1[j].getAbsolutePath());
            System.out.println(fs1[j]);
        }
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < files.size(); i++) {
            nameSet.add(files.get(i));
        }
        try {
            String[] namelist = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(taxaRefBAMDir, plate + "_taxaRefBAM_chr" + chr + ".txt").getAbsolutePath());
                bw.write("Taxa\tReference\tBamPath\n");
                for (int j = 0; j < namelist.length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(namelist[j].split("/")[namelist[j].split("/").length - 1].replace(bamsuffix, "")).append("\t");
                    sb.append(new File(ref, "chr" + chr + ".fa\t").getAbsolutePath());
                    sb.append(new File(namelist[j]).getAbsolutePath());
                    System.out.println(sb.toString());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Writing taxaRefBam files takes " + (endTime - startTime) + "ms");
    }

    public void mkFinalVCF() {
        Set<String> taxaSet = taxaBamsMap.keySet();
        String[] taxa = taxaSet.toArray(new String[taxaSet.size()]);
        Arrays.sort(taxa);
        String outfileS = new File(outputDirS, subDirS[2]).getAbsolutePath();
        outfileS = new File(outfileS, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf").getAbsolutePath();
        try {
            SimpleDateFormat sdf = new SimpleDateFormat("MM/dd/yyyy HH:mm:ss.SSS");
            Date dt = new Date();
            String S = sdf.format(dt);
            BufferedWriter bw = pgl.infra.utils.IOUtils.getTextWriter(outfileS);
            bw.write("##fileformat=VCFv4.1\n");
            bw.write("##fileDate=" + S.split(" ")[0] + "\n");
            bw.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
            bw.write("##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and alternate alleles in the order listed\">\n");
            bw.write("##FORMAT=<ID=GL,Number=G,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1, or  0/0, 0/1, 0/2, 1/1, 1/2, 2/2 if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
            bw.write("##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n");
            bw.write("##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed starting with REF\">\n");
            bw.write("##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n");
            bw.write("##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or AA,AB,AC,BB,BC,CC if 2 alt alleles\">\n");
            bw.write("##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n");
            bw.write("##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n");
            bw.write("##ALT=<ID=DEL,Description=\"Deletion\">\n");
            bw.write("##ALT=<ID=INS,Description=\"Insertion\">\n");
            StringBuilder sb = new StringBuilder("#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT");
            for (int i = 0; i < taxa.length; i++) {
                sb.append("\t").append(taxa[i]);
            }
            bw.write(sb.toString());
            bw.newLine();
            BufferedReader[] brs = new BufferedReader[taxa.length];
            String indiVCFFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            for (int i = 0; i < brs.length; i++) {
                String indiVCFFileS = new File(indiVCFFolderS, taxa[i] + ".chr" + PStringUtils.getNDigitNumber(3, chr) + ".indi.vcf").getAbsolutePath();
                brs[i] = new BufferedReader(new FileReader(indiVCFFileS), 4096);
            }
            BufferedReader br = null;
            if (posAlleleFileS.endsWith(".gz")) {
                br = pgl.infra.utils.IOUtils.getTextGzipReader(posAlleleFileS);
            } else {
                br = pgl.infra.utils.IOUtils.getTextReader(posAlleleFileS);
            }

            String temp = br.readLine();
            List<String> temList = null;
            int cnt = 0;
            String[] genoArray = new String[brs.length];
            while ((temp = br.readLine()) != null) {
                sb.setLength(0);
                temList = PStringUtils.fastSplit(temp);
                sb.append(temList.get(0)).append("\t").append(temList.get(1)).append("\t").append(temList.get(0)).append("-").append(temList.get(1)).append("\t");
                sb.append(temList.get(2)).append("\t").append(temList.get(3)).append("\t.\t.\t");
                for (int i = 0; i < brs.length; i++) {
                    genoArray[i] = brs[i].readLine();
                }
                sb.append(this.getInfo(genoArray, temList.get(3))).append("\tGT:AD:GL");
                for (int i = 0; i < genoArray.length; i++) {
                    sb.append("\t").append(genoArray[i]);
                }
                bw.write(sb.toString());
                bw.newLine();
                cnt++;
                if (cnt % 1000000 == 0) System.out.println(String.valueOf(cnt) + " SNPs output to " + outfileS);
            }
            bw.flush();
            bw.close();
            br.close();
            for (int i = 0; i < brs.length; i++) {
                brs[i].close();
            }
            for (int i = 0; i < taxa.length; i++) {
                String indiVCFFileS = new File(indiVCFFolderS, taxa[i] + ".chr" + PStringUtils.getNDigitNumber(3, chr) + ".indi.vcf").getAbsolutePath();
                new File(indiVCFFileS).delete();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        new File(outputDirS, subDirS[0]).delete();
        new File(outputDirS, subDirS[1]).delete();
        System.out.println("Final VCF is completed at " + outfileS);
    }

    private String getInfo(String[] genoArray, String altList) {
        int dp = 0;
        int nz = 0;
        int nAlt = PStringUtils.fastSplit(altList, ",").size();
        int[] adCnt = new int[1 + nAlt];
        int[] acCnt = new int[1 + nAlt];
        int[][] gnCnt = new int[1 + nAlt][1 + nAlt];
        int ht = 0;
        List<String> tempList = null;
        List<String> temList = null;
        for (int i = 0; i < genoArray.length; i++) {
            if (genoArray[i].startsWith(".")) {
                nz++;
                continue;
            }
            tempList = PStringUtils.fastSplit(genoArray[i], ":");
            temList = PStringUtils.fastSplit(tempList.get(1), ",");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                dp += c;
                adCnt[j] += c;
            }
            temList = PStringUtils.fastSplit(tempList.get(0), "/");
            for (int j = 0; j < temList.size(); j++) {
                int c = Integer.parseInt(temList.get(j));
                acCnt[c]++;
            }
            int index1 = Integer.parseInt(temList.get(0));
            int index2 = Integer.parseInt(temList.get(1));
            gnCnt[index1][index2]++;
            if (index1 != index2) ht++;
        }
        nz = genoArray.length - nz;
        int sum = 0;
        for (int i = 0; i < acCnt.length; i++) {
            sum += acCnt[i];
        }
        float maf = (float) ((double) acCnt[0] / sum);
        if (maf > 0.5) maf = (float) (1 - maf);
        StringBuilder sb = new StringBuilder();
        sb.append("DP=").append(dp).append(";NZ=").append(nz).append(";AD=");
        for (int i = 0; i < adCnt.length; i++) {
            sb.append(adCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";AC=");
        for (int i = 1; i < acCnt.length; i++) {
            sb.append(acCnt[i]).append(",");
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";GN=");
        for (int i = 0; i < gnCnt.length; i++) {
            for (int j = i + 1; j < gnCnt.length; j++) {
                sb.append(gnCnt[i][j]).append(",");
            }
        }
        sb.deleteCharAt(sb.length() - 1);
        sb.append(";HT=").append(ht).append(";MAF=").append(maf);
        return sb.toString();
    }

    public void parseParameters(String infileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(infileS);
        List<String> pLineList = d.getFirstElement();
        taxaRefBamFileS = pLineList.get(0);
        posAlleleFileS = pLineList.get(1);
        posFileS = pLineList.get(2);
        String[] tem = pLineList.get(3).split(":");
        chr = Integer.valueOf(tem[0]);
        if (tem.length == 2) {
            tem = tem[1].split(",");
            this.regionStart = Integer.parseInt(tem[0]);
            this.regionEnd = Integer.parseInt(tem[1]) + 1;
        }
        this.combinedErrorRate = Double.parseDouble(pLineList.get(4));
        samtoolsPath = pLineList.get(5);
        this.nThreads = Integer.parseInt(pLineList.get(6));
        outputDirS = pLineList.get(7);
        new File(outputDirS).mkdir();
        try {
            BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(taxaRefBamFileS);
            String temp = br.readLine();
            List<String> l = new ArrayList<>();
            while ((temp = br.readLine()) != null) {
                l = PStringUtils.fastSplit(temp);
                String key = l.get(0);
                List<String> bamList = new ArrayList<>();
                for (int j = 0; j < l.size() - 2; j++) {
                    bamList.add(l.get(j + 2));
                }
                taxaBamsMap.put(key, bamList);
                taxaRefMap.put(key, l.get(1));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void scanIndiVCFByThreadPool() {
        this.creatFactorialMap();
        pgl.infra.table.RowTable<String> t = new pgl.infra.table.RowTable<>(posAlleleFileS);
        HashMap<Integer, String> posRefMap = new HashMap<>();
        HashMap<Integer, String[]> posAltMap = new HashMap<>();
        int[] positions = new int[t.getRowNumber()];
        for (int i = 0; i < t.getRowNumber(); i++) {
            positions[i] = t.getCellAsInteger(i, 1);
            posRefMap.put(positions[i], t.getCell(i, 2));
            String[] tem = t.getCell(i, 3).split(",");
            posAltMap.put(t.getCellAsInteger(i, 1), tem);
        }
        Set<String> taxaSet = taxaBamsMap.keySet();
        ArrayList<String> taxaList = new ArrayList(taxaSet);
        Collections.sort(taxaList);

//        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(taxaList.size(), this.nThreads);
        LongAdder counter = new LongAdder();
        ExecutorService pool = Executors.newFixedThreadPool(this.nThreads);
        List<Future<IndiVCF>> resultList = new ArrayList<>();
        for (int i = 0; i < taxaList.size(); i++) {
            String indiVCFFolderS = new File(outputDirS, subDirS[1]).getAbsolutePath();
            String indiVCFFileS = new File(indiVCFFolderS, taxaList.get(i) + ".chr" + PStringUtils.getNDigitNumber(3, chr) + ".indi.vcf").getAbsolutePath();
            List<String> bamPaths = taxaBamsMap.get(taxaList.get(i));
            StringBuilder sb = new StringBuilder(samtoolsPath);
            sb.append(" mpileup -A -B -q 20 -Q 20 -f ").append(this.taxaRefMap.get(taxaList.get(i)));
            for (int j = 0; j < bamPaths.size(); j++) {
                sb.append(" ").append(bamPaths.get(j));
            }
            sb.append(" -l ").append(posFileS).append(" -r ");
            sb.append(chr);
            if (this.regionStart != Integer.MIN_VALUE) {
                sb.append(":").append(this.regionStart).append("-").append(regionEnd - 1);
            }
            String command = sb.toString();
            IndiVCF idv = new IndiVCF(command, indiVCFFileS, posRefMap, posAltMap, positions, bamPaths, counter);
            Future<IndiVCF> result = pool.submit(idv);
            resultList.add(result);
        }
        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    class IndiVCF implements Callable<IndiVCF> {
        String command = null;
        String indiVCFFileS = null;
        HashMap<Integer, String> posRefMap = null;
        HashMap<Integer, String[]> posAltMap = null;
        int[] positions = null;
        List<String> bamPaths = null;
        LongAdder counter = null;

        public IndiVCF(String command, String indiVCFFileS, HashMap<Integer, String> posRefMap, HashMap<Integer, String[]> posAltMap, int[] positions, List<String> bamPaths, LongAdder counter) {
            this.command = command;
            this.indiVCFFileS = indiVCFFileS;
            this.posRefMap = posRefMap;
            this.posAltMap = posAltMap;
            this.positions = positions;
            this.bamPaths = bamPaths;
            this.counter = counter;
        }


        @Override
        public IndiVCF call() throws Exception {
            try {
                Runtime rt = Runtime.getRuntime();
                Process p = rt.exec(command);
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));

//                BufferedReader bre = new BufferedReader(new InputStreamReader(p.getErrorStream()));
//                String temp = null;
//                while ((temp = bre.readLine()) != null) {
//                    if (temp.startsWith("[m")) continue;
//                    System.out.println(command);
//                    System.out.println(temp);
//                }
//                bre.close();

                BufferedWriter bw = pgl.infra.utils.IOUtils.getTextWriter(indiVCFFileS);
                String current = br.readLine();
                List<String> currentList = null;
                int currentPosition = -1;
                if (current != null) {
                    currentList = PStringUtils.fastSplit(current);
                    currentPosition = Integer.parseInt(currentList.get(1));
                }
                StringBuilder sb = new StringBuilder();
                for (int i = 0; i < positions.length; i++) {
                    if (current == null) {
                        bw.write("./.");
                        bw.newLine();
                    } else {
                        if (positions[i] == currentPosition) {
                            String ref = posRefMap.get(currentPosition);
                            String[] alts = posAltMap.get(currentPosition);
                            char[] alleleC = new char[alts.length + 1];
                            alleleC[0] = ref.charAt(0);
                            for (int j = 0; j < alts.length; j++) {
                                if (alts[j].startsWith("I") || alts[j].startsWith("<I")) {
                                    alleleC[j + 1] = '+';
                                } else if (alts[j].startsWith("D") || alts[j].startsWith("<D")) {
                                    alleleC[j + 1] = '-';
                                } else {
                                    alleleC[j + 1] = alts[j].charAt(0);
                                }
                            }
                            int[] cnts = new int[alts.length + 1];
                            sb.setLength(0);
                            for (int j = 0; j < bamPaths.size(); j++) {
                                sb.append(currentList.get(4 + j * 3));
                            }
                            StringBuilder ssb = new StringBuilder();
                            int curIndex = 0;
                            for (int j = 0; j < sb.length(); j++) {
                                char cChar = sb.charAt(j);
                                if (cChar == '+') {
                                    ssb.append(sb.subSequence(curIndex, j + 1));
                                    curIndex = j + 2 + Character.getNumericValue(sb.charAt(j + 1));
                                } else if (cChar == '-') {
                                    ssb.append(sb.subSequence(curIndex, j + 1));
                                    curIndex = j + 2 + Character.getNumericValue(sb.charAt(j + 1));
                                }
                            }
                            ssb.append(sb.subSequence(curIndex, sb.length()));
                            sb = ssb;
                            String s = sb.toString().toUpperCase();
                            for (int j = 0; j < s.length(); j++) {
                                char cChar = s.charAt(j);
                                if (cChar == '.' || cChar == ',') {
                                    cnts[0]++;
                                    continue;
                                }
                                for (int k = 1; k < alleleC.length; k++) {
                                    if (cChar == alleleC[k]) cnts[k]++;
                                }
                            }
                            for (int j = 1; j < alleleC.length; j++) {
                                if (alleleC[j] == '+') cnts[0] = cnts[0] - cnts[j];
                                else if (alleleC[j] == '-') cnts[0] = cnts[0] - cnts[j];
                            }
                            String vcf = getGenotype(cnts);
                            bw.write(vcf);
                            bw.newLine();
                            current = br.readLine();
                            if (current != null) {
                                currentList = PStringUtils.fastSplit(current);
                                currentPosition = Integer.parseInt(currentList.get(1));
                            }
                        } else if (positions[i] < currentPosition) {
                            bw.write("./.");
                            bw.newLine();
                        } else {
                            System.out.println("Current position is greater than pileup position. It should not happen. Program quits");
                            System.exit(1);
                        }
                    }
                }
                p.waitFor();
                bw.flush();
                bw.close();
                br.close();
            } catch (Exception ee) {
                ee.printStackTrace();
            }
            counter.increment();
            int cnt = counter.intValue();
            if (cnt % 10 == 0)
                System.out.println("Finished individual genotyping in " + String.valueOf(cnt) + " taxa. Total: " + String.valueOf(taxaBamsMap.size()));
            return this;
        }

    }

    private String getGenotype(int[] cnt) {
        int n = cnt.length * (cnt.length + 1) / 2;
        int[] likelihood = new int[n];
        int sum = 0;
        for (int i = 0; i < cnt.length; i++) sum += cnt[i];
        if (sum == 0) return "./.";
        else if (sum > this.maxFactorial) {
            double portion = (double) this.maxFactorial / sum;
            for (int i = 0; i < cnt.length; i++) {
                cnt[i] = (int) (cnt[i] * portion);
            }
            sum = this.maxFactorial;
        }
        double coe = this.factorialMap.get(sum);
        for (int i = 0; i < cnt.length; i++) coe = coe / this.factorialMap.get(cnt[i]);
        double max = Double.MAX_VALUE;
        int a1 = 0;
        int a2 = 0;
        for (int i = 0; i < cnt.length; i++) {
            for (int j = i; j < cnt.length; j++) {
                int index = (j * (j + 1) / 2) + i;
                double value = Double.MAX_VALUE;
                if (i == j) {
                    value = -Math.log10(coe * Math.pow((1 - 0.75 * this.combinedErrorRate), cnt[i]) * Math.pow(this.combinedErrorRate / 4, (sum - cnt[i])));
                } else {
                    value = -Math.log10(coe * Math.pow((0.5 - this.combinedErrorRate / 4), cnt[i] + cnt[j]) * Math.pow(this.combinedErrorRate / 4, (sum - cnt[i] - cnt[j])));
                }
                if (value < max) {
                    max = value;
                    a1 = i;
                    a2 = j;
                }
                likelihood[index] = (int) Math.round(value);
            }
        }
        StringBuilder sb = new StringBuilder();
        sb.append(a1).append("/").append(a2).append(":");
        for (int i = 0; i < cnt.length; i++) sb.append(cnt[i]).append(",");
        sb.deleteCharAt(sb.length() - 1);
        sb.append(":");
        for (int i = 0; i < likelihood.length; i++) sb.append(likelihood[i]).append(",");
        sb.deleteCharAt(sb.length() - 1);
        return sb.toString();
    }

    private void creatFactorialMap() {
        this.factorialMap = HashIntDoubleMaps.getDefaultFactory().newMutableMap();
        for (int i = 0; i < this.maxFactorial + 1; i++) {
            this.factorialMap.put(i, factorial(i));
        }
    }

    public void mkDir() {
        for (int i = 0; i < subDirS.length; i++) {
            File f = new File(outputDirS, subDirS[i]);
            f.mkdir();
        }
    }

//    private void callsnp() {
//        genotypeDir = this.genotypeDir;
//        String bamDir = new File(workingDirS, subDir[2]).getAbsolutePath();
//        String pos = "no";
//        StringBuilder sb = new StringBuilder();
//        sb.append("java -jar /data1/home/xiaohan/jar/Hapscanner_0518.jar -g ");
//        sb.append(new File(genotypeDir).getAbsolutePath());
//        sb.append(" -b " + new File(bamDir).getAbsolutePath());
//        sb.append(" -o " + new File(HapscannerDir).getAbsolutePath());
//        sb.append(" -bs " + bamsuffix);
//        sb.append(" -gs " + genotypesuffix);
//        sb.append(" -t " + threads);
//        sb.append(" -pos no -p " + plate);
////        sb.append(" -r " + ref);
//        sb.append(" -samtools " + new File(samtools).getAbsolutePath());
//        String command = sb.toString();
//        System.out.println(command);
//        try {
//            File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
//            String[] cmdarry = {"/bin/bash", "-c", command};
//            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//            p.waitFor();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }

    private void sortbam(String[] args) {
//        workingDirS = args[0];
//        samtools = "~/miniconda3/bin/samtools";
        String inputDir = new File(workingDirS, subDir[1]).getAbsolutePath();
        File[] fs = IOUtils.listFilesEndsWith(new File(inputDir).listFiles(), ".bam");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName().replace(".bam", "");
            if (!nameSet.contains(name)) {
                nameSet.add(name);
            }
        }
        StringBuilder sb = new StringBuilder();
        nameSet.parallelStream().forEach(f -> {
            try {
                sb.setLength(0);
                String bam = new File(workingDirS, subDir[1] + "/" + f + ".bam").getAbsolutePath();
                String sortedbam = new File(workingDirS, subDir[2] + "/" + f + ".sorted.bam").getAbsolutePath();
                sb.append(samtools + " sort " + bam + " > " + sortedbam + "\n");
                sb.append(samtools + " index " + sortedbam + "\n");
                sb.append("rm " + bam + "\n");
                String command = sb.toString();
                File dir = new File(new File(workingDirS).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    private void bwa(String[] args) {
//        workingDirS = args[0];
//        threads = args[1];
//        bwa = "~/miniconda3/bin/bwa";
//        bwalib = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz";
//        samtools = "~/miniconda3/bin/samtools";
        String inputDir = new File(workingDirS, subDir[0]).getAbsolutePath();
        File[] fs = IOUtils.listFilesEndsWith(new File(inputDir).listFiles(), "1.fq.gz");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName().replace("1.fq.gz", "");
            nameSet.add(name);
        }
        StringBuilder sb = new StringBuilder();
        nameSet.stream().forEach(f -> {
            sb.setLength(0);
            String fq1 = new File(workingDirS, subDir[0] + "/" + f + "1.fq.gz").getAbsolutePath();
            String fq2 = new File(workingDirS, subDir[0] + "/" + f + "2.fq.gz").getAbsolutePath();
            String sam = new File(workingDirS, subDir[1] + "/" + f.split("_")[0] + ".sam").getAbsolutePath();
            String bam = new File(workingDirS, subDir[1] + "/" + f.split("_")[0] + ".bam").getAbsolutePath();
            String title = "@RG" + "\\" + "tID:" + f.split("_")[0] + "\\" + "tPL:illumina" + "\\" + "t" + f.split("_")[0] + "\\" + "tLB:library";
            System.out.println(title);
            sb.append(bwa + " mem -t " + threads + " -R " + "\'" + title + "\'" + " " + bwalib + " " + fq1 + " " + fq2 + " > " + sam + "\n");
            sb.append(samtools + " view -S -b " + sam + " > " + bam + "\n");
            sb.append("rm " + sam + "\n");
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(workingDirS).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

//    private void parsefastqs() {
//        pgl.app.grt.LibraryInfo li = new LibraryInfo(barcodeFileS, libraryFastqMapFileS, this.cutter1, this.cutter2);
//        String tagBySampleDirS = new File(workingDirS, subDir[0]).getAbsolutePath();
//        TagParser tp = new TagParser(li);
//        tp.parseFastq(tagBySampleDirS);
//
//        File[] fs = new File(tagBySampleDirS).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, "tp");
//        HashSet<String> nameSet = new HashSet();
//        for (int i = 0; i < fs.length; i++) {
//            if (fs[i].isHidden()) {
//                continue;
//            }
//            String Name = fs[i].getName().replace(".tp", "");
//            nameSet.add(Name);
//            System.out.println(Name);
//        }
//        nameSet.parallelStream().forEach(f -> {
//            try {
//                String r1FastqFileS = new File(tagBySampleDirS, f + ".R1.fq").getAbsolutePath();
//                String r2FastqFileS = new File(tagBySampleDirS, f + ".R2.fq").getAbsolutePath();
//                String infile = new File(tagBySampleDirS, f + ".tp").getAbsolutePath();
//                TagAnnotations tc = new TagAnnotations(infile);
//                tc.writeFastqFile(r1FastqFileS, r2FastqFileS);
//                StringBuilder sb = new StringBuilder();
//                sb.append("bgzip " + r1FastqFileS + "\n");
//                sb.append("bgzip " + r2FastqFileS + "\n");
//                String command = sb.toString();
//                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
//                String[] cmdarry = {"/bin/bash", "-c", command};
//                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                p.waitFor();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//    }

    private void parsefastqs() {
        LibraryOfGRT li = new LibraryOfGRT(barcodeFileS, libraryFastqMapFileS, cutter1, cutter2);
        li.splitBarcode(new File(workingDirS, subDir[0]).getAbsolutePath());
        StringBuilder sb = new StringBuilder();
        sb.append("rm NA*");
        try {
            String command = sb.toString();
            File dir = new File(new File(workingDirS, subDir[0]).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        File[] fs = new File(workingDirS, subDir[0]).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "fq");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(new File(String.valueOf(fs[i])).getAbsolutePath());
            System.out.println(new File(String.valueOf(fs[i])).getAbsolutePath());
        }

        nameSet.parallelStream().forEach(f -> {
            try {
                StringBuilder sb1 = new StringBuilder();
                sb1.append("bgzip " + f);
                File dir = new File(new File(workingDirS, subDir[0]).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", sb1.toString()};
                System.out.println(sb1.toString());
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });

    }

    private String getChimericRemovedRead(String cutter1, String cutter2, String read, int barcodeLength) {
        read = read.substring(barcodeLength, read.length());
        int index1 = read.indexOf(cutter1);
        int index2 = read.indexOf(cutter2);
        if (index1 < 0) {
            if (index2 < 0) {
                return read;
            } else {
                return read.substring(0, index2);
            }
        } else {
            if (index2 < 0) {
                return read.substring(0, index1);
            } else {
                if (index1 < index2) return read.substring(0, index1);
                else return read.substring(0, index2);
            }
        }
    }

    public static void main(String[] args) {
        new heterozygosity(args);
    }
}
