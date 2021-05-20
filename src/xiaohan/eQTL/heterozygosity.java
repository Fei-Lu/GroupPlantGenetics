package xiaohan.eQTL;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import org.apache.commons.cli.*;
import pgl.app.grt.LibraryInfo;
import pgl.app.grt.TagAnnotations;
import pgl.app.grt.TagParser;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import tech.tablesaw.api.Row;
import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class heterozygosity {
    int index = 1;
    String barcodeFileS = "/data1/home/xiaohan/nam/20210515GBStest" + index + "/barcode_GBS" + index + ".txt";
    String libraryFastqMapFileS = "/data1/home/xiaohan/nam/20210515GBStest" + index + "/libraryFastqMap_GBS" + index + ".txt";
//
//    String barcodeFileS = "/Users/yxh/Documents/NAM_library/Lib_GBS/source/20180601_GBSLibrarybarcode.txt";
//    String libraryFastqMapFileS = "/Users/yxh/Documents/NAM_library/Lib_GBS/source/LibraryFastqMap.txt";

    String cutter1 = "GGATCC";
    String cutter2 = "CCGG";
    String workingDirS = "/data1/home/xiaohan/nam/20210515GBStest" + index + "/";
//    String workingDirS = "/Users/yxh/Documents/NAM_library/Lib_GBS/source";

    String[] subDir = {"subFastqs", "bam", "sortedbam", "IBS", "heter", "sum"};
    String HapscannerDir = "/data1/home/xiaohan/hapscanner/";
    String kinship = "/data1/home/xiaohan/nam/20210515GBStest" + index + "/KinshipMap.txt";
    Multimap<String, String> kinshipMap = ArrayListMultimap.create();
    HashSet<String> HBSet = new HashSet<>();

    String bwa = "/data1/home/xiaohan/miniconda3/bin/bwa";
    String bwalib = "/data1/publicData/wheat/reference/v1.0/ABD/bwaLib/abd_iwgscV1.fa.gz";
    String samtools = "/data1/home/xiaohan/miniconda3/bin/samtools";

    String ref = "/data1/home/xiaohan/hapscanner/ref";
    String plate = "20210515GBStest" + index;
    String genotypeDir = "/data1/home/xiaohan/nam/Parent/VCF";
    String bamsuffix = ".sorted.bam";
    String genotypesuffix = ".vcf";
    String threads = "10";


    public heterozygosity(String... args) {
        this.mkdir();
//        this.trimmomatics();
//        this.parsefastqs();
//        this.bwa();
//        this.sortbam();
//        this.callsnp();
//        this.interset();
        this.getkinship();
//        this.getheterozygosity();
        this.getIBS();
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
        for (int i = 0; i < 42; i++) {
            chrSet.add(i + 1);
        }

        chrSet.parallelStream().forEach(f -> {
            try {
                int dep1 = 0;
                int dep2 = 0;
                String vcf = new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String genovcf = new File(genotypeDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String outfile = new File(workingDirS, subDir[4] + "/chr" + PStringUtils.getNDigitNumber(3, f) + "_het.txt").getAbsolutePath();
                String outfile1 = new File(workingDirS, subDir[4] + "/chr" + PStringUtils.getNDigitNumber(3, f) + "_plot.txt").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(vcf);
                BufferedReader br1 = IOUtils.getTextReader(genovcf);
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
                String temp = null;
                String temp1 = null;
                temp = br.readLine();
                temp1 = br1.readLine();
                while ((!temp.startsWith("#C"))) {
                    if (temp.startsWith("#C")) break;
                    temp = br.readLine();
                }
                System.out.println(temp);
                while ((!temp1.startsWith("#C"))) {
                    if (temp1.startsWith("#C")) break;
                    temp1 = br1.readLine();
                }
                System.out.println(temp1);
                HashMap<Integer, String> HBMap = new HashMap<>();
                HashMap<String, Integer> PMap = new HashMap<>();
                String[] temps = temp.split("\t");
                String[] temps1 = temp1.split("\t");
                StringBuilder sb = new StringBuilder();
                sb.append("Chr\tPos\t");
                for (int i = 9; i < temps.length; i++) {
                    HBMap.put(i, temps[i]);
                    sb.append(temps[i] + "\t");
                }
                bw1.write(sb.toString().replaceAll("\\s+$", ""));
                bw1.newLine();
                for (int i = 9; i < temps1.length; i++) {
                    PMap.put(temps1[i], i);
                }
                int[][] count = new int[temps.length - 9][4];
                String[] HBnames = Arrays.copyOfRange(temps, 9, temps.length);
                while ((temp = br.readLine()) != null) {
                    temp1 = br1.readLine();
                    temps = temp.split("\t");
                    temps1 = temp1.split("\t");
                    bw1.write(temps[0] + "\t" + temps[1] + "\t");
                    for (int i = 9; i < temps.length; i++) {
                        int het = -1;
//                        System.out.println(HBMap.get(i));
                        String[] P = kinshipMap.get(HBMap.get(i)).toArray(new String[0]);
//                        System.out.println(P[0]);
                        if (temps1[PMap.get(P[0])].startsWith("./.") || temps1[PMap.get(P[1])].startsWith("./.")) {
                            bw1.write("NA\t");
                            continue;
                        }
                        if (temps1[PMap.get(P[0])].startsWith("0/1") || temps1[PMap.get(P[1])].startsWith("0/1")) {
                            bw1.write("NA\t");
                            continue;
                        }
                        if (temps1[PMap.get(P[0])].split(":")[0].equals(temps1[PMap.get(P[1])].split(":")[0])) {
                            bw1.write("NA\t");
                            continue;
                        } else if (!temps1[PMap.get(P[0])].split(":")[0].equals(temps1[PMap.get(P[1])].split(":")[0])) {
                            if (temps[i].startsWith("./.")) {
                                count[i - 9][2] += 1;
                                het = -1;
                            } else {
                                System.out.println(temps[i]);
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
                                    bw1.write(het + "\t");
                                    continue;
                                } else continue;
                            }
                        }
                    }
                    bw1.newLine();
                }
                bw1.flush();
                bw1.close();
                br.close();
                br1.close();
                for (int i = 0; i < HBnames.length; i++) {
                    bw.write(HBnames[i] + "\t" + count[i][0] + "\t" + count[i][1] + "\t" + count[i][2]);
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    private void getIBS() {
        String inputDir = new File(HapscannerDir, "output/" + plate + "/VCF").getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            chrSet.add(i + 1);
        }
        chrSet.stream().forEach(f -> {
            try {
                String vcf1 = new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String vcf2 = new File(genotypeDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String ibsOutfileS = new File(workingDirS,subDir[3]+"/check" + f + ".txt").getAbsolutePath();
                GenotypeGrid g1 = new GenotypeGrid(vcf1, GenoIOFormat.VCF);
                GenotypeGrid g2 = new GenotypeGrid(vcf2, GenoIOFormat.VCF);
                GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
                SumTaxaDivergence std = new SumTaxaDivergence(g);
                std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
                g.getIBSDistanceMatrix();
                RowTable<String> rt = new RowTable<>(ibsOutfileS);
                for (int i = 0; i < rt.getRowNumber(); i++) {
                    for (int j = 0; j < rt.getColumnNumber(); j++) {
                        if(rt.getCell(i,j).startsWith("N")){
                            rt.setCell(i,j,"1");
                        }
                    }
                }
                List<String> HBset = new ArrayList<>();
                List<String> Pset = new ArrayList<>();
                String[] header = rt.getHeader().toArray(new String[0]);
                for (int i = 0; i < header.length; i++) {
                    if (header[i].startsWith("HB")) {
                        HBset.add(header[i]);
                    } else Pset.add(header[i]);
                }
                BufferedWriter bw = IOUtils.getTextWriter(new File(workingDirS, subDir[subDir.length - 1] + "/sum" + f + ".txt").getAbsolutePath());
                for (int i = 0; i < HBset.size(); i++) {
                    String HB = HBset.get(i);
                    Collection<String> P = kinshipMap.get(HB);
                    String[] pa = rt.getColumn(rt.getColumnIndex("Dxy")).toArray(new String[0]);
                    System.out.println(pa[0]);
                    double[] IBS = rt.getColumnAsDoubleArray(rt.getColumnIndex(HB));
                    double[] IBSsub = Arrays.copyOfRange(IBS, HBset.size(), header.length-1);
                    double[] minMin = SNPmappingInGene.minMin(IBSsub);
                    System.out.println(HBset.size());
                    System.out.println(Pset.size());
                    System.out.println(IBSsub.length);
                    System.out.println(header.length);
                    String p1 = header[HBset.size() + (int) minMin[2] +1];
                    String p2 = header[HBset.size() + (int) minMin[3] +1];
                    System.out.println(HB + " of " + p1 + " & " + p2);
                    if (P.contains(p1) && P.contains(p2)) {
                        bw.write(HB + "\tcorrect\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\n");
                        continue;
                    } else {
                        bw.write(HB + "\twrong\t" + p1 + "\t" + p2 + "\t" + minMin[0] + "\t" + minMin[1] + "\n");
                        continue;
                    }
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
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

    private void interset() {
        String inputDir = new File(HapscannerDir, "output/" + plate).getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            chrSet.add(i + 1);
        }
        chrSet.stream().forEach(f -> {
            try {
                String vcf = new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String genovcf = new File(genotypeDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf.gz").getAbsolutePath();
                String subvcf = new File(workingDirS, subDir[3] + "/chr" + PStringUtils.getNDigitNumber(3, f) + "sub.vcf").getAbsolutePath();
                String subgenovcf = new File(workingDirS, subDir[3] + "/chr" + PStringUtils.getNDigitNumber(3, f) + "subgeno.vcf").getAbsolutePath();
                BufferedReader br1 = IOUtils.getTextReader(vcf);
                BufferedWriter bw1 = IOUtils.getTextWriter(subvcf);
                BufferedReader br2 = IOUtils.getTextReader(genovcf);
                BufferedWriter bw2 = IOUtils.getTextWriter(subgenovcf);
                ArrayList<String> sitelist = new ArrayList<>();
                String temp = null;
                String[] temps = null;
                int number = 0;
                int total = 0;
                while ((temp = br1.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        if (temp.startsWith("#C")) {
                            total = (int) ((temp.split("\t").length - 9) * 0.6);
                        }
                        bw1.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    number = 0;
                    for (int i = 9; i < temps.length; i++) {
                        if (temps[i].startsWith("./.")) {
                            number++;
                        }
                    }
                    if (number < total) {
                        sitelist.add(temps[1]);
                        bw1.write(temp + "\n");
                    }
                }
                br1.close();
                bw1.flush();
                bw1.close();

                while ((temp = br2.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw2.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    if (sitelist.contains(temps[1])) {
                        bw2.write(temp + "\n");
                    }
                }
                br2.close();
                bw2.flush();
                bw2.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    private void callsnp() {
        genotypeDir = this.genotypeDir;
        String bamDir = new File(workingDirS, subDir[2]).getAbsolutePath();
        String pos = "no";
        StringBuilder sb = new StringBuilder();
        sb.append("java -jar /data1/home/xiaohan/jar/Hapscanner_0518.jar -g ");
        sb.append(new File(genotypeDir).getAbsolutePath());
        sb.append(" -b " + new File(bamDir).getAbsolutePath());
        sb.append(" -o " + new File(HapscannerDir).getAbsolutePath());
        sb.append(" -bs " + bamsuffix);
        sb.append(" -gs " + genotypesuffix);
        sb.append(" -t " + threads);
        sb.append(" -pos no -p " + plate);
//        sb.append(" -r " + ref);
        sb.append(" -samtools " + new File(samtools).getAbsolutePath());
        String command = sb.toString();
        System.out.println(command);
        try {
            File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void sortbam() {
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
                String command = sb.toString();
                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    private void bwa() {
        String inputDir = new File(workingDirS, subDir[0]).getAbsolutePath();
        File[] fs = IOUtils.listFilesEndsWith(new File(inputDir).listFiles(), ".fq.gz");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < fs.length; i++) {
            String name = fs[i].getName().replace("R1.fq.gz", "");
            nameSet.add(name.replace("R1.fq.gz", ""));
        }
        StringBuilder sb = new StringBuilder();
        nameSet.stream().forEach(f -> {
            sb.setLength(0);
            String fq1 = new File(workingDirS, subDir[0] + "/" + f + "R1.fq.gz").getAbsolutePath();
            String fq2 = new File(workingDirS, subDir[0] + "/" + f + "R2.fq.gz").getAbsolutePath();
            String sam = new File(workingDirS, subDir[1] + "/" + f.split("_")[0] + ".sam").getAbsolutePath();
            String bam = new File(workingDirS, subDir[1] + "/" + f.split("_")[0] + ".bam").getAbsolutePath();
            String title = "@RG" + "\\" + "tID:" + f.split("_")[0] + "\\" + "tPL:illumina" + "\\" + "t" + f.split("_")[0] + "\\" + "tLB:library";
            System.out.println(title);
            sb.append(bwa + " mem -t " + threads + " -R " + "\'" + title + "\'" + " " + bwalib + " " + fq1 + " " + fq2 + " > " + sam + "\n");
            sb.append(samtools + " view -S -b " + sam + " > " + bam + "\n");
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    private void trimmomatics() {

    }

    private void parsefastqs() {
        pgl.app.grt.LibraryInfo li = new LibraryInfo(barcodeFileS, libraryFastqMapFileS, this.cutter1, this.cutter2);
        String tagBySampleDirS = new File(workingDirS, subDir[0]).getAbsolutePath();
        TagParser tp = new TagParser(li);
        tp.parseFastq(tagBySampleDirS);

        File[] fs = new File(tagBySampleDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "tp");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName().replace(".tp", "");
            nameSet.add(Name);
            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            try {
                String r1FastqFileS = new File(tagBySampleDirS, f + ".R1.fq").getAbsolutePath();
                String r2FastqFileS = new File(tagBySampleDirS, f + ".R2.fq").getAbsolutePath();
                String infile = new File(tagBySampleDirS, f + ".tp").getAbsolutePath();
                TagAnnotations tc = new TagAnnotations(infile);
                tc.writeFastqFile(r1FastqFileS, r2FastqFileS);
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip " + r1FastqFileS + "\n");
                sb.append("bgzip " + r2FastqFileS + "\n");
                String command = sb.toString();
                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
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
