package xiaohan.eQTL.pipline;

import xiaohan.rareallele.GeneFeature;
import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

public class vcf {
    public vcf(String[] args) {
//        this.getFilteredVCF();
//        this.get360VCF();
//        this.getTSSbed(args);
        this.getgenesbed(args);
    }

    public void getTSSbed(String... args) {
        String infile = args[0];
        String outfile = args[1];
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        GeneFeature gf = new GeneFeature(infile);
        StringBuilder sb = new StringBuilder();
        int TSS = -1;
        int chr = -1;
        int start = -1;
        String geneName = null;
        try {
            for (int i = 0; i < gf.getGeneNumber(); i++) {
                sb.setLength(0);
                if (gf.getGeneStrand(i) == 1) {
                    TSS = gf.getGeneStart(i);
                } else {
                    TSS = gf.getGeneEnd(i);
                }
                chr = gf.getGeneChromosome(i);
                start = TSS - 1;
                geneName = gf.getGeneName(i);
                sb.append("chr" + chr + "\t" + start + "\t" + TSS + "\t" + geneName);
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getgenesbed(String... args) {
        String infile = args[0];
        String outfile = args[1];
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        GeneFeature gf = new GeneFeature(infile);
        StringBuilder sb = new StringBuilder();
        int chr = -1;
        int start = -1;
        int end = -1;
        String geneName = null;
        try {
            for (int i = 0; i < gf.getGeneNumber(); i++) {
                sb.setLength(0);
                start = gf.getGeneStart(i);
                end = gf.getGeneEnd(i);
                chr = gf.getGeneChromosome(i);
                geneName = gf.getGeneName(i);
                sb.append("chr" + chr + "\t" + start + "\t" + end + "\t" + geneName);
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void get360VCF() {
        String inputDir = "/data2/xiaohan/genotype/genotype_eQTL";
        String outputDir = "/data2/xiaohan/genotype/genotype360";
        HashSet nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(chr);
        }
        nameSet.parallelStream().forEach(f -> {
            String chr = String.valueOf(f);
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputDir, chr + ".346.B18.recode.vcf.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, chr + ".360.B18.recode.vcf").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                HashMap<String, Integer> NameIndexMap = new HashMap<>();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        for (int j = 9; j < temps.length; j++) {
                            System.out.println(temps[j]);
                            NameIndexMap.put(temps[j], j);
                        }
                        for (int j = 0; j < 9; j++) {
                            bw.write(temps[j] + "\t");
                        }
                        for (int j = 1; j <= 360; j++) {
                            if ((j + "").length() == 1) {
                                bw.write("E00" + j + "\t");
                            }
                            if ((j + "").length() == 2) {
                                bw.write("E0" + j + "\t");
                            }
                            if ((j + "").length() == 3) {
                                bw.write("E" + j + "\t");
                            }
                        }
                        bw.write("\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    for (int j = 0; j < 9; j++) {
                        bw.write(temps[j] + "\t");
                    }
                    String Name = null;
                    for (int j = 1; j <= 360; j++) {
                        if ((j + "").length() == 1) {
                            Name = "E00" + j;
                        } else if ((j + "").length() == 2) {
                            Name = "E0" + j;
                        } else if ((j + "").length() == 3) {
                            Name = "E" + j;
                        }
                        System.out.println(Name);
                        if (j % 25 == 0) {
                            bw.write(temps[NameIndexMap.get("E359")] + "\t");
                        } else {
                            bw.write(temps[NameIndexMap.get(Name)] + "\t");
                        }
                    }
                    bw.write("\n");
                }
                br.close();
                bw.flush();
                bw.close();
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip " + chr + ".360.B18.recode.vcf\"");
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File(outputDir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getFilteredVCF() {
        System.out.println("This is the beginning of the program............................");
        String inputdir = "/data2/junxu/genotypeMaf005";
        String outputdir = "/data2/xiaohan/genotype/genotype_eQTL";
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(String.valueOf(chr));
            System.out.println(chr);
        }
        nameSet.stream().forEach(f -> {
            try {
                System.out.println("-----------------------This is calculating file : chr " + f + "------------------------------");
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputdir, f + ".346.B18.recode.vcf.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, f + ".346.B18.recode.vcf").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String[] tems = null;
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 50000 == 0) {
                        System.out.println("This is counting line : " + countline);
                    }
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    int number = 0;
                    for (int i = 0; i < temps.length - 9; i++) {
                        tems = temps[i + 9].split(":");
                        if (tems[0].equals("0/1")) number++;
                    }
                    if (number <= 43) {
                        bw.write(temp + "\n");
                    } else continue;
                }
                br.close();
                bw.flush();
                bw.close();
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip " + f + ".346.B18.recode.vcf");
                String command = sb.toString();
                File dir = new File(new File(outputdir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public static void main(String[] args) {
        new vcf(args);
    }
}
