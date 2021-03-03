package xiaohan.eQTL;

import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.IOFileFormat;
import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;

public class HapscannerParameters {

    public HapscannerParameters(String[] args) throws IOException, InterruptedException {
//        this.parseparameter();
//        File posDir = new File(new File("/data2/xiaohan/genotype/hapscanner/pos").getAbsolutePath());
//        File posAlleleDir = new File(new File("/data2/xiaohan/genotype/hapscanner/posAllele").getAbsolutePath());
//        if (!posDir.exists() || !posAlleleDir.exists()) {
//            this.poswithAllele();
//        }
//        this.parameter(args);
//        this.taxaRefBAM(args);
//        this.getsubgenotype(args);
//        this.getsubRNAgenotype(args);

//        this.getIsec(args);
//        this.getMergedVCF(a
//        rgs);
//        this.getsortedVCF(args);
//        this.getIBdistane(args);
//        this.getDensityIBS(args);
    }


    public void parseparameter() {

    }

    public void getsubRNAgenotype(String[] args) {
        String plate = args[0];
        String parameter = "/data2/xiaohan/genotype/hapscanner/output/" + args[0] + "/VCF";
        String RNADir = "/data2/xiaohan/genotype/hapscanner/output/" + args[0];

        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(String.valueOf(chr));
        }
        nameSet.parallelStream().forEach(f -> {
            BufferedReader br = null;
            if (Integer.parseInt(f) <= 9) {
                br = IOUtils.getTextReader(new File(parameter, "chr00" + f + ".vcf").getAbsolutePath());
            } else {
                br = IOUtils.getTextReader(new File(parameter, "chr0" + f + ".vcf").getAbsolutePath());
            }
            BufferedWriter bw = IOUtils.getTextWriter(new File(RNADir, "chr" + f + ".vcf").getAbsolutePath());
            String temp = null;
            String[] temps = null;
            try {
                while ((temp = br.readLine()) != null) {
                    if(temp.startsWith("#")){
                        bw.write(temp+"\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    int num = 0;
                    int total = temps.length/10;
                    for (int i = 9; i < temps.length; i++) {
                        if(temps[i].split(";")[0].equals("./.")){
                            num ++;
                        }
                    }
                    if(num <= total){
                        bw.write(temp+"\n");
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
            }
        });
    }


    public void getsubgenotype(String[] args) {
        String plate = args[0];
        String parameter = "/data2/xiaohan/genotype/hapscanner/output/" + args[0] + "/VCF";
        String DNADir = "/data2/xiaohan/genotype/hapscanner/output/" + args[0] + "DNA";
        BufferedReader brinfo = IOUtils.getTextReader(new File(parameter, "chr001.vcf").getAbsolutePath());
        BufferedWriter bwsample = IOUtils.getTextWriter(new File(DNADir, "samplelist.txt").getAbsolutePath());
        try {
            String temp = null;
            String[] temps = null;
            HashSet<String> nameSet = new HashSet();
            int num = 0;
            while ((temp = brinfo.readLine()) != null) {
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    for (int i = 9; i < temps.length; i++) {
                        if ((temps[i] + "").length() == 4) {
                            bwsample.write(temps[i] + "\n");
                        } else num++;
                    }
                }
            }
            brinfo.close();
            bwsample.flush();
            bwsample.close();
//
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                StringBuilder sb1 = new StringBuilder();
                sb1.append("bcftools view -S ");
                sb1.append(new File(DNADir, "samplelist.txt").getAbsolutePath());
                sb1.append(" /data2/xiaohan/genotype/genotype360/" + chr + ".346.B18.recode.vcf.gz -Ov > ");
                sb1.append(new File(DNADir, plate + "_oldchr" + chr + ".vcf").getAbsolutePath());
                sb1.append("  && bgzip ");
                sb1.append(new File(DNADir, plate + "_oldchr" + chr + ".vcf").getAbsolutePath());
                String command = sb1.toString();
                System.out.println(command);
                File dir = new File(new File(DNADir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }

            if (num != 0) {
                for (int i = 1; i < num; i++) {
                    for (int j = 0; j < 42; j++) {
                        int chr = j + 1;
                        StringBuilder sb1 = new StringBuilder();
                        sb1.append("bcftools view -S ");
                        sb1.append(new File("/data2/xiaohan/genotype/hapscanner/output/", "360.txt").getAbsolutePath());
                        sb1.append(" /data2/xiaohan/genotype/genotype360/" + chr + ".346.B18.recode.vcf.gz -Ov > ");
                        sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath());
                        sb1.append("sed -i \'s/" + "E360" + "/E360-" + num + "/g\' ");
                        sb1.append(sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath()));
                        sb1.append("bgzip ");
                        sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath());
                        String command = sb1.toString();
                        System.out.println(command);
                        File dir = new File(new File(DNADir).getAbsolutePath());
                        String[] cmdarry = {"/bin/bash", "-c", command};
                        Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                        p.waitFor();
                    }
                }
                for (int i = 0; i < 42; i++) {
                    int chr = i + 1;
                    StringBuilder sb2 = new StringBuilder();
                    sb2.append("vcf-concat ");
                    sb2.append(plate + "_oldchr" + chr + ".vcf.gz ");
                    for (int j = 1; j < num; j++) {
                        sb2.append(plate + "_oldchr" + chr + "_" + num + ".vcf.gz");
                    }
                    sb2.append(" > " + plate + "_chr" + chr + ".vcf \n");
                    sb2.append("bgzip " + plate + "_chr" + chr + ".vcf \n");
                    String command = sb2.toString();
                    System.out.println(command);
                    File dir = new File(new File(DNADir).getAbsolutePath());
                    String[] cmdarry = {"/bin/bash", "-c", command};
                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                    p.waitFor();
                }
            }

            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                StringBuilder sb2 = new StringBuilder();
                sb2.append("mv " + plate + "_oldchr" + chr + ".vcf.gz " + plate + "_chr" + chr + ".vcf.gz ");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getIsec(String[] args) {
        String plate = args[0];
        String infileDir1 = "/data2/xiaohan/genotype/hapscanner/output/" + plate;
        String infileDir2 = "/data2/xiaohan/genotype/hapscanner/output/" + plate + "DNA";
        File f = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec");
        f.mkdir();
        String outputDir = f.getAbsolutePath();
        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                String infile1 = new File(infileDir1, "chr" + chr + ".vcf").getAbsolutePath();
                String infile2 = new File(infileDir2, plate + "_chr" + chr + ".vcf.gz").getAbsolutePath();
                StringBuilder sb = new StringBuilder();
                sb.append("bedtools intersect -a " + infile1 + " -b " + infile2 + " -wa > " + plate + "_chr" + chr + "_RNA.vcf\n");
                sb.append("bedtools intersect -a " + infile2 + " -b " + infile1 + " -wa > " + plate + "_chr" + chr + "_DNA.vcf\n");
//                sb.append("rm "+infile1+"\n");
//                sb.append("rm "+infile2+"\n");
//                new File(infile1).delete();
//                new File(infile2).delete();
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File("/data2/xiaohan/genotype/hapscanner/output/").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getsortedVCF(String[] args) {
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        try {
            StringBuilder sb1 = new StringBuilder();
            sb1.append("cat " + plate + "_RNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + plate + "_RNA.sorted.vcf ");
            String command = sb1.toString();
            System.out.println(command);
            File dir = new File(infileDir).getAbsoluteFile();
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();

            StringBuilder sb2 = new StringBuilder();
            sb2.append("cat " + plate + "_DNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + plate + "_DNA.sorted.vcf ");
            command = sb2.toString();
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getMergedVCF(String[] args) {
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        try {
            StringBuilder sb1 = new StringBuilder();
            sb1.append("vcf-concat ");
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                sb1.append(plate + "_chr" + chr + "_RNA.vcf ");
            }
            sb1.append(" > " + plate + "_RNA.vcf");
            String command = sb1.toString();
            System.out.println(command);
            File dir = new File(infileDir).getAbsoluteFile();
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();

            StringBuilder sb2 = new StringBuilder();
            sb2.append("vcf-concat ");
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                sb2.append(plate + "_chr" + chr + "_DNA.vcf ");
            }
            sb2.append(" > " + plate + "_DNA.vcf");
            command = sb2.toString();
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void getIBdistane(String[] args) {
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        String infileS1 = new File(infileDir, plate + "_RNA.sorted.vcf").getAbsolutePath();
        String infileS2 = new File(infileDir, plate + "_DNA.sorted.vcf").getAbsolutePath();
//        String infileS1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/DNA.all.sort.vcf";
//        String infileS2 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/RNA.all.sort.vcf";
        String ibsOutfileS = new File(infileDir, "check.txt").getAbsolutePath();
//        String ibsOutfileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/check.txt";
        GenotypeGrid g1 = new GenotypeGrid(infileS1, GenoIOFormat.VCF);
        GenotypeGrid g2 = new GenotypeGrid(infileS2, GenoIOFormat.VCF);
        GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
        SumTaxaDivergence std = new SumTaxaDivergence(g);
        std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
        g.getIBSDistanceMatrix();
    }

    public void getDensityIBS(String[] args) {
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        String infile = new File(infileDir, "check.txt").getAbsolutePath();
        String outfile = new File(infileDir, "IBSdensity.txt").getAbsolutePath();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        int countlines = 0;
        try {
            bw.write("IBSdistance");
            bw.newLine();
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("E")) {
                    continue;
                }
                temps = temp.split("\t");
                if (countlines < 96) {
                    countlines++;
                    bw.write(temps[countlines + 96]);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void parameter(String[] args) {
        String plate = args[0];
        String outputdir = "/data2/xiaohan/genotype/hapscanner/";
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("mkdir " + plate + "&& mkdir" + plate + "DNA");
            String command = sb.toString();
            System.out.println(command);
            File dir = new File(new File("/data2/xiaohan/genotype/hapscanner/output/").getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, "parameter/" + plate + "_parameter_chr" + chr + ".txt").getAbsolutePath());
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
                bw.write("/data2/xiaohan/genotype/hapscanner/taxaRefBAM/" + plate + "_taxaRefBAM_chr" + chr + ".txt\n" +
                        "\n");
                bw.write("#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library. \n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n");
                bw.write("/data2/xiaohan/genotype/hapscanner/posAllele/posAllele_chr" + chr + ".txt\n" +
                        "\n");
                bw.write("#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n");
                bw.write("/data2/xiaohan/genotype/hapscanner/pos/pos_chr" + chr + ".txt\n" +
                        "\n");
                bw.write("#Parameter 4: The chromosome which will be scanned.\n" +
                        chr + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        "/usr/local/bin/samtools\n" +
                        "\n" +
                        "#Parameter 7: Number of threads\n" +
                        "8\n" +
                        "\n" +
                        "#Parameter 8: The directory of output\n");
                bw.write("/data2/xiaohan/genotype/hapscanner/output/" + plate + "\n");
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void taxaRefBAM(String[] args) {
        String parameter = "/data2/junxu/dataTest/test/" + args[0] + "/sams";
        File[] fs = new File(parameter).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, "_Aligned.out.sorted.bam");
        HashSet<String> nameSet = new HashSet();
        String plate = null;
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().replace("_Aligned.out.sorted.bam", ""));
        }
        try {
            String outfileDir = "/data2/xiaohan/genotype/hapscanner/";
            String[] namelist = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir, "taxaRefBAM/"+plate + "_taxaRefBAM_chr" + chr + ".txt").getAbsolutePath());
                bw.write("Taxa\tReference\tBamPath\n");
                for (int j = 0; j < namelist.length; j++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(namelist[j]).append("\t");
                    sb.append("/data2/xiaohan/genotype/hapscanner/ref/chr" + chr + ".fa\t");
                    sb.append(new File(parameter, namelist[j] + "_Aligned.out.sorted.bam").getAbsolutePath());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void poswithAllele() {
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            String infile = "/data2/xiaohan/genotype/genotype_eQTL/" + chr + ".346.B18.recode.vcf.gz";
            String outfile = "/data2/xiaohan/genotype/hapscanner/pos/pos_chr" + chr + ".txt";
            String outfile1 = "/data2/xiaohan/genotype/hapscanner/posAllele/posAllele_chr" + chr + ".txt";
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
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        new HapscannerParameters(args);
    }
}
