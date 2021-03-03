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
        this.parseparameter();
        File posDir = new File(new File("/data2/xiaohan/genotype/hapscanner/pos").getAbsolutePath());
        File posAlleleDir = new File(new File("/data2/xiaohan/genotype/hapscanner/posAllele").getAbsolutePath());
        if (!posDir.exists() || !posAlleleDir.exists()) {
            this.poswithAllele();
        }
        this.parameter(args);
        this.taxaRefBAM(args);
        this.getsubgenotype(args);

        this.getSamePosition(args);
        this.getIsec(args);
        this.getMergedVCF(args);
        this.getsortedVCF(args);
        this.getIBdistane(args);
        this.getDensityIBS(args);
    }


    public void parseparameter() {

    }

    public void getsubgenotype(String[] args) {
        String parameter = "/data2/junxu/dataTest/test/" + args[0] + "/sams";
        File[] fs = new File(parameter).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, "_Aligned.out.sorted.bam");
        HashSet<String> nameSet = new HashSet();
        String plate = args[0];
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().replace("_" + plate + "_Aligned.out.sorted.bam", ""));
        }
        String[] namelist = nameSet.toArray(new String[nameSet.size()]);
        Arrays.sort(namelist);

        String genotypeDir = "/data2/xiaohan/genotype/hapscanner/output/" + plate + "DNA";
        StringBuilder sb = new StringBuilder();
        try {
            for (int i = 0; i < namelist.length; i++) {
                sb.append(namelist[i] + "\n");
            }
            BufferedWriter bw = IOUtils.getTextWriter(new File(genotypeDir, "samplelist.txt").getAbsolutePath());
            bw.write(sb.toString());
            bw.flush();
            bw.close();

            BufferedWriter bw1 = IOUtils.getTextWriter(new File(genotypeDir, "substract&merge.txt").getAbsolutePath());
            for (int i = 0; i < 42; i++) {
                int chr = i+1;
                StringBuilder sb1 = new StringBuilder();
                sb1.append("bcftools view -S ");
                sb1.append(new File(genotypeDir, "samplelist.txt").getAbsolutePath());
                sb1.append(" /data2/xiaohan/genotype/genotype_eQTL/"+chr+".346.B18.recode.vcf.gz -Ov > ");
                sb1.append(new File(genotypeDir,plate+"_chr"+chr+".vcf").getAbsolutePath());
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File(genotypeDir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getIsec(String[] args){
        String plate = args[0];
        String infileDir1 = "/data2/xiaohan/genotype/hapscanner/output/"+plate;
        String infileDir2 = "/data2/xiaohan/genotype/hapscanner/output/"+plate + "DNA";
        File f = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec");
        f.mkdir();
        String outputDir = f.getAbsolutePath();
        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                String infile1 = new File(infileDir1, "").getAbsolutePath();
                String infile2 = new File(infileDir2, "").getAbsolutePath();
                StringBuilder sb = new StringBuilder();
                sb.append("bedtools intersect -a "+infile1+" -b "+infile2+" -wa > chr"+chr+"_A.vcf\n");
                sb.append("bedtools intersect -a "+infile2+" -b "+infile1+" -wa > chr"+chr+"_B.vcf\n");
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File("/data2/xiaohan/genotype/hapscanner/output/").getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void getsortedVCF(String[] args){
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        try {
            StringBuilder sb1 = new StringBuilder();
            sb1.append("cat RNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > RNA.sorted.vcf ");
            String command = sb1.toString();
            System.out.println(command);
            File dir = new File(infileDir).getAbsoluteFile();
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();

            StringBuilder sb2 = new StringBuilder();
            sb2.append("cat DNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > DNA.sorted.vcf ");
            command = sb2.toString();
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void getMergedVCF(String[] args){
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        try {
            StringBuilder sb1 = new StringBuilder();
            sb1.append("vcf-concat ");
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                sb1.append(chr + ".vcf ");
            }
            sb1.append(" > RNA.vcf");
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
                sb2.append(chr + ".vcf ");
            }
            sb2.append(" > DNA.vcf");
            command = sb2.toString();
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    public void getSamePosition(String[] args) {
        String plate = args[0];
        String inputDirRNA = "/data1/home/xiaohan/rareallele/Hapscanner/outputDir/mis01VCF";
        String inputDirDNA = "/data2/xiaohan/DNAgenotype";
        String outputDirRNA = "/data2/xiaohan/RNAgenotype";
        String outputDirDNA = "/data2/xiaohan/DNAgenotype";
//        String inputDirRNA = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/changename/inputDirRNA";
//        String inputDirDNA = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/changename/inputDirDNA";
//        String outputDirRNA = "/Users/yxh/Documents/RareAllele/004te..st/SiPASpipeline/input/changename/outputDirRNA";
//        String outputDirDNA = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/changename/outputDirDNA";
        String tempRNA = null;
        String[] tempRNAs = null;
        String tempDNA = null;
        String[] tempDNAs = null;
        String tempRNA1 = null;
        String[] tempRNAs1 = null;
        String tempDNA1 = null;
        String[] tempDNAs1 = null;
        for (int i = 0; i < 42; i++) {
            int chrNumber = i + 1;
//            int chrNumber = 1;
            try {
                //
//                BufferedReader brRNA = IOUtils.getTextReader(new File(inputDirRNA, "chr" + chrNumber + ".mis01.recode.vcf").getAbsolutePath());
//                HashSet<String> PositionSetRNA = new HashSet();
                HashSet<String> PositionSetDNA = new HashSet();
                HashMap<String, String> positionRNA = new HashMap();
                HashMap<String, String> positionDNA = new HashMap();
//                while ((tempRNA = brRNA.readLine()) != null) {
//                    if (tempRNA.startsWith("#")) continue;
//                    tempRNAs = tempRNA.split("\t");
//                    String position1 = tempRNAs[1];
//                    PositionSetRNA.add(position1);
//                    positionRNA.put(position1, tempRNA);
//                }
//                brRNA.close();
                BufferedReader brDNA = IOUtils.getTextReader(new File(inputDirDNA, "DNAchr" + chrNumber + ".vcf").getAbsolutePath());
                while ((tempDNA = brDNA.readLine()) != null) {
                    if (tempDNA.startsWith("#")) {
                        continue;
                    }
                    tempDNAs = tempDNA.split("\t");
                    String position2 = tempDNAs[1];
                    PositionSetDNA.add(position2);
                    positionDNA.put(position2, tempDNA);
                }
                brDNA.close();
//                //
                BufferedReader brRNA1 = IOUtils.getTextReader(new File(inputDirRNA, "chr" + chrNumber + ".mis01.recode.vcf").getAbsolutePath());
                BufferedWriter bwRNA = IOUtils.getTextWriter(new File(outputDirRNA, "RNAchr" + chrNumber + ".vcf").getAbsolutePath());
                int count = 0;
                System.out.println("This is writing chr" + chrNumber + " File……………………………………………………………………………………………………………………");
                while ((tempRNA1 = brRNA1.readLine()) != null) {
                    if (tempRNA1.startsWith("#")) {
                        bwRNA.write(tempRNA1);
                        bwRNA.newLine();
                        continue;
                    }
                    tempRNAs1 = tempRNA1.split("\t");
                    String positionR = tempRNAs1[1];
                    if (PositionSetDNA.contains(positionR)) {
                        bwRNA.write(tempRNA1);
                        bwRNA.newLine();
                    }
                    if (count % 5000 == 1) {
                        System.out.println(count);
                    }
                }
                brRNA1.close();
                bwRNA.flush();
                bwRNA.close();
//                BufferedReader brDNA1 = IOUtils.getTextGzipReader(new File(inputDirDNA, chrNumber + ".92.B18.recode.vcf.gz").getAbsolutePath());
//                BufferedWriter bwDNA = IOUtils.getTextWriter(new File(outputDirDNA, "DNAchr" + chrNumber + ".vcf").getAbsolutePath());
//                int count = 0;
//                System.out.println("This is writing chr"+chrNumber+" File……………………………………………………………………………………………………………………");
//                while ((tempDNA1 = brDNA1.readLine()) != null) {
//                    count ++;
//                    if (tempDNA1.startsWith("#")) {
//                        bwDNA.write(tempDNA1);
//                        bwDNA.newLine();
//                        continue;
//                    }
//                    tempDNAs1 = tempDNA1.split("\t");
//                    String positionD = tempDNAs1[1];
//                    if (PositionSetRNA.contains(positionD)) {
//                        bwDNA.write(tempDNA1);
//                        bwDNA.newLine();
//                    }
//                    if(count % 5000 ==1){
//                        System.out.println(count);
//                    }
//                }
//                brDNA1.close();
//                bwDNA.flush();bwDNA.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }

    }

    public void getIBdistane(String[] args) {
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "Isec").getAbsolutePath();
        String infileS1 = new File(infileDir,"RNA.sorted.vcf").getAbsolutePath();
        String infileS2 = new File(infileDir,"DNA.sorted.vcf").getAbsolutePath();
//        String infileS1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/DNA.all.sort.vcf";
//        String infileS2 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/RNA.all.sort.vcf";
        String ibsOutfileS = new File(infileDir,"check.txt").getAbsolutePath();
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
        String infile = new File(infileDir,"check.txt").getAbsolutePath();
        String outfile = new File(infileDir,"IBSdensity.txt").getAbsolutePath();
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
                if (countlines < 92) {
                    countlines++;
                    bw.write(temps[countlines + 92]);
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
        String outputdir = args[1];
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
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, "parameters" + plate + "_parameter_chr" + chr + ".txt").getAbsolutePath());
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
            String outfileDir = args[1];
            String[] namelist = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir, plate + "_taxaRefBAM_chr" + chr + ".txt").getAbsolutePath());
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
