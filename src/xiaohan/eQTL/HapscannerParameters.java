package xiaohan.eQTL;

import com.sun.xml.internal.ws.api.message.ExceptionHasMessage;
import format.table.RowTable;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.IOFileFormat;
import smile.stat.Stat;
import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutionException;

public class HapscannerParameters {

    String parameterDir = "/data2/xiaohan/genotype/hapscanner/parameter/";
    String posDir = "/data2/xiaohan/genotype/hapscanner/pos/";
    String posAlleleDir = "/data2/xiaohan/genotype/hapscanner/posAllele/";
    String taxaRefBAMDir = "/data2/xiaohan/genotype/hapscanner/taxaRefBAM/";
    String outputDir = "/data2/xiaohan/genotype/hapscanner/output/";

//    String genotypeDir = "/data2/junxu/genotypeMaf005/";
//    String genotypesuffix = ".360.vcf.gz";
//    String DNAprefix = "";

    String genotypeDir = "/data2/junxu/genotype";
    String genotypesuffix = ".346.B18.recode.vcf.gz";
    String DNAprefix = "B18-";

    String BamDir = "/data2/junxu/dataTest/test/";
    String JarDir = "/data1/home/xiaohan/jar/";





    public HapscannerParameters(String[] args) throws IOException, InterruptedException {
//        this.parseparameter(args);
//        File posDir = new File(new File(this.posDir).getAbsolutePath());
//        File posAlleleDir = new File(new File(this.posAlleleDir).getAbsolutePath());
//        if (!posDir.exists() || !posAlleleDir.exists()) {
//        this.poswithAllele();
//        }
//        this.parameter(args);
//        this.taxaRefBAM(args);
//        this.Hapscanner(args);
//        this.Hapscanner();
//
//        this.changename(args);


//        this.getsubgenotype(args);
//        this.getsubRNAgenotype(args);

        this.getIsec(args);
        this.getMergedVCF(args);
        this.getsortedVCF(args);
        this.changeRNAVCFname(args);
        this.getIBdistane(args);
        this.getDensityIBS(args);
    }


    public void parseparameter(String[] args) {
        System.out.println("This is parsing parameter ******************************************************************");
        String plate = args[0];
        System.out.println("Dealing plate "+plate);
        String RNAdir = new File(this.outputDir, plate).getAbsolutePath();
        String DNAdir = new File(this.outputDir, plate + "/DNA").getAbsolutePath();
        String Isecdir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
        File f2 = new File(RNAdir);
        f2.mkdir();
        File f = new File(DNAdir);
        f.mkdir();
        File f3 = new File(Isecdir);
        f3.mkdir();
        try{
            if(!f.exists() || !f2.exists() || !f3.exists()){
                System.out.println("Not complete mkdir ***");
            }
        }
        catch (Exception e){

        }
    }

    public void getsubRNAgenotype(String[] args) {
        System.out.println("This is getting subRNA genotype ************************************************************");
        String plate = args[0];
        String parameter = this.outputDir + args[0] + "/changeName";
        String RNADir = this.outputDir + args[0];

        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(String.valueOf(chr));
        }

        BufferedWriter bw1 = IOUtils.getTextWriter(new File(RNADir, "header.txt").getAbsolutePath());
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
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        if (f.equals("1")) {
                            bw1.write(temp + "\n");
                        }
                        continue;
                    }
                    temps = temp.split("\t");
                    int num = 0;
                    int total = (int) (temps.length * 0.4);
                    for (int i = 9; i < temps.length; i++) {
                        if (temps[i].split(";")[0].equals("./.")) {
                            num++;
                        }
                    }
                    if (num <= total) {
                        bw.write(temp + "\n");
                    }
                }
                br.close();
                bw.flush();
                bw.close();
                bw1.flush();
                bw1.close();
            } catch (Exception e) {
            }
        });
    }

    public void changename(String[] args) {
        System.out.println("This is changing Names *********************************************************************");
        String plate = args[0];
        String parameter = this.outputDir + args[0] + "/VCF";
        String changeName = this.outputDir + args[0] + "/changeName";
        File f = new File(changeName);
        f.mkdir();
        File[] fs = new File(parameter).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName());
        }
        nameSet.parallelStream().forEach(c -> {
            try {
                StringBuilder sb = new StringBuilder();
                sb.append(" bcftools reheader -s /data2/xiaohan/genotype/hapscanner/output/" + plate + "_E.txt");
                sb.append(" " + new File(parameter, c).getAbsolutePath());
                sb.append(" -o " + new File(changeName, c).getAbsolutePath() + "\n");
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(parameter).getAbsoluteFile();
                String[] cmdarry1 = {"/bin/bash", "-c", command};
                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
                p1.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });


    }

//
//    public void getsubgenotype(String[] args) {
//        System.out.println("This is getting subGenotype ****************************************************************");
//        String plate = args[0];
//        String parameter = "/data2/xiaohan/genotype/hapscanner/output/" + args[0] + "/changeName";
//        String DNADir = "/data2/xiaohan/genotype/hapscanner/output/" + args[0] + "/DNA";
//        BufferedReader brinfo = IOUtils.getTextReader(new File(parameter, "chr001.vcf").getAbsolutePath());
//        BufferedWriter bwsample = IOUtils.getTextWriter(new File(DNADir, "samplelist.txt").getAbsolutePath());
//        HashSet<String> chrSet = new HashSet<>();
//        int num = 0;
//        for (int i = 0; i < 42; i++) {
//            int chr = i + 1;
//            chrSet.add(String.valueOf(chr));
//        }
//        try {
//            String temp = null;
//            String[] temps = null;
//
//            while ((temp = brinfo.readLine()) != null) {
//                if (temp.startsWith("#C")) {
//                    temps = temp.split("\t");
//                    for (int i = 9; i < temps.length; i++) {
//                        if ((temps[i] + "").length() == 4) {
//                            bwsample.write("B18-" + temps[i] + "\n");
//                        } else num++;
//                    }
//                }
//            }
//            brinfo.close();
//            bwsample.flush();
//            bwsample.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
////
//        chrSet.parallelStream().forEach(f -> {
//            try {
//                String chr = f;
//                StringBuilder sb1 = new StringBuilder();
//                sb1.append("bcftools view -S ");
//                sb1.append(new File(DNADir, "samplelist.txt").getAbsolutePath());
//                sb1.append(" /data2/junxu/genotype/" + chr + ".346.B18.recode.vcf.gz -Ov > ");
//                sb1.append(new File(DNADir, plate + "_oldchr" + chr + ".vcf").getAbsolutePath());
//                sb1.append("  && bgzip ");
//                sb1.append(new File(DNADir, plate + "_oldchr" + chr + ".vcf").getAbsolutePath());
//                String command = sb1.toString();
//                System.out.println(command);
//                File dir = new File(new File(DNADir).getAbsolutePath());
//                String[] cmdarry = {"/bin/bash", "-c", command};
//                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                p.waitFor();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//
//        try {
//            if (num != 0) {
//                for (int i = 1; i <= num; i++) {
//                    for (int j = 0; j < 42; j++) {
//                        int chr = i + 1;
//                        StringBuilder sb1 = new StringBuilder();
//                        sb1.append("bcftools view -S ");
//                        sb1.append(new File("/data2/xiaohan/genotype/hapscanner/output/", "360.txt").getAbsolutePath());
////                        sb1.append(" /data2/xiaohan/genotype/genotype360/" + chr + ".346.B18.recode.vcf.gz -Ov > ");
//                        sb1.append(" /data2/junxu/genotype/" + chr + ".346.B18.recode.vcf.gz -Ov > ");
//                        sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath());
//                        sb1.append("sed -i \'s/" + "E360" + "/E360-" + num + "/g\' ");
//                        sb1.append(sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath()));
//                        sb1.append("bgzip ");
//                        sb1.append(new File(DNADir, plate + "_oldchr" + chr + "_" + num + ".vcf\n").getAbsolutePath());
//                        String command = sb1.toString();
//                        System.out.println(command);
//                        File dir = new File(new File(DNADir).getAbsolutePath());
//                        String[] cmdarry = {"/bin/bash", "-c", command};
//                        Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                        p.waitFor();
//                    }
//                }
//                for (int i = 0; i < 42; i++) {
//                    int chr = i + 1;
//                    StringBuilder sb2 = new StringBuilder();
//                    sb2.append("vcf-merge ");
//                    sb2.append(plate + "_oldchr" + chr + ".vcf.gz ");
//                    for (int j = 1; j < num; j++) {
//                        sb2.append(plate + "_oldchr" + chr + "_" + num + ".vcf.gz");
//                    }
//                    sb2.append(" > " + plate + "_chr" + chr + ".vcf \n");
//                    sb2.append("bgzip " + plate + "_chr" + chr + ".vcf \n");
//                    String command = sb2.toString();
//                    System.out.println(command);
//                    File dir = new File(new File(DNADir).getAbsolutePath());
//                    String[] cmdarry = {"/bin/bash", "-c", command};
//                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                    p.waitFor();
//                }
//            } else {
//                for (int i = 0; i < 42; i++) {
//                    int chr = i + 1;
//                    StringBuilder sb2 = new StringBuilder();
//                    sb2.append("mv " + plate + "_oldchr" + chr + ".vcf.gz " + plate + "_chr" + chr + ".vcf.gz ");
//                    File dir = new File(new File(DNADir).getAbsolutePath());
//                    String command = sb2.toString();
//                    String[] cmdarry = {"/bin/bash", "-c", command};
//                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                    p.waitFor();
//                }
//            }
//        } catch (
//                Exception e) {
//            e.printStackTrace();
//        }
//
//    }

    public void getIsec(String[] args) {
        System.out.println("This is getting intersection ***************************************************************");
        String plate = args[0];
        String infileDir1 = new File(this.outputDir, plate).getAbsolutePath();
        String infileDir2 = new File(this.genotypeDir).getAbsolutePath();
        String outputDir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
        HashSet<String> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            chrSet.add(String.valueOf(chr));
        }

        BufferedReader br = IOUtils.getTextGzipReader(new File(this.genotypeDir,"1"+genotypesuffix).getAbsolutePath());
        BufferedWriter bw = IOUtils.getTextWriter(new File(this.outputDir,plate+"/DNAheader.txt").getAbsolutePath());
        String temp = null;
        try{
            while((temp = br.readLine())!=null){
                if(temp.startsWith("#")){
                    bw.write(temp+"\n");
                }
                if(!temp.startsWith("#"))break;
            }
            br.close();
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }
        chrSet.parallelStream().forEach(chr -> {
            try {
                String infile1 = new File(infileDir1, "chr" + chr + ".vcf").getAbsolutePath();
                String infile2 = new File(infileDir2, chr + genotypesuffix).getAbsolutePath();
                StringBuilder sb = new StringBuilder();
                sb.append("bedtools intersect -a " + infile1 + " -b " + infile2 + " -wa > " + plate + "_chr" + chr + "_RNA_noheader.vcf\n");
                sb.append("bedtools intersect -a " + infile2 + " -b " + infile1 + " -wa > " + plate + "_chr" + chr + "_DNA_noheader.vcf\n");

//                sb.append("rm "+infile1+"\n");
//                sb.append("rm "+infile2+"\n");
//                new File(infile1).delete();
//                new File(infile2).delete();
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File(outputDir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();

                sb.setLength(0);
                sb.append("cat ");
                sb.append(new File(infileDir1, "header.txt").getAbsolutePath());
                sb.append(" " + plate + "_chr" + chr + "_RNA_noheader.vcf > " + plate + "_chr" + chr + "_RNA.vcf\n");
                sb.append("cat ");
                sb.append(new File(infileDir1, "DNAheader.txt").getAbsolutePath());
                sb.append(" " + plate + "_chr" + chr + "_DNA_noheader.vcf > " + plate + "_chr" + chr + "_DNA.vcf\n");
                command = sb.toString();
                System.out.println(command);
                String[] cmdarry1 = {"/bin/bash", "-c", command};
                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
                p1.waitFor();

                sb.setLength(0);
                sb.append("rm " + plate + "_chr" + chr + "_RNA_noheader.vcf\n");
                sb.append("rm " + plate + "_chr" + chr + "_DNA_noheader.vcf\n");
                command = sb.toString();
                System.out.println(command);
                String[] cmdarry2 = {"/bin/bash", "-c", command};
                Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir);
                p2.waitFor();

            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void changeRNAVCFname(String[] args) {
        String plate = args[0];
        String infileDir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
        try {
            BufferedReader br = IOUtils.getTextReader(new File(infileDir, plate + "_RNA.sortedtemp.vcf").getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(infileDir, plate + "_RNA.sorted.vcf").getAbsolutePath());
            String temp = null;
            String[] temps = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    bw.write(temp + "\n");
                    continue;
                }
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < 9; i++) {
                        bw.write(temps[i] + "\t");
                    }
                    for (int i = 9; i < temps.length; i++) {
                        bw.write("RNA" + temps[i] + "\t");
                    }
                    bw.write("\n");
                    continue;
                }
                bw.write(temp + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getsortedVCF(String[] args) {
        System.out.println("This is sorting merged VCF files ***********************************************************");
        String plate = args[0];
        String infileDir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
        try {
            StringBuilder sb1 = new StringBuilder();
            sb1.append("cat " + plate + "_RNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + plate + "_RNA.sortedtemp.vcf");
            String command = sb1.toString();
            System.out.println(command);
            File dir = new File(infileDir).getAbsoluteFile();
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();

            StringBuilder sb2 = new StringBuilder();
            sb2.append("cat " + plate + "_DNA.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + plate + "_DNA.sorted.vcf");
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
        System.out.println("This is merging VCF files (RNA and DNA) ****************************************************");
        String plate = args[0];
        String infileDir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
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
        System.out.println("This is getting IBS matrix *****************************************************************");
        String plate = args[0];
        String infileDir = new File("/data2/xiaohan/genotype/hapscanner/output/", plate + "/Isec").getAbsolutePath();
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

//    public void getDensityIBS(String[] args) {
//        System.out.println("This is writing file of RNA and DNA IBS plot ***********************************************");
//        String plate = args[0];
//        String infileDir = new File(this.outputDir, plate + "Isec").getAbsolutePath();
//        String infile = new File(infileDir, "check.txt").getAbsolutePath();
//        String outfile = new File(infileDir, "IBSdensity.txt").getAbsolutePath();
////        String infile = "/data1/home/xiaohan/jar/check.txt";
////        String outfile = "/data1/home/xiaohan/jar/IBSdensity.txt";
//        BufferedReader br = IOUtils.getTextReader(infile);
//        BufferedWriter bw = IOUtils.getTextWriter(outfile);
//        String temp = null;
//        String[] temps = null;
//        int countlines = 0;
//        try {
//            bw.write("IBSdistance");
//            bw.newLine();
//            while ((temp = br.readLine()) != null) {
//                if (!temp.startsWith("E")) {
//                    continue;
//                }
//                temps = temp.split("\t");
//                if (countlines < 96) {
//                    countlines++;
//                    bw.write(temps[countlines + 96]);
//                    bw.newLine();
//                }
//            }
//            br.close();
//            bw.flush();
//            bw.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//    }


    public void getDensityIBS(String[] args) {
        System.out.println("This is writing file of RNA and DNA IBS plot ***********************************************");
        String plate = args[0];
        String infileDir = new File(this.outputDir, plate + "/Isec").getAbsolutePath();
//        String infor = new File(this.outputDir,plate+"/chr001.vcf").getAbsolutePath();
        String infile = new File(infileDir, "check.txt").getAbsolutePath();
        String outfile = new File(infileDir, "IBSdensity.txt").getAbsolutePath();
        String outfile1 = new File(infileDir,"IBSheatmap.txt").getAbsolutePath();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader brinfo = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
        String temp = null;
        String[] temps = null;
        int countlines = 0;
        try {
            HashMap<String, Integer> nameIndexMap = new HashMap<>();
            HashSet<String> samples = new HashSet<>();
            while ((temp = brinfo.readLine()) != null) {
                if (temp.startsWith("Dxy")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        nameIndexMap.put(temps[i], i);
                        if (temps[i].startsWith("RNA")) {
                            samples.add(temps[i]);
                        }
                    }
                }
            }
            brinfo.close();
            String[] samplelist = samples.toArray(new String[samples.size()]);

            bw.write("RNA\tDNA\tIBSdistance\n");
            xiaohan.eQTL.RowTable<String> t = new xiaohan.eQTL.RowTable<>(infile);
            for (int i = 0; i < samplelist.length; i++) {
                String RNA = samplelist[i];
                String DNA = this.DNAprefix+RNA.substring(3,7);
                int RNAindex = nameIndexMap.get(RNA);
                int DNAindex = nameIndexMap.get(DNA)-1;
                bw.write(RNA + "\t" + DNA + "\t");
                bw.write(t.getCell(DNAindex, RNAindex));
                bw.write("\n");
            }
            bw.flush();
            bw.close();

            bw1.write("RNA\tDNA\tIBSdistance\n");
            for (int i = 0; i < samplelist.length; i++) {
                String RNA = samplelist[i];
                for (int j = 0; j < samplelist.length; j++) {
                    String RNAtemp = samplelist[j];
                    String DNA = this.DNAprefix+RNAtemp.substring(3,7);
                    int RNAindex = nameIndexMap.get(RNA);
                    int DNAindex = nameIndexMap.get(DNA)-1;
                    bw1.write(RNA + "\t" + DNA + "\t");
                    bw1.write(t.getCell(DNAindex, RNAindex));
                    bw1.write("\n");
                }
            }
            bw1.flush();
            bw1.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void parameter(String[] args) {
        System.out.println("This is writing parameters files ***********************************************************");
        String plate = args[0];
        String outputdir = this.parameterDir;
        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, plate + "_parameter_chr" + chr + ".txt").getAbsolutePath());
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
                        "/usr/bin/samtools\n" +
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

    public void Hapscanner(String[] args) {
        System.out.println("This is running hapscanner *****************************************************************");
        String plate = args[0];
        BufferedWriter bw = IOUtils.getTextWriter("command_" + plate + ".txt");
        File[] fs = new File(String.valueOf(new File(this.parameterDir).getAbsoluteFile())).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, plate + "_");
        Arrays.sort(fs);
        try {
            for (int i = 0; i < fs.length; i++) {
                bw.write("java -jar TIGER_v1.0.1.jar -a HapScanner -p ");
                bw.write(String.valueOf(fs[i]));
                bw.write("\n");
            }
            bw.flush();
            bw.close();
            StringBuilder sb = new StringBuilder();
            sb.append("sh command1.sh command_" + plate + ".txt");
            String command = sb.toString();
            System.out.println(command);
            File dir = new File(new File(this.JarDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
//        fList.stream().forEach(p -> {
//            new HapScannercp(String.valueOf(p));
//        });
    }

//    public void Hapscanner(String[] args) {
//        System.out.println("This is running Hapscanner *****************************************************************");
////        List<File> fList = fList;
//        fList.stream().forEach(p -> {
//            try {
//                StringBuilder sb = new StringBuilder();
//                sb.append("java -jar /data1/home/xiaohan/jar/TIGER_v1.0.1.jar -a HapScaner -p ");
//                sb.append(p);
//                String command = sb.toString();
//                System.out.println(command);
//                File dir = new File("/data1/home/xiaohan/jar/").getAbsoluteFile();
//                String[] cmdarry1 = {"/bin/bash", "-c", command};
//                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
//                p1.waitFor();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//    }

    public void taxaRefBAM(String[] args) {
        System.out.println("This is writing taxaRefBam files *************************************************************************");
        String parameter = BamDir + args[0] + "/sams";
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
//
            String outfileDir = this.taxaRefBAMDir;
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
        File dir = new File(this.posDir).getAbsoluteFile();
        dir.mkdir();
        File dir2 = new File(this.posAlleleDir).getAbsoluteFile();
        dir2.mkdir();
        System.out.println("This is writing pos file ****************************************************");
        HashSet<String> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            chrSet.add(String.valueOf(chr));
        }
        chrSet.parallelStream().forEach(chr -> {
            System.out.println("Chr " + chr);
//            String infile = "/data2/xiaohan/genotype/genotype_eQTL/" + chr + ".346.B18.recode.vcf.gz";
            String infile = "/data2/junxu/genotype/" + chr + ".346.B18.recode.vcf.gz";
            String outfile = new File(this.posDir,"pos_chr" + chr + ".txt").getAbsolutePath();
            String outfile1 = new File(this.posAlleleDir,"posAllele_chr" + chr + ".txt").getAbsolutePath();
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
    }


    public static void main(String[] args) throws IOException, InterruptedException {
        new HapscannerParameters(args);
    }
}
