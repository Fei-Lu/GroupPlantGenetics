package xiaohan.eQTL;

import org.apache.commons.cli.*;
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
import java.util.concurrent.*;

public class Hapscann {

    Options options = new Options();

    //    String hapscanDir = "/data2/xiaohan/genotype/hapscanner";
    String hapscanDir = null;
    String outputDir = null;
    String parameterDir = null;
    String taxaRefBAMDir = null;
    String posDir = null;
    String posAlleleDir = null;

    String genotypeDir = null;
    String BamDir = null;
    String genotypesuffix = null;
    String bamsuffix = null;
    String jarDir = null;
    String samtools = null;
    String faDir = null;


    String plate = null;
    String threads = null;

    String introduction = this.createIntroduction();
    HelpFormatter optionFormat = new HelpFormatter();

    //    boolean changeName = false;
    String overwrite = "no";

    public Hapscann(String[] args) throws IOException, InterruptedException {
        this.createOptions();
        this.parseparameter(args);
        this.parameter();
        this.taxaRefBAM();
//        this.Hapscanner(args);
//        this.Hapscanner();
        this.Hapscanner1();
//
//        this.getsubRNAgenotypeor();
//        this.getIsec();
//        this.getMergedVCF();
//        this.getsortedVCF();
//        this.changeRNAVCFname();
//        this.getIBdistane();
//        this.getDensityIBS();
//        this.filtersample();
    }

    public void Hapscanner1(){
        File[] fs = new File(this.parameterDir).listFiles();
        fs = IOUtils.listFilesStartsWith(fs, plate + "_");
        Arrays.sort(fs);
        for (int i = 0; i < fs.length; i++) {
            String infile = fs[i].getAbsolutePath();
        }
    }

    static class Command implements Callable<Command> {
        String command = null;
        File dir = null;
        public Command (String command, File dir){
            this.command = command;
            this.dir = dir;
        }
        @Override
        public Command call() throws Exception {
            try {
                System.out.println(command);
                String[] cmdarry1 = {"/bin/bash", "-c", command};
                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
                p1.waitFor();
            }catch (Exception e){
                e.printStackTrace();
            }
            return this;
        }
    }

    public void filtersample() {
        long startTime = System.currentTimeMillis();   //获取开始时间
        System.out.println("This is filtering samples ******************************************************************");
        String outfileDir = new File(outputDir, "Isec").getAbsolutePath();
        String infile = new File(outfileDir, "IBSdensity.txt").getAbsolutePath();
        String infor = new File(outfileDir, "check.txt").getAbsolutePath();
        HashSet<String> sampleSet = new HashSet<>();
        HashSet<String> notsampleSet = new HashSet<>();
        HashSet<String> addingsampleSet = new HashSet<>();
        HashMap<String, String> RNADNAmap = new HashMap<>();
        try {
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getTextReader(infile);

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("RNA\t")) continue;
                temps = temp.split("\t");
                RNADNAmap.put(temps[0], temps[1]);
                if (Double.parseDouble(temps[2]) < 0.1) {
                    sampleSet.add(temps[0]);
                } else {
                    notsampleSet.add(temps[0]);
                }
            }
            br.close();

            System.out.println("original RNA samples");
            String[] temp1 = sampleSet.toArray(new String[sampleSet.size()]);
            for (int i = 0; i < temp1.length; i++) {
                System.out.println(temp1[i]);
            }

            BufferedReader brinfo = IOUtils.getTextReader(infor);


            xiaohan.eQTL.RowTable<String> t = new xiaohan.eQTL.RowTable<>(infor);
            //index of RNA samples and DNA samples : 1,2,3,...,RNAsamplecount,...,total
            HashMap<String, Integer> nameIndexMap = new HashMap<>();
            int RNAsamplecount = 0;
            int DNAsamplecount = 0;
            List<String> header = t.getHeader();
            String[] headers = header.toArray(new String[header.size()]);
            for (int i = 0; i < headers.length; i++) {
                if (headers[i].equals("Dxy")) continue;
                nameIndexMap.put(headers[i], i);
                if (headers[i].startsWith("RNA")) {
                    RNAsamplecount++;
                } else {
                    DNAsamplecount++;
                }
            }

            String[] notsamplelist = notsampleSet.toArray(new String[notsampleSet.size()]);

            for (int i = 0; i < notsamplelist.length; i++) {
                System.out.print("***This is dealing sample: " + notsamplelist[i] + " with original IBS : ");
                System.out.println(t.getCell(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(RNADNAmap.get(notsamplelist[i])) - 1));
                String sample = notsamplelist[i];
                int index = nameIndexMap.get(sample);

                double[] DNAIBS = new double[DNAsamplecount];
                for (int j = 0; j < DNAsamplecount; j++) {
                    DNAIBS[j] = Double.parseDouble(t.getCell(j + RNAsamplecount, index));
                }

                double min = Arrays.stream(DNAIBS).min().getAsDouble();

                if (min > 0.1) {
                    System.out.print("  Discard this sample: " + notsamplelist[i] + " with minimal IBS = ");
                    System.out.println(min);
                    continue;
                }

                System.out.print("  RNA sample has minimal IBS as ");
                System.out.println(min);

                HashSet<String> DNAcandidateSet = new HashSet<>();
                for (int j = 0; j < DNAIBS.length; j++) {
                    if (DNAIBS[j] == min) {
                        DNAcandidateSet.add(t.getCell(j + RNAsamplecount, 0));
                        System.out.println("    adding DNA candidate :" + t.getCell(j + RNAsamplecount, 0));
                    }
                }
                String[] DNAlist = DNAcandidateSet.toArray(new String[DNAcandidateSet.size()]);

                for (int j = 0; j < DNAlist.length; j++) {
                    int DNAindex = j + 1;
                    System.out.println("    " + DNAindex + ".This is examing DNA sample: " + DNAlist[j]);
                    if (j == 0) {
                        double[] RNAIBS = new double[RNAsamplecount];
                        for (int k = 0; k < RNAIBS.length; k++) {
                            RNAIBS[k] = Double.parseDouble(t.getCell(k, nameIndexMap.get(DNAlist[j])));
                        }

                        double min2 = Arrays.stream(RNAIBS).min().getAsDouble();
                        System.out.println("    DNAsample No." + DNAindex + " : " + DNAlist[j] + " has min IBS " + min2);

                        if (min == min2 && min < 0.1) {
                            System.out.println("    Two mins equal : adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
                            sampleSet.add(notsamplelist[i]);
                            addingsampleSet.add(notsamplelist[i]);
                            RNADNAmap.put(notsamplelist[i], DNAlist[j]);
                        }
                        if (min != min2) {
                            double residual = min - min2;
                            double abs = Math.abs(residual);
                            if (abs < 0.002) {
                                System.out.println("    Two mins not equal but abs of residual < 0.002 adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
                                sampleSet.add(notsamplelist[i]);
                                addingsampleSet.add(notsamplelist[i]);
                                RNADNAmap.put(notsamplelist[i], DNAlist[j]);
                            } else {
                                System.out.println("    Discard DNA sample for better RNA sample with IBS " + min2);
                            }
                        }
                        continue;
                    }

                    double IBStemp = t.getCellAsDouble(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(DNAlist[j]) - 1);
                    double IBSbefore = t.getCellAsDouble(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(RNADNAmap.get(notsamplelist[i])) - 1);
                    if (IBStemp == IBSbefore) {
                        sampleSet.remove(notsamplelist[i]);
                    } else if (IBSbefore < IBStemp) {
                        continue;
                    } else {
                        double[] RNAIBS = new double[RNAsamplecount];
                        for (int k = 0; k < RNAIBS.length; k++) {
                            RNAIBS[k] = Double.parseDouble(t.getCell(k, nameIndexMap.get(DNAlist[j])));
                        }

                        double min2 = Arrays.stream(RNAIBS).min().getAsDouble();
                        System.out.println("    DNAsample No." + DNAindex + " : " + DNAlist[j] + " has min IBS " + min2);

                        if (min == min2 && min < 0.1) {
                            System.out.println("    Two mins equal : adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
                            sampleSet.add(notsamplelist[i]);
                            RNADNAmap.put(notsamplelist[i], DNAlist[j]);
                        }
                        if (min != min2) {
                            double residual = min - min2;
                            double abs = Math.abs(residual);
                            if (abs < 0.002) {
                                System.out.println("    Two mins not equal but abs < 0.002 adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
                                sampleSet.add(notsamplelist[i]);
                                RNADNAmap.put(notsamplelist[i], DNAlist[j]);
                            } else {
                                System.out.println("    Discard DNA sample for better RNA sample with IBS " + min2);
                            }
                        }
                    }
                }
            }

            System.out.println("");
            System.out.println("Eventually we changed sample ... to new DNA sampe");
            for (String str : addingsampleSet) {
                System.out.println(str + "\t" + RNADNAmap.get(str));
            }


            BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir, "phenolist.txt").getAbsolutePath());
            String[] samplelist = sampleSet.toArray(new String[sampleSet.size()]);
            for (int i = 0; i < samplelist.length; i++) {
                bw.write(samplelist[i] + "\t" + RNADNAmap.get(samplelist[i]));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Filtering samples takes " + (endTime - startTime) + "ms");
        System.out.println("End of program.");
    }

    public void createOptions() {
        options = new Options();
//        options.addOption("g", true, "genotype Dir");
        options.addOption("b", true, "sorted bam file Dir");
        options.addOption("o", true, "output File Dir");
        options.addOption("p", true, "plate");
//        options.addOption("pos", true, "Whether or not to overwrite pos and posAllele Dir");
        options.addOption("t", true, "threads");
        options.addOption("j", true, "jar dir");
        options.addOption("bs", true, "bamsuffix");
//        options.addOption("gs", true, "genotypesuffix");
        options.addOption("samtools",true,"samtools path");
    }

    public void parseparameter(String[] args) {
        System.out.println("This is parsing parameter ******************************************************************");
        CommandLineParser parser = new DefaultParser();
        try {
            CommandLine line = parser.parse(options, args);
//            genotypeDir = line.getOptionValue("g");
            BamDir = line.getOptionValue("b");
            plate = line.getOptionValue("p");
            hapscanDir = line.getOptionValue("o");
//            overwrite = line.getOptionValue("pos");
            threads = line.getOptionValue("t");
            jarDir = line.getOptionValue("j");
            bamsuffix = line.getOptionValue("bs");
//            genotypesuffix = line.getOptionValue("gs");
            samtools = line.getOptionValue("samtools");
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(0);
        }


        bamsuffix = bamsuffix.replace(" ","");
        genotypesuffix = genotypesuffix.replace(" ","");
        System.out.println(bamsuffix);
        System.out.println(genotypesuffix);


        outputDir = hapscanDir + "/" + plate;

//        hapscanDir = outputDir.replace(outputDir.split("/")[outputDir.split("/").length-1],"");
        if (genotypeDir == null) {
            System.out.println("Genotype Dir doesn't exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (BamDir == null) {
            System.out.println("Bam Dir doesn't exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (hapscanDir == null) {
            System.out.println("output Dir doesn't exist");
            this.printIntroductionAndUsage();
            System.exit(0);
        }
        if (outputDir == null) {
            File f = new File(outputDir);
            f.mkdir();
        }

        System.out.println("Dealing plate " + plate);


        System.out.println("Start making dirs");
        taxaRefBAMDir = new File(hapscanDir, "taxaRefBAM").getAbsolutePath();
        File taxadir = new File(taxaRefBAMDir);
        if (!taxadir.exists()) taxadir.mkdir();
        System.out.println("mkdir taxaRefBAM");

        parameterDir = new File(hapscanDir, "parameter").getAbsolutePath();
        File paradir = new File(parameterDir);
        if (!paradir.exists()) paradir.mkdir();
        System.out.println("mkdir parameter");

        posDir = new File(hapscanDir, "pos").getAbsolutePath();
        posAlleleDir = new File(hapscanDir, "posAllele").getAbsolutePath();

        File posdir = new File(new File(hapscanDir, "pos").getAbsolutePath());
        File posAlleledir = new File(new File(hapscanDir, "posAllele").getAbsolutePath());

        if (overwrite.equals("yes") || overwrite.equals("y")) {
            posdir.mkdir();
            posAlleledir.mkdir();
            this.poswithAllele();
        } else if (overwrite.equals("no") || overwrite.equals("n")) {
            if (!posdir.exists() || !posAlleledir.exists()) {
                posdir.mkdir();
                posAlleledir.mkdir();
                this.poswithAllele();
            }
        }

//        String RNAdir = new File(this.outputDir, "RNA").getAbsolutePath();
//        String Isecdir = new File(this.outputDir, "Isec").getAbsolutePath();
//        File f2 = new File(RNAdir);
//        if (!f2.exists()) {
//            f2.mkdir();
//        }
//        File f3 = new File(Isecdir);
//        if (!f3.exists()) {
//            f3.mkdir();
//        }
//        try {
//            if (!f2.exists() || !f3.exists()) {
//                System.out.println("Not complete mkdir ***");
//            }
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//
    }

    public void printIntroductionAndUsage() {
        System.out.println("Incorrect options input. Program stops.");
        System.out.println(introduction);
        optionFormat.printHelp("Hapscan.jar", options);
    }

    public String createIntroduction() {
        StringBuilder sb = new StringBuilder();
        sb.append("Welcome to Use Hapscanner to genotyping and calculate IBS matrix.\n");
        sb.append("It uses six options to run this programs. : \n");
        sb.append("    -g  your genotype files dir is. \n");
        sb.append("    -b  your sorted bam files dir is. \n");
        sb.append("    -o  your output files dir is. \n");
        sb.append("    -p your output files dir prefix is. \n");
        sb.append("    -pos Whether or not to overwrite pos and posAllele files. e.g. \"yes\" or \"y\" for overwrite \"no\" or \"n\" for not. \n");
        sb.append("    -t number of threads to running this program. \n");
        sb.append("Command line like this: \n");
        sb.append("java -jar HapscannerParameter.jar -g ../genotype -b ../bams -o /data1/home/xiaohan/hapscanner -prefix 41 -pos n > log.txt &\n");
        sb.append("For more details, you can visit .\n");
        return sb.toString();
    }

    public void getsubRNAgenotypeor() {
        long startTime = System.currentTimeMillis();   //获取开始时间
        System.out.println("This is getting subRNA genotype ************************************************************");
        String RNADir = new File(outputDir, "RNA").getAbsolutePath();
        String originalRNA = new File(outputDir, "VCF").getAbsolutePath();
        HashSet<String> nameSet = new HashSet<>();
        File[] fs = new File(originalRNA).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName());
        }
        BufferedWriter bw1 = IOUtils.getTextWriter(new File(RNADir, "header.txt").getAbsolutePath());
        nameSet.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(new File(originalRNA, f).getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(RNADir, f).getAbsolutePath());
            String temp = null;
            String[] temps = null;
            try {
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        if (f.equals("chr001.vcf")) {
                            bw1.write(temp + "\n");
                        }
                        continue;
                    }
                    temps = temp.split("\t");
                    int num = 0;
//                    int total = (int) (temps.length * 0.4);
                    int total = temps.length;
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
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Getting subRNA genotype vcf takes " + (endTime - startTime) + "ms");
    }

//    public void changename(String[] args) {
//        System.out.println("This is changing Names *********************************************************************");
//        String plate = args[0];
//        String parameter = this.outputDir + args[0] + "/VCF";
//        String changeName = this.outputDir + args[0] + "/changeName";
//        File f = new File(changeName);
//        f.mkdir();
//        File[] fs = new File(parameter).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
//        HashSet<String> nameSet = new HashSet<>();
//        for (int i = 0; i < fs.length; i++) {
//            if (fs[i].isHidden()) continue;
//            nameSet.add(fs[i].getName());
//        }
//        nameSet.parallelStream().forEach(c -> {
//            try {
//                StringBuilder sb = new StringBuilder();
//                sb.append(" bcftools reheader -s /data2/xiaohan/genotype/hapscanner/output/" + plate + "_E.txt");
//                sb.append(" " + new File(parameter, c).getAbsolutePath());
//                sb.append(" -o " + new File(changeName, c).getAbsolutePath() + "\n");
//                String command = sb.toString();
//                System.out.println(command);
//                File dir = new File(parameter).getAbsoluteFile();
//                String[] cmdarry1 = {"/bin/bash", "-c", command};
//                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
//                p1.waitFor();
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//
//
//    }

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

    public void getIsec() {
        long startTime = System.currentTimeMillis();   //获取开始时间
        System.out.println("This is getting intersection ***************************************************************");
        String infileDir1 = new File(outputDir, "RNA").getAbsolutePath();
        String infileDir2 = new File(genotypeDir).getAbsolutePath();
        String IsecDir = new File(outputDir, "Isec").getAbsolutePath();
        HashSet<String> nameSet = new HashSet<>();
        File[] fs = new File(infileDir1).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName());
        }
        BufferedReader br = IOUtils.getTextGzipReader(new File(genotypeDir, "1" + genotypesuffix).getAbsolutePath());
        BufferedWriter bw = IOUtils.getTextWriter(new File(infileDir1, "DNAheader.txt").getAbsolutePath());
        String temp = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp + "\n");
                }
                if (!temp.startsWith("#")) break;
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        nameSet.parallelStream().forEach(f -> {
            try {
                int chr = Integer.parseInt(String.valueOf(f).replace(".vcf", "").replace("chr", ""));
                String infile1 = new File(infileDir1, f).getAbsolutePath();
                String infile2 = new File(infileDir2, chr + genotypesuffix).getAbsolutePath();
                StringBuilder sb = new StringBuilder();
                sb.append("bedtools intersect -a " + infile1 + " -b " + infile2 + " -wa > " + plate + "_chr" + chr + "_RNA_noheader.vcf\n");
                sb.append("bedtools intersect -a " + infile2 + " -b " + infile1 + " -wa > " + plate + "_chr" + chr + "_DNA_noheader.vcf\n");

                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File(outputDir, "Isec").getAbsolutePath());
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
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Getting intersect takes " + (endTime - startTime) + "ms");
    }

    public void changeRNAVCFname() {
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
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

    public void getsortedVCF() {
        System.out.println("This is sorting merged VCF files ***********************************************************");
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
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

    public void getMergedVCF() {
        System.out.println("This is merging VCF files (RNA and DNA) ****************************************************");
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
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


    public void getIBdistane() {
        System.out.println("This is getting IBS matrix *****************************************************************");
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
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


    public void getDensityIBS() {
        System.out.println("This is writing file of RNA and DNA IBS plot ***********************************************");
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
//        String infor = new File(this.outputDir,plate+"/chr001.vcf").getAbsolutePath();
        String infile = new File(infileDir, "check.txt").getAbsolutePath();
        String outfile = new File(infileDir, "IBSdensity.txt").getAbsolutePath();
        String outfile1 = new File(infileDir, "IBSheatmap.txt").getAbsolutePath();
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
                String DNA = null;
                if (RNA.contains("JM22")) {
                    DNA = "E025";
                } else {
                    DNA = RNA.substring(3, 7);
                }
                int RNAindex = nameIndexMap.get(RNA);
                int DNAindex = nameIndexMap.get(DNA) - 1;
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
                    String DNA = RNA.substring(3, 7);
                    int RNAindex = nameIndexMap.get(RNA);
                    int DNAindex = nameIndexMap.get(DNA) - 1;
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


    public void parameter() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing parameters files ***********************************************************");
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
                bw.write(new File(hapscanDir, "taxaRefBAM/" + plate + "_taxaRefBAM_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library. \n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n");
                bw.write(new File(hapscanDir, "posAllele/posAllele_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n");
                bw.write(new File(hapscanDir, "pos/pos_chr" + chr + ".txt\n").getAbsolutePath() +
                        "\n");
                bw.write("#Parameter 4: The chromosome which will be scanned.\n" +
                        chr + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        new File(samtools).getAbsolutePath()+"\n" +
                        "\n" +
                        "#Parameter 7: Number of threads\n" +
                        threads + "\n" +
                        "\n" +
                        "#Parameter 8: The directory of output\n");
                bw.write(outputDir + "\n");
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis();
        System.out.println("******* Writing parameter files takes " + (endTime - startTime) + "ms");
    }

//    public void Hapscanner() {
//        File[] fs = new File(parameterDir).listFiles();
//        fs = IOUtils.listFilesStartsWith(fs,plate);
//        for (int i = 0; i < fs.length; i++) {
//            System.out.println(fs[i].getAbsolutePath());
//            new HapScannercp(fs[i].getAbsolutePath());
//        }
//    }

    public void Hapscanner() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is running hapscanner *****************************************************************");
        String JarDir = jarDir;
        BufferedWriter bw = IOUtils.getTextWriter("command_" + plate + ".txt");
        File[] fs = new File(this.parameterDir).listFiles();
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
            sb.append("java -jar ThreadPool.jar command_" + plate + ".txt 1 ");
            String command = sb.toString();
            System.out.println(command);
            File dir = new File(new File(JarDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Hapscanner genotyping takes " + (endTime - startTime) + "ms");
    }

    public List<String> listfile(String file) {
        List<String> files = new ArrayList<>();
        File[] fileDir = null;
        if (new File(file).isDirectory()) {
            fileDir = new File(file).listFiles();
            for (int i = 0; i < fileDir.length; i++) {
                if (fileDir[i].isDirectory()) {
                    listfile(new File(String.valueOf(fileDir[i])).getAbsolutePath());
                } else {
                    files.add(fileDir[i].getAbsolutePath());
                }
            }
        } else {
            files.add(file);
        }
        return files;
    }

    public void taxaRefBAM() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing taxaRefBam files *************************************************************************");
        ArrayList<String> files = new ArrayList<>();
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
//                    int length = namelist[j].split("/").length;
//                    sb.append(namelist[j].split("/")[length - 1].replace(bamsuffix, "")).append("\t");
                    sb.append(namelist[j].split("/")[namelist[j].split("/").length-1].replace(bamsuffix, "")).append("\t");
                    sb.append("/data1/home/xiaohan/hapscanner/ref/chr" + chr + ".fa\t");
                    sb.append(new File(namelist[j]).getAbsolutePath());
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


    public static void main(String[] args) throws IOException, InterruptedException {
        new Hapscann(args);
    }
}
