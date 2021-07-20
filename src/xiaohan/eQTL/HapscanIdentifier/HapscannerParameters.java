package xiaohan.eQTL.HapscanIdentifier;

import com.koloboke.collect.map.IntDoubleMap;
import com.koloboke.collect.map.hash.HashIntDoubleMaps;
import org.apache.commons.cli.*;
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

public class HapscannerParameters {

    Options options = new Options();

    String hapscanDir = null;
    String outputDir = null;
    String parameterDir = null;
    String taxaRefBAMDir = null;
    String posDir = null;
    String posAlleleDir = null;
    String referenceDir = null;

    String genotypeDir = null;
    String BamDir = null;
    String genotypesuffix = ".360.vcf.gz";

    int chrNumber = -1;

    String plate = null;
    String threads = null;

    String introduction = this.createIntroduction();
    HelpFormatter optionFormat = new HelpFormatter();

    //    boolean changeName = false;
    String overwrite = "no";

    //hapscan
    //The taxaRefBam file containing information of taxon and its corresponding bam files. The bam file should have .bai file in the same folder
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
    private BufferedReader br1;

    public HapscannerParameters(String[] args) throws IOException, InterruptedException {
//        this.createOptions();
//        this.parseparameter(args);
        this.parseparameters(args[0]);
        chrNumber = Integer.parseInt(args[1]);
//        this.mkdir(args);
        this.parameter();
        this.taxaRefBAM();
        this.Hapscanner();

        this.getsubRNAgenotypeor1();
        this.getIsec1();
        this.getMerge();
//        this.getMergedVCF();
//        this.getsortedVCF();
//        this.changeRNAVCFname();
        this.getIBdistane();
        this.getDensityIBS();
        this.filtersample();
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


            RowTable<String> t = new RowTable<>(infor);
            //index of RNA samples and DNA samples : 1,2,3,...,RNAsamplecount,...,total
            String[] header = t.getHeader().toArray(new String[0]);
            List<String> RNASet = new ArrayList<>();
            List<String> DNASet = new ArrayList<>();
            for (int i = 0; i < header.length; i++) {
                if (header[i].startsWith("RNA")) {
                    RNASet.add(header[i]);
                } else if (!header[i].startsWith("Dxy")) {
                    DNASet.add(header[i]);
                }
            }

            HashMap<String, Integer> nameIndexMap = new HashMap<>();
            BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir, "correction.txt").getAbsolutePath());
            BufferedWriter bw1 = IOUtils.getTextWriter(new File(outfileDir, "list1.txt").getAbsolutePath());
            BufferedWriter bw2 = IOUtils.getTextWriter(new File(outfileDir, "list2.txt").getAbsolutePath());
            BufferedWriter bw3 = IOUtils.getTextWriter(new File(outfileDir, "list3.txt").getAbsolutePath());


            for (int i = 0; i < RNASet.size(); i++) {
                String RNA = RNASet.get(i);
                String DNA = null;
                if (RNA.contains("JM22")) {
                    DNA = "E025";
                } else if (RNA.contains("CS")) {
                    DNA = "E360";
                } else {
                    DNA = RNA.substring(3, 7);
                }

                double[] IBS = t.getColumnAsDoubleArray(t.getColumnIndex(RNA));
                double[] subIBS = Arrays.copyOfRange(IBS, RNASet.size(), header.length - 1);
                // 对应值
                double IBSvalue = Double.parseDouble(t.getCell(t.getColumnIndex(RNA) - 1, t.getColumnIndex(DNA)));
                double[] Min = SNPmappingInGene.Min(subIBS);

                // 最小的DNA
                String DNATrue = header[RNASet.size() + (int) Min[1] + 1];
                // 最小的IBS
                double IBSDNAvalue = Min[0];

                bw.write(RNA + "\t" + DNA + "\t" + IBSvalue + "\t" + DNATrue + "\t" + IBSDNAvalue + "\t");

                if (DNA.equals(DNATrue)) {
                    bw.write("TRUE" + "\n");
                    bw1.write(RNA + "\t" + DNA + "\n");
                } else {
                    bw.write("False" + "\n");
                    if (Math.abs(IBSDNAvalue - IBSvalue) < 0.01) {
                        bw2.write(RNA + "\t" + DNA + "\n");
                    } else if (IBSDNAvalue < 0.1) {
                        bw3.write(RNA + "\t" + DNATrue + "\n");
                    }
                }

//                if (IBSvalue < 0.1) {
//                    bw.write(RNA + "\t" + DNA + "\n");
//                } else if (DNA.equals(DNATrue)) {
//                    bw.write(RNA + "\t" + DNA + "\t" + IBSDNAvalue + "\n");
//                } else {
//                    double[] IBSDNA = t.getColumnAsDoubleArray(t.getColumnIndex(DNA));
//                    double[] subIBSDNA = Arrays.copyOfRange(IBS, 0, RNASet.size());
//                    double[] Min2 = SNPmappingInGene.Min(subIBSDNA);
//                    // 最小的DNA对应的最小的RNA
//                    String RNATrue = header[RNASet.size() + (int) Min[1] + 1];
//                    // 最小的IBS
//                    double IBSRNAvalue = Min[0];
//                    if (RNA.equals(RNATrue)) {
//                        bw.write(RNA + "\t" + DNATrue + "\n");
//                    }
//                }
            }
//            int RNAsamplecount = 0;
//            int DNAsamplecount = 0;
//            List<String> header = t.getHeader();
//            String[] headers = header.toArray(new String[header.size()]);
//            for (int i = 0; i < headers.length; i++) {
//                if (headers[i].equals("Dxy")) continue;
//                nameIndexMap.put(headers[i], i);
//                if (headers[i].startsWith("RNA")) {
//                    RNAsamplecount++;
//                } else {
//                    DNAsamplecount++;
//                }
//            }
//
//            String[] notsamplelist = notsampleSet.toArray(new String[notsampleSet.size()]);
//
//            for (int i = 0; i < notsamplelist.length; i++) {
//                System.out.print("***This is dealing sample: " + notsamplelist[i] + " with original IBS : ");
//                System.out.println(t.getCell(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(RNADNAmap.get(notsamplelist[i])) - 1));
//                String sample = notsamplelist[i];
//                int index = nameIndexMap.get(sample);
//
//                double[] DNAIBS = new double[DNAsamplecount];
//                for (int j = 0; j < DNAsamplecount; j++) {
//                    DNAIBS[j] = Double.parseDouble(t.getCell(j + RNAsamplecount, index));
//                }
//
//                double min = Arrays.stream(DNAIBS).min().getAsDouble();
//
//                if (min > 0.1) {
//                    System.out.print("  Discard this sample: " + notsamplelist[i] + " with minimal IBS = ");
//                    System.out.println(min);
//                    continue;
//                }
//
//                System.out.print("  RNA sample has minimal IBS as ");
//                System.out.println(min);
//
//                HashSet<String> DNAcandidateSet = new HashSet<>();
//                for (int j = 0; j < DNAIBS.length; j++) {
//                    if (DNAIBS[j] == min) {
//                        DNAcandidateSet.add(t.getCell(j + RNAsamplecount, 0));
//                        System.out.println("    adding DNA candidate :" + t.getCell(j + RNAsamplecount, 0));
//                    }
//                }
//                String[] DNAlist = DNAcandidateSet.toArray(new String[DNAcandidateSet.size()]);
//
//                for (int j = 0; j < DNAlist.length; j++) {
//                    int DNAindex = j + 1;
//                    System.out.println("    " + DNAindex + ".This is examing DNA sample: " + DNAlist[j]);
//                    if (j == 0) {
//                        double[] RNAIBS = new double[RNAsamplecount];
//                        for (int k = 0; k < RNAIBS.length; k++) {
//                            RNAIBS[k] = Double.parseDouble(t.getCell(k, nameIndexMap.get(DNAlist[j])));
//                        }
//
//                        double min2 = Arrays.stream(RNAIBS).min().getAsDouble();
//                        System.out.println("    DNAsample No." + DNAindex + " : " + DNAlist[j] + " has min IBS " + min2);
//
//                        if (min == min2 && min < 0.1) {
//                            System.out.println("    Two mins equal : adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
//                            sampleSet.add(notsamplelist[i]);
//                            addingsampleSet.add(notsamplelist[i]);
//                            RNADNAmap.put(notsamplelist[i], DNAlist[j]);
//                        }
//                        if (min != min2) {
//                            double residual = min - min2;
//                            double abs = Math.abs(residual);
//                            if (abs < 0.002) {
//                                System.out.println("    Two mins not equal but abs of residual < 0.002 adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
//                                sampleSet.add(notsamplelist[i]);
//                                addingsampleSet.add(notsamplelist[i]);
//                                RNADNAmap.put(notsamplelist[i], DNAlist[j]);
//                            } else {
//                                System.out.println("    Discard DNA sample for better RNA sample with IBS " + min2);
//                            }
//                        }
//                        continue;
//                    }
//
//                    double IBStemp = t.getCellAsDouble(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(DNAlist[j]) - 1);
//                    double IBSbefore = t.getCellAsDouble(nameIndexMap.get(notsamplelist[i]), nameIndexMap.get(RNADNAmap.get(notsamplelist[i])) - 1);
//                    if (IBStemp == IBSbefore) {
//                        sampleSet.remove(notsamplelist[i]);
//                    } else if (IBSbefore < IBStemp) {
//                        continue;
//                    } else {
//                        double[] RNAIBS = new double[RNAsamplecount];
//                        for (int k = 0; k < RNAIBS.length; k++) {
//                            RNAIBS[k] = Double.parseDouble(t.getCell(k, nameIndexMap.get(DNAlist[j])));
//                        }
//
//                        double min2 = Arrays.stream(RNAIBS).min().getAsDouble();
//                        System.out.println("    DNAsample No." + DNAindex + " : " + DNAlist[j] + " has min IBS " + min2);
//
//                        if (min == min2 && min < 0.1) {
//                            System.out.println("    Two mins equal : adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
//                            sampleSet.add(notsamplelist[i]);
//                            RNADNAmap.put(notsamplelist[i], DNAlist[j]);
//                        }
//                        if (min != min2) {
//                            double residual = min - min2;
//                            double abs = Math.abs(residual);
//                            if (abs < 0.002) {
//                                System.out.println("    Two mins not equal but abs < 0.002 adding sample " + notsamplelist[i] + " and replace DNA sample " + RNADNAmap.get(notsamplelist[i]) + " as " + DNAlist[j]);
//                                sampleSet.add(notsamplelist[i]);
//                                RNADNAmap.put(notsamplelist[i], DNAlist[j]);
//                            } else {
//                                System.out.println("    Discard DNA sample for better RNA sample with IBS " + min2);
//                            }
//                        }
//                    }
//                }
//            }
//
//            System.out.println("");
//            System.out.println("Eventually we changed sample ... to new DNA sampe");
//            for (String str : addingsampleSet) {
//                System.out.println(str + "\t" + RNADNAmap.get(str));
//            }
//
//
//            String[] samplelist = sampleSet.toArray(new String[sampleSet.size()]);
//            for (int i = 0; i < samplelist.length; i++) {
//                bw.write(samplelist[i] + "\t" + RNADNAmap.get(samplelist[i]));
//                bw.newLine();
//            }
            bw.flush();
            bw.close();
            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();
            bw3.flush();
            bw3.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Filtering samples takes " + (endTime - startTime) + "ms");
        System.out.println("End of program.");
    }

//    public void createOptions() {
//        options = new Options();
//        options.addOption("g", true, "genotype Dir");
//        options.addOption("b", true, "sorted bam file Dir");
//        options.addOption("o", true, "output File Dir");
//        options.addOption("p", true, "plate");
//        options.addOption("pos", true, "Whether or not to overwrite pos and posAllele Dir");
//        options.addOption("t", true, "threads");
//    }

//    public void parseparameter(String[] args) {
//        System.out.println("This is parsing parameter ******************************************************************");
//        CommandLineParser parser = new DefaultParser();
//        try {
//            CommandLine line = parser.parse(options, args);
//            genotypeDir = line.getOptionValue("g");
//            BamDir = line.getOptionValue("b");
//            plate = line.getOptionValue("p");
//            outputDir = line.getOptionValue("o");
//            overwrite = line.getOptionValue("pos");
//            threads = line.getOptionValue("t");
//        } catch (Exception e) {
//            e.printStackTrace();
//            System.exit(0);
//        }
//        if (genotypeDir == null) {
//            System.out.println("Genotype Dir doesn't exist");
//            this.printIntroductionAndUsage();
//            System.exit(0);
//        }
//        if (BamDir == null) {
//            System.out.println("Bam Dir doesn't exist");
//            this.printIntroductionAndUsage();
//            System.exit(0);
//        }
//        if (outputDir == null) {
//            System.out.println("Output Dir doesn't exist");
//            this.printIntroductionAndUsage();
//            System.exit(0);
//        }
//        System.out.println("Dealing plate " + plate);
//
//        System.out.println("Start making dirs");
//        taxaRefBAMDir = new File(hapscanDir, "taxaRefBAM").getAbsolutePath();
//        File taxadir = new File(taxaRefBAMDir);
//        if (!taxadir.exists()) taxadir.mkdir();
//        System.out.println("mkdir taxaRefBAM");
//
//        parameterDir = new File(hapscanDir, "parameter").getAbsolutePath();
//        File paradir = new File(parameterDir);
//        if (!paradir.exists()) paradir.mkdir();
//        System.out.println("mkdir parameter");
//
//        posDir = new File(hapscanDir, "pos").getAbsolutePath();
//        posAlleleDir = new File(hapscanDir, "posAllele").getAbsolutePath();
//
//        File posdir = new File(new File(hapscanDir, "pos").getAbsolutePath());
//        File posAlleledir = new File(new File(hapscanDir, "posAllele").getAbsolutePath());
//
//        if (overwrite.equals("yes") || overwrite.equals("y")) {
//            posdir.mkdir();
//            posAlleledir.mkdir();
//            this.poswithAllele();
//        } else if (overwrite.equals("no") || overwrite.equals("n")) {
//            if (!posdir.exists() || !posAlleledir.exists()) {
//                posdir.mkdir();
//                posAlleledir.mkdir();
//                this.poswithAllele();
//            }
//        }
//
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
//    }

    public void getsubRNAgenotypeor1() {
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
        nameSet.parallelStream().forEach(f -> {
            BufferedReader br = IOUtils.getTextReader(new File(originalRNA, f).getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(RNADir, f.replace("chr", "RNA_chr")).getAbsolutePath());
            String temp = null;
            String[] temps = null;
            StringBuilder sb = new StringBuilder();
            try {
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("##")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    if (temp.startsWith("#C")) {
                        for (int i = 0; i < 9; i++) {
                            sb.append(temps[i] + "\t");
                        }
                        for (int i = 9; i < temps.length; i++) {
                            sb.append("RNA" + temps[i] + "\t");
                        }
                        bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
                        continue;
                    }
                    int num = 0;
//                    int total = (int) (temps.length * 0.4);
                    int total = temps.length - 9;
                    for (int i = 9; i < temps.length; i++) {
                        if (temps[i].split(";")[0].equals("./.")) {
                            num++;
                        }
                    }
                    if (num <= total * 0.4) {
                        bw.write(temp + "\n");
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
            }
        });
        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Getting subRNA genotype vcf takes " + (endTime - startTime) + "ms");
    }

//    public void getsubRNAgenotypeor() {
//        long startTime = System.currentTimeMillis();   //获取开始时间
//        System.out.println("This is getting subRNA genotype ************************************************************");
//        String RNADir = new File(outputDir, "RNA").getAbsolutePath();
//        String originalRNA = new File(outputDir, "VCF").getAbsolutePath();
//        HashSet<String> nameSet = new HashSet<>();
//        File[] fs = new File(originalRNA).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
//        for (int i = 0; i < fs.length; i++) {
//            nameSet.add(fs[i].getName());
//        }
//        BufferedWriter bw1 = IOUtils.getTextWriter(new File(RNADir, "header.txt").getAbsolutePath());
//        nameSet.parallelStream().forEach(f -> {
//            BufferedReader br = IOUtils.getTextReader(new File(originalRNA, f).getAbsolutePath());
//            BufferedWriter bw = IOUtils.getTextWriter(new File(RNADir, f).getAbsolutePath());
//            String temp = null;
//            String[] temps = null;
//            try {
//                while ((temp = br.readLine()) != null) {
//                    if (temp.startsWith("#")) {
//                        bw.write(temp + "\n");
//                        if (f.equals("chr001.vcf")) {
//                            bw1.write(temp + "\n");
//                        }
//                        continue;
//                    }
//                    temps = temp.split("\t");
//                    int num = 0;
////                    int total = (int) (temps.length * 0.4);
//                    int total = temps.length-9;
//                    for (int i = 9; i < temps.length; i++) {
//                        if (temps[i].split(";")[0].equals("./.")) {
//                            num++;
//                        }
//                    }
//                    if (num <= total * 0.4) {
//                        bw.write(temp + "\n");
//                    }
//                }
//                br.close();
//                bw.flush();
//                bw.close();
//                bw1.flush();
//                bw1.close();
//            } catch (Exception e) {
//            }
//        });
//        long endTime = System.currentTimeMillis(); //获取结束时间
//        System.out.println("******* Getting subRNA genotype vcf takes " + (endTime - startTime) + "ms");
//    }

    public void getMerge() {
        String infileDirRNA = new File(outputDir, "RNA").getAbsolutePath();
        String infileDirDNA = new File(outputDir, "Isec").getAbsolutePath();
        BufferedWriter bwRNA = IOUtils.getTextWriter(new File(infileDirRNA, "RNAall.vcf").getAbsolutePath());
        BufferedWriter bwDNA = IOUtils.getTextWriter(new File(infileDirDNA, "DNAall.vcf").getAbsolutePath());
        String temp = null;
        try {
            for (int i = 0; i < chrNumber; i++) {
                int chr = i + 1;
                BufferedReader br = IOUtils.getTextReader(new File(infileDirRNA, "RNA_chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf").getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (chr != 1 && temp.startsWith("#")) {
                        continue;
                    } else {
                        bwRNA.write(temp + "\n");
                    }
                }
                br.close();
            }
            bwRNA.flush();
            bwRNA.close();
            for (int i = 0; i < chrNumber; i++) {
                int chr = i + 1;
                BufferedReader br = IOUtils.getTextReader(new File(infileDirDNA, "DNA_chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf").getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (chr != 1 && temp.startsWith("#")) {
                        continue;
                    } else {
                        bwDNA.write(temp + "\n");
                    }
                }
                br.close();
            }
            bwDNA.flush();
            bwDNA.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
//        HashSet<String> RNAnameSet = new HashSet<>();
//        HashSet<String> DNAnameSet = new HashSet<>();

//        File[] fs = new File(infileDirRNA).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
//        for (int i = 0; i < fs.length; i++) {
//            RNAnameSet.add(fs[i].getName());
//        }
//
//        File[] fs1 = new File(infileDirDNA).listFiles();
//        fs1 = IOUtils.listFilesEndsWith(fs1, ".vcf");
//        for (int i = 0; i < fs1.length; i++) {
//            DNAnameSet.add(fs1[i].getName());
//        }

//        StringBuilder sb = new StringBuilder();
//        StringBuilder sb1 = new StringBuilder();
//
//        sb.append("vcf-concat ");
////        for (int i = 0; i < chrNumber; i++) {
////            int chr = i + 1;
////            sb.append(new File(infileDirRNA, "RNA_chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf").getAbsolutePath() + " ");
////        }
//        sb.append(infileDirRNA);
//        sb.append("/*.vcf ");
//        sb.append(" > " + new File(infileDirRNA, "RNAall.vcf").getAbsolutePath() + "\n");
//
//
//        sb.append("vcf-concat ");
////        for (int i = 0; i < chrNumber; i++) {
////            int chr = i + 1;
////            sb.append(new File(infileDirDNA, "DNA_chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf").getAbsolutePath() + " ");
////        }
//        sb.append(infileDirDNA);
//        sb.append("/*.vcf ");
//        sb.append(" > " + new File(infileDirDNA, "DNAall.vcf").getAbsolutePath() + "\n");
//
//        try {
//            String command = sb.toString();
//            System.out.println(command);
//            File dir = new File(new File(outputDir).getAbsolutePath());
//            String[] cmdarry = {"/bin/bash", "-c", command};
//            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//            p.waitFor();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
    }

    public void getIsec1() {
        long startTime = System.currentTimeMillis();   //获取开始时间
        System.out.println("This is getting intersection ***************************************************************");
        String infileDir1 = new File(outputDir, "RNA").getAbsolutePath();
        String infileDir2 = new File(genotypeDir).getAbsolutePath();
        String infileDir3 = new File(outputDir, "Isec").getAbsolutePath();
        HashSet<String> nameSet = new HashSet<>();
        File[] fs = new File(infileDir1).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName());
        }
        nameSet.parallelStream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                int chr = Integer.parseInt(String.valueOf(f).replace(".vcf", "").replace("RNA_chr", ""));
                String infile1 = new File(infileDir1, f).getAbsolutePath();
                String infile2 = new File(infileDir2, chr + genotypesuffix).getAbsolutePath();
                String outfile = new File(infileDir3, f.replace("RNA", "DNA")).getAbsolutePath();
                HashSet<String> positions = new HashSet<>();
                BufferedReader br = IOUtils.getTextReader(infile1);
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("#")) continue;
                    positions.add(temps[1]);
                }
                br.close();
                BufferedReader br1 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                while ((temp = br1.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    if (positions.contains(temps[1])) {
                        bw.write(temp + "\n");
                    }
                }
                br1.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

//    public void getIsec() {
//        long startTime = System.currentTimeMillis();   //获取开始时间
//        System.out.println("This is getting intersection ***************************************************************");
//        String infileDir1 = new File(outputDir, "RNA").getAbsolutePath();
//        String infileDir2 = new File(genotypeDir).getAbsolutePath();
//        String IsecDir = new File(outputDir, "Isec").getAbsolutePath();
//        HashSet<String> nameSet = new HashSet<>();
//        File[] fs = new File(infileDir1).listFiles();
//        fs = IOUtils.listFilesEndsWith(fs, ".vcf");
//        for (int i = 0; i < fs.length; i++) {
//            nameSet.add(fs[i].getName());
//        }
//        BufferedReader br = IOUtils.getTextGzipReader(new File(genotypeDir, "1" + genotypesuffix).getAbsolutePath());
//        BufferedWriter bw = IOUtils.getTextWriter(new File(infileDir1, "DNAheader.txt").getAbsolutePath());
//        String temp = null;
//        try {
//            while ((temp = br.readLine()) != null) {
//                if (temp.startsWith("#")) {
//                    bw.write(temp + "\n");
//                }
//                if (!temp.startsWith("#")) break;
//            }
//            br.close();
//            bw.flush();
//            bw.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        nameSet.parallelStream().forEach(f -> {
//            try {
//                int chr = Integer.parseInt(String.valueOf(f).replace(".vcf", "").replace("chr", ""));
//                String infile1 = new File(infileDir1, f).getAbsolutePath();
//                String infile2 = new File(infileDir2, chr + genotypesuffix).getAbsolutePath();
//                StringBuilder sb = new StringBuilder();
//                sb.append("bedtools intersect -a " + infile1 + " -b " + infile2 + " -wa > " + plate + "_chr" + chr + "_RNA_noheader.vcf\n");
//                sb.append("bedtools intersect -a " + infile2 + " -b " + infile1 + " -wa > " + plate + "_chr" + chr + "_DNA_noheader.vcf\n");
//
//                String command = sb.toString();
//                System.out.println(command);
//                File dir = new File(new File(outputDir, "Isec").getAbsolutePath());
//                String[] cmdarry = {"/bin/bash", "-c", command};
//                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                p.waitFor();
//
//                sb.setLength(0);
//                sb.append("cat ");
//                sb.append(new File(infileDir1, "header.txt").getAbsolutePath());
//                sb.append(" " + plate + "_chr" + chr + "_RNA_noheader.vcf > " + plate + "_chr" + chr + "_RNA.vcf\n");
//                sb.append("cat ");
//                sb.append(new File(infileDir1, "DNAheader.txt").getAbsolutePath());
//                sb.append(" " + plate + "_chr" + chr + "_DNA_noheader.vcf > " + plate + "_chr" + chr + "_DNA.vcf\n");
//                command = sb.toString();
//                System.out.println(command);
//                String[] cmdarry1 = {"/bin/bash", "-c", command};
//                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
//                p1.waitFor();
//
//                sb.setLength(0);
//                sb.append("rm " + plate + "_chr" + chr + "_RNA_noheader.vcf\n");
//                sb.append("rm " + plate + "_chr" + chr + "_DNA_noheader.vcf\n");
//                command = sb.toString();
//                System.out.println(command);
//                String[] cmdarry2 = {"/bin/bash", "-c", command};
//                Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir);
//                p2.waitFor();
//
//            } catch (Exception e) {
//                e.printStackTrace();
//            }
//        });
//        long endTime = System.currentTimeMillis(); //获取结束时间
//        System.out.println("******* Getting intersect takes " + (endTime - startTime) + "ms");
//    }

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
            sb1.append("cat " + new File(outputDir, "RNA/RNAall.vcf").getAbsolutePath() + " | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + new File(outputDir, "RNA/RNAall.sorted.vcf").getAbsolutePath());
            String command = sb1.toString();
            System.out.println(command);
            File dir = new File(infileDir).getAbsoluteFile();
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();

            StringBuilder sb2 = new StringBuilder();
            sb2.append("cat " + new File(outputDir, "Isec/DNAall.vcf").getAbsolutePath() + " | awk '$1 ~ /^#/ {print $0;next} {print $0 | \"sort -k1,1 -k2,2n\"}' > " + new File(outputDir, "RNA/RNAall.sorted.vcf").getAbsolutePath());
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
            for (int i = 0; i < chrNumber; i++) {
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
            for (int i = 0; i < chrNumber; i++) {
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
        String infileDir = new File(outputDir).getAbsolutePath();
        String infileS1 = new File(infileDir, "RNA/RNAall.vcf").getAbsolutePath();
        String infileS2 = new File(infileDir, "Isec/DNAall.vcf").getAbsolutePath();
//        String infileS1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/DNA.all.sort.vcf";
//        String infileS2 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/RNA.all.sort.vcf";
        String ibsOutfileS = new File(infileDir, "Isec/check.txt").getAbsolutePath();
//        String ibsOutfileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/input/check.txt";
        GenotypeGrid g1 = new GenotypeGrid(infileS1, GenoIOFormat.VCF);
        GenotypeGrid g2 = new GenotypeGrid(infileS2, GenoIOFormat.VCF);
        GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
        SumTaxaDivergence std = new SumTaxaDivergence(g);
        std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
        g.getIBSDistanceMatrix();
    }

    public void getDensityIBS() {
        System.out.println("This is writing file of RNA and DNA IBS plot ***********************************************");
        String infileDir = new File(outputDir, "Isec").getAbsolutePath();
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
            LinkedHashSet<String> samples = new LinkedHashSet<>();
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
            RowTable<String> t = new RowTable<>(infile);
            for (int i = 0; i < samplelist.length; i++) {
                String RNA = samplelist[i];
                String DNA = null;
                if (RNA.contains("JM22")) {
                    DNA = "E025";
                } else if (RNA.contains("CS")) {
                    DNA = "E360";
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

    public void mkDirHapscan() {
        for (int i = 0; i < subDirS.length; i++) {
            File f = new File(outputDirS, subDirS[i]);
            f.mkdir();
        }
    }

    public void parseParametersHapscan(String infileS) {
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

    public void Hapscanner() {
//        File[] fs = new File(this.parameterDir).listFiles();
//        fs = IOUtils.listFilesStartsWith(fs, plate + "_");
//        Arrays.sort(fs);
//
//        for (int i = 0; i < fs.length; i++) {
//            this.parseParametersHapscan(fs[i].getAbsolutePath());
//            this.mkDirHapscan();
//            this.scanIndiVCFByThreadPool();
//            this.mkFinalVCF();
//        }
        for (int i = 0; i < chrNumber; i++) {
            int chr = i + 1;
            File fs = new File(new File(parameterDir, plate + "_parameter_chr" + chr + ".txt").getAbsolutePath());
            this.parseParametersHapscan(fs.getAbsolutePath());
            this.mkDirHapscan();
            this.scanIndiVCFByThreadPool();
            this.mkFinalVCF();
        }
    }

    public void poswithAllele() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing pos file ****************************************************");
        HashSet<String> nameSet = new HashSet<>();
        File[] fs = new File(genotypeDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "360.vcf.gz");
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

    public void taxaRefBAM() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing taxaRefBam files *************************************************************************");
        String[] dirs = new File(this.BamDir).list();
        ArrayList<String> files = new ArrayList<>();
        if (plate.equals("S4leaf")) {
            for (int i = 0; i < dirs.length; i++) {
                String name = dirs[i].toString();
                String name1 = name.substring(name.length() - 2, name.length());
                if (Integer.parseInt(name1) > 70) {
                    if (!dirs[i].endsWith("DS_Store") && !dirs[i].startsWith("countTable")) {
                        File[] fs1 = new File(BamDir, dirs[i] + "/sams").listFiles();
                        fs1 = IOUtils.listFilesEndsWith(fs1, "_Aligned.out.sorted.bam");
                        for (int j = 0; j < fs1.length; j++) {
                            files.add(fs1[j].getAbsolutePath());
                            System.out.println(fs1[j]);
                        }
                    }
                }
            }
        } else {
            for (int i = 0; i < dirs.length; i++) {
                if (!dirs[i].endsWith("DS_Store") && !dirs[i].startsWith("countTable")) {
                    File[] fs1 = new File(BamDir, dirs[i] + "/sams").listFiles();
                    fs1 = IOUtils.listFilesEndsWith(fs1, "_Aligned.out.sorted.bam");
                    for (int j = 0; j < fs1.length; j++) {
                        files.add(fs1[j].getAbsolutePath());
                        System.out.println(fs1[j]);
                    }
                }
            }
        }

        HashSet<String> nameSet = new HashSet();
        for (
                int i = 0; i < files.size(); i++) {
            nameSet.add(files.get(i));
//            System.out.println(files.get(i));
        }
        try {
            String[] namelist = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < chrNumber; i++) {
                int chr = i + 1;
                BufferedWriter bw = IOUtils.getTextWriter(new File(taxaRefBAMDir, plate + "_taxaRefBAM_chr" + chr + ".txt").getAbsolutePath());
                bw.write("Taxa\tReference\tBamPath\n");
                for (int j = 0; j < namelist.length; j++) {
                    StringBuilder sb = new StringBuilder();
                    int length = namelist[j].split("/").length;
                    sb.append(namelist[j].split("/")[length - 1].replace("_Aligned.out.sorted.bam", "").replace("B18-", "")).append("\t");
                    sb.append(new File(referenceDir, "chr" + chr + ".fa").getAbsolutePath() + "\t");
                    sb.append(new File(namelist[j]).getAbsolutePath());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (
                Exception e) {
            e.printStackTrace();
        }

        long endTime = System.currentTimeMillis(); //获取结束时间
        System.out.println("******* Writing taxaRefBam files takes " + (endTime - startTime) + "ms");

    }

    public void parameter() {
        long startTime = System.currentTimeMillis();
        System.out.println("This is writing parameters files ***********************************************************");
        try {
            for (int i = 0; i < chrNumber; i++) {
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
                bw.write(new File(taxaRefBAMDir, plate + "_taxaRefBAM_chr" + chr + ".txt").getAbsolutePath() + "\n" +
                        "\n");
                bw.write("#Parameter 2: The posAllele file (with header), the format is Chr\\tPos\\tRef\\tAlt (from VCF format). The positions come from genetic variation library. \n" +
                        "#A maximum of 2 alternative alleles are supported, which is seperated by \",\", e.g. A,C.\n" +
                        "#Deletion and insertion are supported, denoted as \"D\" and \"I\".\n");
                bw.write(new File(posAlleleDir, "posAllele_chr" + chr + ".txt").getAbsolutePath() + "\n" +
                        "\n");
                bw.write("#Parameter 3: The pos files (without header), the format is Chr\\tPos. The positions come from haplotype library, which is used in mpileup.\n");
                bw.write(new File(posDir, "pos_chr" + chr + ".txt").getAbsolutePath() + "\n" +
                        "\n");
                bw.write("#Parameter 4: The chromosome which will be scanned.\n" +
                        chr + "\n" +
                        "\n" +
                        "#Parameter 5: Combined error rate of sequencing and misalignment. Heterozygous read mapping are more likely to be genotyped as homozygote when the combined error rate is high.\n" +
                        "0.05\n" +
                        "\n" +
                        "#Parameter 6: The path of samtools\n" +
                        samtoolsPath + "\n" +
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

    public void mkdir(String... args) {
        chrNumber = Integer.parseInt(args[1]);
    }


    public void parseparameters(String infileS) {
        Dyad<List<String>, List<String>> d = AppUtils.getParameterList(infileS);
        List<String> pLineList = d.getFirstElement();
        genotypeDir = pLineList.get(0);
        BamDir = pLineList.get(1);
        plate = pLineList.get(2);
        outputDir = pLineList.get(3) + "/" + plate;
        System.out.println(outputDir);

        parameterDir = pLineList.get(4);
        taxaRefBAMDir = pLineList.get(5);
        posDir = pLineList.get(6);
        posAlleleDir = pLineList.get(7);
        referenceDir = pLineList.get(8);

        samtoolsPath = pLineList.get(9);
        threads = pLineList.get(10);
//        chrNumber = Integer.parseInt(pLineList.get(11));


        File posdir = new File(new File(posDir).getAbsolutePath());
        File posAlleledir = new File(new File(posAlleleDir).getAbsolutePath());
        if (!posdir.exists() || !posAlleledir.exists()) {
            posdir.mkdir();
            posAlleledir.mkdir();
            this.poswithAllele();
        }

        System.out.printf("WTF");

        File output = new File(new File(outputDir).getAbsolutePath());
        output.mkdir();

        File RNAdir = new File(new File(outputDir, "RNA").getAbsolutePath());
        File Isecdir = new File(new File(outputDir, "Isec").getAbsolutePath());
        System.out.println(new File(outputDir, "RNA").getAbsolutePath());

        RNAdir.mkdir();
        Isecdir.mkdir();

        if (RNAdir.isDirectory()) {
            System.out.println("Yes");
            if (!RNAdir.exists()) {
                System.out.println("Ooops! Not made!");
                RNAdir.mkdir();
                if (RNAdir.exists()) {
                    System.out.println("Ohyeah! Made it");
                }
            } else {
                System.out.println("Already exist!");
            }
        }


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

    public static void main(String[] args) throws IOException, InterruptedException {
        new HapscannerParameters(args);
    }
}
