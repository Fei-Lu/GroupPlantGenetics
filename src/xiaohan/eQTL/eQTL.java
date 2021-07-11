package xiaohan.eQTL;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.koloboke.collect.impl.hash.Hash;
import com.koloboke.collect.impl.hash.MutableLHashDoubleSet;
import com.koloboke.collect.set.DoubleSet;
import daxing.common.MD5;
import daxing.load.ancestralSite.Standardization;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import xiaohan.eQTL.javapractice.javaPractice;
import xiaohan.rareallele.GeneFeature;
import xiaohan.rareallele.IOUtils;
import xiaohan.rareallele.VCFsplit;
import xiaohan.rareallele.pheno;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

public class eQTL {

    int score = 0;
    String kinship = "/data1/home/xiaohan/nam/KinshipMap.txt";
    Multimap<String, String> kinshipMap = ArrayListMultimap.create();
    HashSet<String> HBSet = new HashSet<>();

    public eQTL(String[] args) throws IOException, InterruptedException {
         /*
//        simulation
//        */
//        this.simulationData(args);

        /*
//        Meta-Tissue Analysis pipeline
//         */
//        this.metasoft(args);
//
         /*
//        Hapscanner parameter
//         */
//        this.hapscanner(args);

        /*
        effect size
         */
//        this.effectsize(args);

        /*
        triads
         */
//        this.getTraidsPattern();
//        this.patternIdentify();

        /*
        VCF
         */
//        this.vcf(args);

        /*
        pheno
         */
//        this.pheno(args);
//        this.covaraties(args);

        /*
        annotation
         */
//        this.annotation(args);
//        this.test();
//        this.getExpressionWithaFC(args[0], args[1]);

        /*
        meta-tissue analysis
         */
//        this.multiTissue(args);
//
        /*
        hetero
         */
//        this.getheterogeneity(args[0],args[1]);

        /*
        heterozygosity
         */
//        this.heterozygosity(args);
//        this.GBSsimulation(args);
//
//
//        this.intergenicPattern();
//        this.intergenicPatternEnrichment();
//        this.intergenicPatternTransposon();
//        this.getTransposonLength();
//        this.getTransposonClass();
//        this.command();
//        this.depth();
//        this.getsigcis();

        /*
        java with faster java practice
         */
//        this.javaPractice(args);
//        this.ERCCRoc();
//        this.getkinship();
//        this.getquality();
//        this.ThreadPool(args[0],args[1]);
//        this.subextractGeneRegion(args[0], args[1]);
//        this.hetandMaf(args[0], args[1]);
//        this.testforfunctional();
//        this.temp(args);
//        this.md5check();
//        this.getkinship();
//        this.getIBS();
//        this.getmedZ();
        this.testIBS(args);
//        this.testdensity(args);
//        this.filter(args);
//        this.getmedZ();
//        this.getsingleZ();
//        this.print();
    }

    public void filter(String[] args) {
        long startTime = System.currentTimeMillis();   //获取开始时间
        System.out.println("This is filtering samples ******************************************************************");
        String infile = new File(args[3]).getAbsolutePath();
        String infor = new File(args[2]).getAbsolutePath();
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


            BufferedWriter bw = IOUtils.getTextWriter(new File(args[5]).getAbsolutePath());
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

    public void testdensity(String[] args) {
        System.out.println("This is writing file of RNA and DNA IBS plot ***********************************************");
        String infile = new File(args[2]).getAbsolutePath();
        String outfile = new File(args[3]).getAbsolutePath();
        String outfile1 = new File(args[4]).getAbsolutePath();
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

    public void print() {
        int[] a = new int[]{42, 47, 57, 59, 69, 73, 74, 77, 79, 81};
        for (int i = 0; i < a.length; i++) {
            for (int j = 1; j < 3; j++) {
//                System.out.println("nohup zcat PM-XS01KF2020120443-08遗传所19个外来文库测序过滤任务单-加测/ANNO_XS01KF2020120443_PM-XS01KF2020120443-08_BHV337DSXY_2021-03-02_22-22-47/Cleandata/20201228SIPASR" + a[i] + "/ 20201228SIPASR" + a[i] + "_R"+j+".fq.gz /data2/junxu/Data/ANNO_XS01KF2020120443_PM-XS01KF2020120443-11_AHVYTFDSXY_2021-04-11_21-30-57/Cleandata/20201228SIPASR"+a[i]+"/20201228SIPASR"+a[i]+"_R"+j+".fq.gz | gzip > /data2/xiaohan/MergedData/20201228SIPASR"+a[i]+"/20201228SIPASR"+a[i]+"_R"+j+".fq.gz &");
                System.out.println(String.valueOf(a[i]));
            }
        }
    }

    public void testIBS(String[] args) {
        String infile1 = args[0];
        String infile2 = args[1];
        String outfile = args[2];
        String vcf1 = new File(infile1).getAbsolutePath();
        String vcf2 = new File(infile2).getAbsolutePath();
        GenotypeGrid g1 = null;
        GenotypeGrid g2 = null;
        if (infile1.endsWith("gz")) {
            g1 = new GenotypeGrid(vcf1, GenoIOFormat.VCF_GZ);
        } else g1 = new GenotypeGrid(vcf1, GenoIOFormat.VCF);
        if (infile2.endsWith("gz")) {
            g2 = new GenotypeGrid(vcf2, GenoIOFormat.VCF_GZ);
        } else g2 = new GenotypeGrid(vcf2, GenoIOFormat.VCF);
        String ibsOutfileS = new File(outfile).getAbsolutePath();
        GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
        SumTaxaDivergence std = new SumTaxaDivergence(g);
        std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
        g.getIBSDistanceMatrix();
    }

    public void getsingleZ() {
        String infile = "/data2/xiaohan/pheno/rareMutation/all_normalized_expression.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infile);
        String outfile = "/data2/xiaohan/pheno/rareMutation/outliers_medz_picked.txt";
        String temp = null;
        List<String> temps = null;
        HashSet<String> tissueSet = new HashSet<>();
        HashSet<String> geneSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        int countline1 = 0;
        int countline = 0;
        int columnnumber = 0;
        int index = -1;
        String max = null;
        String error = null;
        try {
            String[] header = null;
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 5000 == 0) {
                    System.out.println(countline);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) {
                    header = new String[temps.size() - 2];
                    for (int i = 2; i < temps.size(); i++) {
                        header[i - 2] = temps.get(i);
                    }
                    columnnumber = temps.size() - 2;
                    System.out.println(columnnumber);
                    continue;
                }
                tissueSet.add(temps.get(0));
            }
            br.close();
            String[] tissues = tissueSet.toArray(new String[0]);
            BufferedWriter[] bw = new BufferedWriter[tissues.length];
            BufferedWriter[] bw1 = new BufferedWriter[tissues.length];
            for (int i = 0; i < bw.length; i++) {
                bw[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/rareMutation/" + tissues[i] + "_outliers_singlez_picked.txt").getAbsolutePath());
                bw1[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/rareMutation/" + tissues[i] + "_outliers_counts.txt").getAbsolutePath());
            }
            while ((temp = br1.readLine()) != null) {
                countline1++;
                if (countline1 % 5 == 0) {
                    System.out.println(countline1);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) {
                    sb.setLength(0);
                    for (int i = 1; i < temps.size(); i++) {
                        sb.append(temps.get(i) + "\t");
                    }
                    for (int i = 0; i < bw1.length; i++) {
                        bw1[i].write(sb.toString());
                        bw1[i].newLine();
                    }
                    continue;
                }
                index = Arrays.binarySearch(tissues, temps.get(0));

            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void getmedZ() {
        String infile = "/data2/xiaohan/rareMutation/pheno/all_normalized_expression.txt";
//        String infile = "/Users/yxh/Documents/RareAllele/script/001pheno_processing/all_normalized_expression.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br1 = IOUtils.getTextReader(infile);
        String outfile = "/data2/xiaohan/rareMutation/pheno/outliers_medz_picked.txt";
        String outfilecount = "/data2/xiaohan/rareMutation/pheno/outliers_medz_counts.txt";
//        String outfile = "/Users/yxh/Documents/RareAllele/script/001pheno_processing/outliers_medz_picked.txt";
        String temp = null;
        List<String> temps = null;
        HashSet<String> tissueSet = new HashSet<>();
        HashSet<String> geneSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        int countline1 = 0;
        int countline = 0;
        int columnnumber = 0;
        String error = null;

        //picked variables
        String gene = null;
        String name = null;
        String tissue = null;
        String max = null;
        try {
            String[] header = null;
            String headers = null;
            System.out.println("Starting initiation ----------------------------------");
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 5000 == 0) {
                    System.out.println(countline);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) {
                    header = new String[temps.size() - 2];
                    for (int i = 2; i < temps.size(); i++) {
                        header[i - 2] = temps.get(i);
                    }
                    sb.setLength(0);
                    sb.append("GENE\t");
                    for (int i = 0; i < header.length; i++) {
                        sb.append(header[i] + "\t");
                    }
                    headers = sb.toString().replaceAll("\\s+$", "");
                    columnnumber = temps.size() - 2;
                    System.out.println(columnnumber);
                    continue;
                }
                tissueSet.add(temps.get(0));
                geneSet.add(temps.get(1));
            }
            br.close();
            String[][][] expression = new String[geneSet.size()][columnnumber][tissueSet.size()];
            String[] tissues = tissueSet.toArray(new String[0]);
            String[] genes = geneSet.toArray(new String[0]);
            HashMap<String, Integer> tissuesIndex = new HashMap<>();
            HashMap<String, Integer> genesIndex = new HashMap<>();
            for (int i = 0; i < tissues.length; i++) {
                tissuesIndex.put(tissues[i], i);
                for (int j = 0; j < genes.length; j++) {
                    genesIndex.put(genes[j], j);
                    for (int k = 0; k < columnnumber; k++) {
                        expression[j][k][i] = "NA";
                    }
                }
            }
            System.out.println("Starting building expression file ----------------------------------");
            while ((temp = br1.readLine()) != null) {
                countline1++;
                if (countline1 % 5 == 0) {
//                    System.out.println(countline1);
                }
                temps = PStringUtils.fastSplit(temp);
                if (temp.startsWith("Tissue")) continue;
                for (int i = 0; i < columnnumber; i++) {
                    expression[genesIndex.get(temps.get(1))][i][tissuesIndex.get(temps.get(0))] = temps.get(i + 2);
                }
            }
            System.out.println("Finished building expression file !");
            br1.close();
            String[][] expressionmedian = new String[geneSet.size()][columnnumber];
            String[][] expressionN = new String[geneSet.size()][columnnumber];
            String[] array = new String[tissues.length];
            int len = -1;
            double median = -1;
            System.out.println("Starting calculating n and median ----------------------------------");
            for (int i = 0; i < genes.length; i++) {
                for (int j = 0; j < columnnumber; j++) {
                    array = expression[i][j];
                    expressionN[i][j] = medn(array);
//                    System.out.println(medn(array));
                    expressionmedian[i][j] = String.valueOf(medmedian(array));
                }
            }
            // writing medz picked file
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            System.out.println("Starting picking outliers ----------------------------------");
            for (int i = 0; i < genes.length; i++) {
                sb.setLength(0);
                sb.append(genes[i] + "\t");
                String[] outliers = pickoutlier(expressionmedian[i]);
//                System.out.println(genes[i] + "_" + expressionmedian[i][0] +"_"+ outliers[0] +"_"+ outliers[1]);
                sb.append(header[Integer.parseInt(outliers[1])] + "\t");
                sb.append(expressionN[i][Integer.parseInt(outliers[1])] + "\t");
                sb.append(outliers[0]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            // writing medz counts file
            BufferedWriter bwcount = IOUtils.getTextWriter(outfilecount);
            bwcount.write(headers+"\n");
            for (int i = 0; i < genes.length; i++) {
                sb.setLength(0);
                sb.append(genes[i]+"\t");
                for (int j = 0; j < columnnumber; j++) {
                    sb.append(expressionN[i][j]+"\t");
                }
                bwcount.write(sb.toString().replaceAll("\\s+$","")+"\n");
            }
            bwcount.flush();
            bwcount.close();
            // writing single counts file and picked file
            String[] array1 = new String[columnnumber];
            String[] outliers = new String[2];
            BufferedWriter[] bw1 = new BufferedWriter[tissues.length];
            BufferedWriter[] bw2 = new BufferedWriter[tissues.length];
            for (int i = 0; i < tissues.length; i++) {
                bw1[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/rareMutation/" + tissues[i] + "_outliers_counts.txt").getAbsolutePath());
                bw2[i] = IOUtils.getTextWriter(new File("/data2/xiaohan/pheno/rareMutation/" + tissues[i] + "_outliers_singlez_picked.txt").getAbsolutePath());
                bw1[i].write(headers + "\n");
                bw2[i].write("GENE\tIND\tDFS\tZ\n");
            }
            for (int i = 0; i < genes.length; i++) {
                for (int j = 0; j < tissues.length; j++) {
                    gene = genes[i];
                    sb.setLength(0);
                    sb.append(gene+"\t");
                    for (int k = 0; k < columnnumber; k++) {
                        array1[k] = expression[i][k][j];
                        if (array1[k].equals("NA")) {
                            sb.append("0\t");
                        } else {
                            sb.append("1\t");
                        }
                    }
                    outliers = pickoutlier(array1);
                    if (outliers[0].equals("NA") && outliers[1].equals("NA")) continue;
                    name = header[Integer.parseInt(outliers[1])];
                    tissue = tissues[j];
                    max = expression[i][Integer.parseInt(outliers[1])][j];
                    bw1[j].write(sb.toString().replaceAll("\\s+$", "") + "\n");
                    bw2[j].write(gene + "\t" + name + "\t" + tissue + "\t" + max + "\n");
                }
            }
            for (int i = 0; i < tissues.length; i++) {
                bw1[i].flush();
                bw1[i].close();
                bw2[i].flush();
                bw2[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
//            System.out.println(countline1);
//            System.out.println(countline);
//            System.out.println(columnnumber);
            System.out.println(gene);
        }
    }

    private String medn(String[] array) {
        int len = 0;
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                len++;
            }
        }
        String length = String.valueOf(len);
        return length;
    }

    private String medmedian(String[] array) {
        HashSet<Double> DataSet = new HashSet<>();
        String median = null;
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                DataSet.add(Double.parseDouble(array[i]));
            }
        }
        if (DataSet.isEmpty()) {
            median = "NA";
        } else {
            Object[] data = DataSet.toArray();
            if (data.length % 2 == 1) {
                median = String.valueOf((double) data[data.length / 2]);
            } else {
                median = String.valueOf(((double) data[data.length / 2] + (double) data[data.length / 2 - 1]) / 2);
            }
        }
        return median;
    }

    private String[] pickoutlier(String[] array) {
        HashMap<String, Integer> IndexValues = new HashMap<>();
        for (int i = 0; i < array.length; i++) {
            IndexValues.put(array[i], i);
        }
        String max = null;
        String rank = null;
        HashSet<Double> DataSet = new HashSet<>();
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                DataSet.add(Double.parseDouble(array[i]));
            }
        }
        if (DataSet.isEmpty()) {
            max = "NA";
            rank = "NA";
        } else {
            Object[] data = DataSet.toArray();
            double[] abs = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                abs[i] = Math.abs((double) data[i]);
            }
            int indexmax = (int) getMaxIndex(abs)[1];
            max = String.valueOf((double) data[indexmax]);
            rank = String.valueOf(IndexValues.get(max));
        }
        String[] results = new String[]{max, rank};
        return results;
    }

    public static double[] getMaxIndex(double[] arr) {
        if (arr == null || arr.length == 0) {
            return null;//如果数组为空 或者是长度为0 就返回null
        }
        int maxIndex = 0;//假设第一个元素为最大值 那么下标设为0
        double[] arrnew = new double[2];//设置一个长度为2的数组用作记录，第一个元素存储最大值，第二个元素存储下标
        for (int i = 0; i < arr.length - 1; i++) {
            if (arr[maxIndex] < arr[i + 1]) {
                maxIndex = i + 1;
            }
        }
        arrnew[0] = arr[maxIndex];
        arrnew[1] = maxIndex;
        return arrnew;
    }

    private String expTissuNumber(String[] array) {
        String max = null;
        HashSet<Double> DataSet = new HashSet<>();
        for (int i = 0; i < array.length; i++) {
            if (!array[i].equals("NA")) {
                DataSet.add(Double.parseDouble(array[i]));
            }
        }
        if (DataSet.isEmpty()) {
            max = "0";
        } else {
            max = String.valueOf(DataSet.size());
        }
        return max;
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

    private void getIBS() {
        String inputDir = new File("/data1/home/xiaohan/hapscanner/output/test4-1/VCF").getAbsolutePath();
        HashSet<Integer> chrSet = new HashSet<>();
        StringBuilder sb = new StringBuilder();
        try {
            for (int m = 0; m < 21; m++) {
                int chr1 = m * 2 + 1;
                int chr2 = m * 2 + 2;
                String chrABDName = chrUtils.getNamechrABD(m + 1);
//                sb.setLength(0);
//                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr1) + ".vcf").getAbsolutePath() + "\n");
//                sb.append("bgzip " + new File(inputDir, "chr" + PStringUtils.getNDigitNumber(3, chr2) + ".vcf").getAbsolutePath() + "\n");
//                System.out.println(sb.toString());
//                File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
//                String[] cmdarry = {"/bin/bash", "-c", sb.toString()};
//                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
//                p.waitFor();

                String vcf1 = new File(inputDir, chrABDName + ".vcf.gz").getAbsolutePath();
                String vcf2 = new File("/data1/home/xiaohan/hapscanner/output/Parent/VCF", chrABDName + ".vcf.gz").getAbsolutePath();

                GenotypeGrid g1 = new GenotypeGrid(vcf1, GenoIOFormat.VCF_GZ);
                GenotypeGrid g2 = new GenotypeGrid(vcf2, GenoIOFormat.VCF_GZ);

                String ibsOutfileS = new File("/data1/home/xiaohan/nam/test4-1/IBS/check_" + chrABDName + ".txt").getAbsolutePath();

                GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
                SumTaxaDivergence std = new SumTaxaDivergence(g);
                std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
                g.getIBSDistanceMatrix();
                xiaohan.eQTL.RowTable<String> rt = new xiaohan.eQTL.RowTable<>(ibsOutfileS);
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

                BufferedWriter bw = IOUtils.getTextWriter(new File("/data1/home/xiaohan/nam/test4-1/sum/sum" + chrABDName + ".txt").getAbsolutePath());
                for (int i = 0; i < HBset.size(); i++) {
                    String HB = HBset.get(i);
                    String HBnickname = HB.substring(0, 5);
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

    public void md5check() {
        String inputFile = "/Volumes/Lulab4T_87/PM-XS01KF2020120443-13遗传所6个外来文库测序过滤任务单/ANNO_XS01KF2020120443_PM-XS01KF2020120443-13_BHW2JWDSXY_2021-03-30_14-05-46/md5.txt";
        String inputDir = "/Volumes/Lulab4T_87/PM-XS01KF2020120443-13遗传所6个外来文库测序过滤任务单/ANNO_XS01KF2020120443_PM-XS01KF2020120443-13_BHW2JWDSXY_2021-03-30_14-05-46/";
        BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(inputFile);
        try {
            String temp = null;
            String temp1 = null;
            String file = null;
            while ((temp = br.readLine()) != null) {
                temp1 = temp.split(" ")[0];
                file = new File(inputDir, temp.split(" ")[2]).getAbsolutePath();
//                System.out.println(file);
                if (MD5.checkMD5(file, temp1)) {
                    System.out.println("True");
                } else System.out.printf("False");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void temp(String[] args) {
        String barcodeFileS = "/data1/home/xiaohan/nam/GBStest0630/barcode.txt";
        String libraryFastqMapFileS = "/data1/home/xiaohan/nam/GBStest0630/" + args[0] + "libraryFastqMap.txt";
        String cutter1 = "GCGGCCGC";
        String cutter2 = "CCGG";
        String outputDir = args[1];
        LibraryOfGRTcopy li = new LibraryOfGRTcopy(barcodeFileS, libraryFastqMapFileS, cutter1, cutter2);
        li.splitBarcode(new File(outputDir).getAbsolutePath());

    }

    public void hetandMaf(String infileDir, String outfileDir) {
        HashSet<Integer> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(chr);
        }
        nameSet.parallelStream().forEach(f -> {
            try {
                String infile = new File(infileDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
                String outfile = new File(outfileDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".txt").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(infile);
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                String temp = null;
                String[] temps = null;
                bw.write("chr\tpos\thet\tmaf\tdepth\n");
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        continue;
                    }
                    temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                    String info = temps[7];
                    String[] infos = info.split(";");
                    String het = infos[5].split("=")[1];
                    String maf = infos[6].split("=")[1];
                    String depth = infos[2].split("=")[1];
                    bw.write(temps[0] + "\t" + temps[1] + "\t" + het + "\t" + maf + "\t" + depth + "\n");
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }


    public void subextractGeneRegion(String infileDir, String outfileDir) {
        String annotationfile = "/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3";
        HashSet<Integer> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(chr);
        }
        nameSet.parallelStream().forEach(f -> {
            String infile = new File(infileDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
            String outfile = new File(outfileDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            String temp = null;
            String[] temps = null;
            try {
                GeneFeature gf = new GeneFeature(annotationfile);
                HashSet<String> geneSet = new HashSet<>();
                for (int i = 0; i < gf.getGeneNumber(); i++) {
                    if (gf.getGeneChromosome(i) == f) {
                        geneSet.add(gf.getGeneName(i));
                    }
                }
                String[] genelist = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genelist);
                System.out.println("Finished gathering genelist");

                int[][] geneRange = new int[genelist.length][3];
                for (int i = 0; i < genelist.length; i++) {
                    geneRange[i][0] = gf.getGeneChromosome(gf.getGeneIndex(genelist[i]));
                    geneRange[i][1] = gf.getGeneStart(gf.getGeneIndex(genelist[i]));
                    geneRange[i][2] = gf.getGeneEnd(gf.getGeneIndex(genelist[i]));
                }

                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    int snppos = Integer.parseInt(PStringUtils.fastSplit(temp).toArray(new String[0])[1]);
                    int[] index = SNPmappingInGene.binarySearch(geneRange, snppos);
                    if (index[0] == -1) {
                        continue;
                    } else {
                        bw.write(temp + "\n");
                    }
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getExpressionWithaFC(String plate, String outfile) {
        //reading file with effect size
        String[] plates = plate.split(",");
        String[] infor = new String[plates.length];
        String[] infile = new String[plates.length];
        for (int i = 0; i < plates.length; i++) {
            infor[i] = "/data2/xiaohan/pheno/" + plates[i] + "/" + plates[i] + "_expressionMedian.txt";
            infile[i] = "/data2/xiaohan/tensorResult/effect/" + plates[i] + "/all.cis.top.DE.log2.txt";
        }

        BufferedReader[] br = new BufferedReader[plates.length];
        for (int i = 0; i < infile.length; i++) {
            System.out.println(infile[i]);
            if (infile[i].endsWith("gz")) {
                br[i] = IOUtils.getTextGzipReader(infile[i]);
            } else {
                br[i] = IOUtils.getTextReader(infile[i]);
            }
        }

        String temp = null;
        String[] temps = null;
        try {
            //Get geneset united
            HashSet<String> geneSet = new HashSet<>();
            HashMap<String, Integer> geneindexMap = new HashMap<>();
            for (int i = 0; i < infile.length; i++) {
                int pindex = 0;
                int aFCindex = 0;
                temp = br[i].readLine();
                temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                for (int k = 0; k < temps.length; k++) {
                    if (temps[k].startsWith("pheno")) {
                        pindex = k;
                    }
                    if (temps[k].equals("log2_aFC")) {
                        aFCindex = k;
                    }
                }
                while ((temp = br[i].readLine()) != null) {
                    if (temp.startsWith("p")) continue;
                    temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                    geneSet.add(temps[pindex]);
                }
                br[i].close();
            }

            String[] genelist = geneSet.toArray(new String[0]);
            Arrays.sort(genelist);

            //Initiating a new countTable
            String[][] countTable = new String[genelist.length][infor.length * 2 + 1];
            for (int i = 0; i < genelist.length; i++) {
                for (int j = 0; j < infor.length * 2 + 1; j++) {
                    countTable[i][j] = String.valueOf(0);
                }
            }

            for (int i = 0; i < genelist.length; i++) {
                geneindexMap.put(genelist[i], i);
                countTable[i][0] = genelist[i];
            }

            //Reading gene expression median files
            for (int i = 0; i < infor.length; i++) {
                if (infor[i].endsWith("gz")) {
                    br[i] = IOUtils.getTextGzipReader(infor[i]);
                } else {
                    br[i] = IOUtils.getTextReader(infor[i]);
                }
            }

            for (int i = 0; i < infor.length; i++) {
                while ((temp = br[i].readLine()) != null) {
                    if (!temps[0].startsWith("T")) continue;
                    temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                    //locating gene index with infor index to gene expression
                    if (geneSet.contains(temps[0])) {
                        countTable[geneindexMap.get(temps[0])][i + 1] = temps[1];
                    }
                }
                br[i].close();
            }

            //Reading top variants effect size files
            for (int i = 0; i < infile.length; i++) {
                if (infile[i].endsWith("gz")) {
                    br[i] = IOUtils.getTextGzipReader(infile[i]);
                } else {
                    br[i] = IOUtils.getTextReader(infile[i]);
                }
            }

            for (int i = 0; i < infile.length; i++) {
                int pindex = 0;
                int aFCindex = 0;
                temp = br[i].readLine();
                temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                for (int k = 0; k < temps.length; k++) {
                    if (temps[k].startsWith("pheno")) {
                        pindex = k;
                    }
                    if (temps[k].equals("log2_aFC")) {
                        aFCindex = k;
                    }
                }
                while ((temp = br[i].readLine()) != null) {
                    if (temp.startsWith("p")) continue;
                    temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                    countTable[geneindexMap.get(temps[pindex])][i + 1 + infor.length] = temps[aFCindex];
                }
                br[i].close();
            }

            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < genelist.length; i++) {
                for (int j = 0; j < infor.length * 2 + 1; j++) {
                    sb.append(countTable[i][j] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
                sb.setLength(0);
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void getheterogeneity(String infile, String outfile) {
        BufferedReader br = null;
        BufferedReader br1 = null;
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(infile);
            br1 = IOUtils.getTextGzipReader(infile);
        } else {
            br = IOUtils.getTextReader(infile);
            br1 = IOUtils.getTextReader(infile);
        }
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        Multimap<String, String> eQTLeGeneMap = ArrayListMultimap.create();
        Multimap<String, String> eGeneeQTLMap = ArrayListMultimap.create();
        try {
            int index = 0;
            int pindex = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("p") || temp.startsWith("In")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals("variant_id")) {
                            index = i;
//                            System.out.println(index);
                        }
                        if (temps[i].equals("phenotype_id")) {
                            pindex = i;
                        }
                    }
                    bw.write("phenotype_id\tvariant_id\teQTLs/pereGene\teGenes/pereQTL");
                    bw.newLine();
                    continue;
                }
                temps = temp.split("\t");
                String gene = temps[pindex];
                String variant = temps[index];
                eGeneeQTLMap.put(gene, variant);
                eQTLeGeneMap.put(variant, gene);
//                System.out.println(gene);
//                System.out.println(variant);
            }
            br.close();
            StringBuilder sb = new StringBuilder();
            while ((temp = br1.readLine()) != null) {
                if (temp.startsWith("p") || temp.startsWith("In")) continue;
                temps = temp.split("\t");
                String[] eQTLs = eGeneeQTLMap.get(temps[pindex]).toArray(new String[0]);
                String[] eGenes = eQTLeGeneMap.get(temps[index]).toArray(new String[0]);
                sb.setLength(0);
                sb.append(temps[pindex] + "\t" + temps[index] + "\t");
                sb.append(eQTLs.length + "\t");
                sb.append(eGenes.length);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void ThreadPool(String infile, String thread) {
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        ExecutorService pool = Executors.newFixedThreadPool(Integer.parseInt(thread));
        File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                String command = temp;
                covaraties.Command com = new covaraties.Command(command, dir);
                Future<covaraties.Command> chrom = pool.submit(com);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void covaraties(String[] args) {
        new covaraties(args);
    }

    public void ERCCRoc() {
        String prefix = "Truseq";
        int sample = 3;
        BufferedWriter bw = IOUtils.getTextWriter("/Users/yxh/Documents/eQTL/SiPAS/ERCC/20210318/5M/" + prefix + "/" + prefix + "PR.txt");
        RowTable<String> rt = new RowTable<>("/Users/yxh/Documents/eQTL/SiPAS/ERCC/20210318/5M/" + prefix + "/" + prefix + "Roc.txt");
        try {
            bw.write("Gene\tMix1_C\tMix1_precision\tMix1_recall\tMix2_C\tMix2_precision\tMix2_recall\n");
            for (int i = 0; i < rt.getRowNumber(); i++) {
                double[] TP = new double[sample * 2];
                double[] FP = new double[sample * 2];
                double[] FN = new double[sample * 2];
                for (int j = 1; j <= 2; j++) {
                    for (int k = 1; k <= sample; k++) {
                        String index = j + "_" + k;
                        String observed = rt.getCell(i, rt.getColumnIndex("log2mix" + index));
                        String expected = rt.getCell(i, rt.getColumnIndex("predictmix" + index));
                        double[] TPFPFN = xiaohan.eQTL.simulationData.getPrecisionandRecall(expected, observed);
                        int ix = (j - 1) * sample + k - 1;
                        TP[ix] = TPFPFN[0];
                        FP[ix] = TPFPFN[1];
                        FN[ix] = TPFPFN[2];
                    }
                }
                StringBuilder sb = new StringBuilder();
                sb.append(rt.getCell(i, rt.getColumnIndex("Gene")) + "\t");
                sb.append(rt.getCell(i, rt.getColumnIndex("Mix1_C")) + "\t");
                DecimalFormat decfor = new DecimalFormat("0.00000000");
                double precision1 = 0.00000000;
                double recall1 = 0.00000000;
                double precision2 = 0.00000000;
                double recall2 = 0.00000000;
                double TP1 = 0.00000000;
                double TP2 = 0.00000000;
                double FP1 = 0.00000000;
                double FP2 = 0.00000000;
                double FN1 = 0.00000000;
                double FN2 = 0.00000000;
                for (int j = 0; j < sample; j++) {
                    TP1 += TP[j];
                    FP1 += FP[j];
                    FN1 += FN[j];
                }
                for (int j = sample; j < sample * 2; j++) {
                    TP2 += TP[j];
                    FP2 += FP[j];
                    FN2 += FN[j];
                }
                precision1 = (double) TP1 / (TP1 + FP1);
                recall1 = (double) TP1 / (TP1 + FN1);
                precision2 = (double) TP2 / (TP2 + FP2);
                recall2 = (double) TP2 / (TP2 + FN2);
                sb.append(precision1 + "\t" + recall1 + "\t");
                sb.append(rt.getCell(i, rt.getColumnIndex("Mix2_C")) + "\t");
                sb.append(precision2 + "\t" + recall2);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void test(String... args) {
        VCFsplit.main(args);
    }

    public void javaPractice(String[] args) {
        new javaPractice(args);
    }

    public void pheno(String[] args) {
        new pheno(args);
    }

    public void annotation(String[] args) {
        new annotation(args);
    }

    public void vcf(String[] args) {
        new vcf(args);
    }

    public void effectsize(String[] args) {
        new effectsize(args);
    }

    public void multiTissue(String[] args) {
        new multiTissue(args);
    }

    public void depth() {
        String infileDir = "/data1/home/xiaohan/hapscanner/testforerror/VCF";
        String outfile = "/data1/home/xiaohan/nam/subsample/";
        HashSet<Integer> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(chr);
        }
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(infileDir, "chr" + PStringUtils.getNDigitNumber(3, f) + ".vcf").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outfile, "chr" + PStringUtils.getNDigitNumber(3, f) + ".depth").getAbsolutePath());
                StringBuilder sb = new StringBuilder();
                String temp = null;
                String[] temps = null;
                String tems = null;
                int depth = 0;
                while ((temp = br.readLine()) != null) {
                    sb.setLength(0);
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        sb.append(temps[0] + "\t" + temps[1] + "\t");
                        for (int i = 9; i < temps.length; i++) {
                            sb.append(temps[i] + "\t");
                        }
                        bw.write(sb.toString().replaceAll("\\s+$", ""));
                        bw.newLine();
                        continue;
                    }
                    temps = temp.split("\t");
                    sb.append(temps[0] + "\t" + temps[1] + "\t");
                    for (int i = 9; i < temps.length; i++) {
                        tems = temps[i].split(":")[1];
                        if (tems.equals(".")) {
                            depth = 0;
                        } else {
                            depth = Integer.parseInt(tems.split(",")[0]) + Integer.parseInt(tems.split(",")[1]);
                        }
                        sb.append(depth + "\t");
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void heterozygosity(String[] args) throws IOException, InterruptedException {
        new heterozygosity(args);
    }

    public void GBSsimulation(String[] args) throws IOException, InterruptedException {
        new GBSsimulation(args);
    }

    public void command() {
        String infileDir = "/data1/home/xiaohan/hapscanner/testforerror/VCF";
        String infor = "/data1/home/xiaohan/nam/subsample/basename.txt";
        String outfileDir = "/data1/home/xiaohan/nam/subsample/";
        HashSet<String> chrSet = new HashSet<>();
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            chrSet.add(String.valueOf(chr));
        }
        BufferedWriter bw = IOUtils.getTextWriter(new File(outfileDir, "sum.txt").getAbsolutePath());
        chrSet.stream().forEach(f -> {
            try {

                BufferedReader brinfor = IOUtils.getTextReader(infor);
                String temp = null;
                String[] temps = null;
                HashSet<String> nameSet = new HashSet<>();
                while ((temp = brinfor.readLine()) != null) {
                    temps = temp.split("-");
                    String name = temp.replace("-" + temps[temps.length - 1], "");
                    if (!nameSet.contains(name)) {
                        nameSet.add(name);
                    }
                }
                brinfor.close();
                String[] namelist = nameSet.toArray(new String[0]);
                Arrays.sort(namelist);
                int[][] count = new int[namelist.length][4];
                for (int i = 0; i < namelist.length; i++) {
                    for (int j = 0; j < 4; j++) {
                        count[i][j] = 0;
                    }
                }

                BufferedReader br = IOUtils.getTextReader(new File(infileDir, "chr" + PStringUtils.getNDigitNumber(3, Integer.parseInt(f)) + ".vcf").getAbsolutePath());
                HashMap<String, Integer> nameMap = new HashMap<>();
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 5000 == 0) {
                        System.out.println(f + ": " + countline);
                    }
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        for (int i = 9; i < temps.length; i++) {
                            nameMap.put(temps[i], i);
                        }
                        continue;
                    }
                    temps = temp.split("\t");
                    for (int i = 0; i < namelist.length; i++) {
//                        System.out.println(namelist[i]);
                        if (!temps[nameMap.get(namelist[i] + "-low")].split(":")[0].equals("./.")) {
                            count[i][0] += 1;
                        }
                        if (!temps[nameMap.get(namelist[i] + "-high")].split(":")[0].equals("./.")) {
                            count[i][1] += 1;
                        }
                        if (!temps[nameMap.get(namelist[i] + "-low")].split(":")[0].equals("./.") && !temps[nameMap.get(namelist[i] + "-high")].split(":")[0].equals("./.")) {
                            count[i][2] += 1;
                            if (!temps[nameMap.get(namelist[i] + "-low")].split(":")[0].equals(temps[nameMap.get(namelist[i] + "-high")].split(":")[0])) {
                                count[i][3] += 1;
                            }
                        }
                    }
                }
                br.close();
                bw.write("chr\ttaxa\tlow\thigh\tcomman\terror\n");
                for (int i = 0; i < namelist.length; i++) {
                    bw.write(f + "\t");
                    bw.write(namelist[i] + "\t");
                    for (int j = 0; j < 4; j++) {
                        bw.write(count[i][j] + "\t");
                    }
                    bw.newLine();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void metasoft(String[] infile) {
        new MetaTissue(infile);
    }

    public void simulationData(String[] args) {
        new simulationData(args);
    }

    public void getTransposonClass() {
        String gff3 = "/data1/home/xiaohan/Transposon/chr";
        String temp = null;
        String[] temps = null;
        HashSet<String> nameSet = new HashSet<>();
        try {
            for (int m = 0; m < 42; m++) {
                int chr = m + 1;
                BufferedReader br = IOUtils.getTextReader(new File(gff3, "chr" + chr + "_Transposon.gff3").getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    String[] tems = temps[8].split(";");
                    String classname = tems[3].split("=")[1].split("_")[0];
                    nameSet.add(classname);
                }
            }
            String[] names = nameSet.toArray(new String[nameSet.size()]);
            for (int i = 0; i < names.length; i++) {
                System.out.println(names[i]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void hapscanner(String[] args) throws IOException, InterruptedException {
        new HapscannerParameters(args);
//        new Hapscann(args);
    }

    public void getTransposonLength() {
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            String infile = "/data1/home/xiaohan/Transposon/chr/chr" + chr + "_Transposon.gff3";
            String outputDir = "/data1/home/xiaohan/Transposon/length/";
            BufferedReader br = IOUtils.getTextReader(infile);
            String temp = null;
            try {
                Multimap<String, String> TransMap = ArrayListMultimap.create();
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    String length = String.valueOf(Integer.parseInt(temp.split("\t")[4]) - Integer.parseInt(temp.split("\t")[3]));
                    String[] tems = temp.split("\t")[8].split(";");
                    String name = tems[3].split("=")[1].split("_")[0];
                    if (name.equals("RIX")) {
                        TransMap.put("DTC", length);
                    } else if (name.equals("RLC")) {
                        TransMap.put("DTH", length);
                    } else if (name.equals("RLG")) {
                        TransMap.put("DTM", length);
                    } else if (name.equals("RLX")) {
                        TransMap.put("DTT", length);
                    } else if (name.equals("XXX")) {
                        TransMap.put("XXX", length);
                    } else continue;
                }
                Collection<String> DTCmap = TransMap.get("DTC");
                String[] DTC = DTCmap.toArray(new String[DTCmap.size()]);
                Collection<String> DTHmap = TransMap.get("DTH");
                String[] DTH = DTHmap.toArray(new String[DTHmap.size()]);
                Collection<String> DTMmap = TransMap.get("DTM");
                String[] DTM = DTMmap.toArray(new String[DTMmap.size()]);
                Collection<String> DTTmap = TransMap.get("DTT");
                String[] DTT = DTTmap.toArray(new String[DTTmap.size()]);
                Collection<String> XXXmap = TransMap.get("XXX");
                String[] XXX = XXXmap.toArray(new String[XXXmap.size()]);
                BufferedWriter bwDTC = IOUtils.getTextWriter(new File(outputDir, chr + "RIX.txt").getAbsolutePath());
                BufferedWriter bwDTH = IOUtils.getTextWriter(new File(outputDir, chr + "RLC.txt").getAbsolutePath());
                BufferedWriter bwDTM = IOUtils.getTextWriter(new File(outputDir, chr + "RLG.txt").getAbsolutePath());
                BufferedWriter bwDTT = IOUtils.getTextWriter(new File(outputDir, chr + "RLX.txt").getAbsolutePath());
                BufferedWriter bwXXX = IOUtils.getTextWriter(new File(outputDir, chr + "XXX.txt").getAbsolutePath());
                for (int i = 0; i < DTC.length; i++) {
                    bwDTC.write(DTC[i]);
                    bwDTC.newLine();
                }
                bwDTC.flush();
                bwDTC.close();
                for (int i = 0; i < DTH.length; i++) {
                    bwDTH.write(DTH[i]);
                    bwDTH.newLine();
                }
                bwDTH.flush();
                bwDTH.close();
                for (int i = 0; i < DTM.length; i++) {
                    bwDTM.write(DTM[i]);
                    bwDTM.newLine();
                }
                bwDTM.flush();
                bwDTM.close();
                for (int i = 0; i < DTT.length; i++) {
                    bwDTT.write(DTT[i]);
                    bwDTT.newLine();
                }
                bwDTT.flush();
                bwDTT.close();
                for (int i = 0; i < XXX.length; i++) {
                    bwXXX.write(XXX[i]);
                    bwXXX.newLine();
                }
                bwXXX.flush();
                bwXXX.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }


    public void patternIdentify() throws IOException {
        String infileDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expression_pattern/Science";
        File[] fs = new File(infileDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "-region.txt");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName().split("region")[0];
            nameSet.add(Name);
//            System.out.println(Name);
        }
        BufferedWriter bw = IOUtils.getTextWriter("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expression_pattern/Science/expression_pattern.txt");
        nameSet.stream().forEach(f -> {
            try {
                BufferedReader br = IOUtils.getTextReader(new File(infileDir, f + "region.txt").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                int[] count = {0, 0, 0, 0, 0, 0, 0};
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("Chr")) continue;
                    if (!temps[2].equals("NaN") && !temps[5].equals("NaN") && !temps[8].equals("NaN")) {
                        String pattern = temps[9];
                        switch (pattern) {
                            case "M000":
                                count[0]++;
                                break;
                            case "M111":
                                count[0]++;
                                break;
                            case "M100":
                                count[1]++;
                                break;
                            case "M010":
                                count[2]++;
                                break;
                            case "M001":
                                count[3]++;
                                break;
                            case "M011":
                                count[4]++;
                                break;
                            case "M101":
                                count[5]++;
                                break;
                            case "M110":
                                count[6]++;
                                break;
                        }
                    }
                }
                bw.write(f + "\t");
                for (int i = 0; i < count.length; i++) {
                    bw.write(count[i] + "\t");
                }
                bw.newLine();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
        bw.flush();
        bw.close();
    }


    public void getTraidsPattern() {
        String ABDfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/TheABD_Science.txt";
        String exprfile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable/DEnorm7_87chr1-42_donor02.txt";
        String name = null;
        for (int m = 0; m < 87; m++) {
            int index = m + 4;
            HashMap<String, String> exprMap = new HashMap<>();
            HashSet<String> triadsSet = new HashSet<>();
            try {
                BufferedReader br = IOUtils.getTextReader(exprfile);
                String temp = null;
                String[] temps = null;
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("Chr")) {
                        name = temps[index];
                        continue;
                    }
                    String geneName = temps[3];
                    String expr = temps[index];
                    triadsSet.add(temps[3]);
                    exprMap.put(geneName, expr);
                }
                br.close();
                String outfile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expression_pattern/Science/" + name + ".txt";
                BufferedWriter bw = IOUtils.getTextWriter(outfile1);
                String triad = null;
                String[] triads = null;
                BufferedReader br1 = IOUtils.getTextReader(ABDfile);
                while ((triad = br1.readLine()) != null) {
                    triads = triad.split("\t");
                    if (triad.startsWith("A")) continue;
                    for (int i = 0; i < 3; i++) {
                        String geneName = triads[i];
                        if (triadsSet.contains(geneName)) {
                            String geneexpr = exprMap.get(geneName);
//                            System.out.println(geneexpr);
                            bw.write(geneName + "\t" + geneexpr + "\t");
                        } else {
                            bw.write(geneName + "\t" + 0 + "\t");
                        }
                    }
                    bw.newLine();
                }
                br1.close();
                bw.flush();
                bw.close();
                BufferedReader brpattern = IOUtils.getTextReader(outfile1);
                String pattern = null;
                String[] patterns = null;
                String outfile2 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expression_pattern/Science/" + name + "-region.txt";
                BufferedWriter bwpattern = IOUtils.getTextWriter(outfile2);
                while ((pattern = brpattern.readLine()) != null) {
                    patterns = pattern.split("\t");
                    double exprA = Double.parseDouble(patterns[1]);
                    double exprB = Double.parseDouble(patterns[3]);
                    double exprD = Double.parseDouble(patterns[5]);
                    double expr = exprA + exprB + exprD;
                    double ratioA = (double) exprA / expr;
                    double ratioB = (double) exprB / expr;
                    double ratioD = (double) exprD / expr;
                    double[] ratiodABD = {ratioA, ratioB, ratioD};
                    String region = Standardization.getNearestPointIndex(ratiodABD).getRegion();
                    StringBuilder sb = new StringBuilder();
                    sb.append(patterns[0]).append("\t").append(patterns[1]).append("\t").append(ratioA).append("\t");
                    sb.append(patterns[2]).append("\t").append(patterns[3]).append("\t").append(ratioB).append("\t");
                    sb.append(patterns[4]).append("\t").append(patterns[5]).append("\t").append(ratioD).append("\t");
                    sb.append(region);
                    bwpattern.write(sb.toString());
                    bwpattern.newLine();
                }
                brpattern.close();
                bwpattern.flush();
                bwpattern.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }


    public void intergenicPattern() {
//        String inputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/5nominalsThred/result/nominals.txt.gz";
//        String outputFileS="/data1/home/junxu/eQTL/FastQTL2/vst/nominals.distributation.txt";
//        GeneFeature gf = new GeneFeature("/data1/home/junxu/wheat_v1.1_Lulab.gff3");
        int subGenome = -1;
//        String pattern = null;
        String mode = null;
//        pattern = "ABD";
//        subGenome = 0;
//        mode = pattern + "/" + subGenome;
        mode = "5M/all";
//        String homogenefile = "/data1/home/xiaohan/rareallele/rankcorrelation/infor/TheABD.txt";
//        String inputFileS = "/data2/xiaohan/tensorQTL/6.txt";
//        String outputFileUp = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".up.distribution.txt";
//        String outputFileDown = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".down.distribution.txt";
//        String outputFileInter = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".inter.distribution.txt";
//        String outputFileUpEf = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".up.ef.txt";
//        String outputFileDownEf = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".down.ef.txt";
//        String outputFileInterEf = "/data2/xiaohan/tensorQTL/summary/countSig/" + mode + ".inter.ef.txt";
        String homogenefile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/TheABD.txt";
        String inputFileS = "/Users/yxh/Documents/eQTL/data_explain/5M/all.shuf.cis_qtl_pairs_sig.txt";
        String outputFileUp = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".up.distribution.txt";
        String outputFileDown = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".down.distribution.txt";
        String outputFileInter = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".inter.distribution.txt";
        String outputFileUpEf = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".up.ef.txt";
        String outputFileDownEf = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".down.ef.txt";
        String outputFileInterEf = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".inter.ef.txt";
        GeneFeature gf = new GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
//        int a=gf.getGeneStart(gf.getGeneIndex("TraesCS1A02G004400"));
//        int b=gf.getGeneEnd(gf.getGeneIndex("TraesCS1A02G004400"));
        int[] countUp = new int[5000];
        double[] efCUp = new double[5000];
        int up = 0;
        int[] countDown = new int[5000];
        double[] efCDown = new double[5000];
        int down = 0;
        int[] countInter = new int[100];
        double[] efCInter = new double[100];
        int Inter = 0;
        String temp = null;
        int pos = 0;
        int start = 0;
        int end = 0;
        int length = 0;
        String geneName = null;
        double ef = 0;
        try {
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bwUp = IOUtils.getTextWriter(outputFileUp);
            BufferedWriter bwDown = IOUtils.getTextWriter(outputFileDown);
            BufferedWriter bwInter = IOUtils.getTextWriter(outputFileInter);
            BufferedWriter bwUpEf = IOUtils.getTextWriter(outputFileUpEf);
            BufferedWriter bwDownEf = IOUtils.getTextWriter(outputFileDownEf);
            BufferedWriter bwInterEf = IOUtils.getTextWriter(outputFileInterEf);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Index")) continue;
                geneName = temp.split("\t")[1];
//                if (String.valueOf(geneName.charAt(8)).equals("D")) {
//                if(getHomoGene.ishomoGene(homogenefile,geneName).equals("A")){
//                if (temp.split("\t")[11].equals("nan")) {
//                    ef = 0;
//                } else {
//                    ef = Math.abs(Double.valueOf(temp.split("\t")[11]));
//                }
                int i = gf.getGeneIndex(geneName);
                if (gf.isWithinThisGene(i, Integer.valueOf(temp.split("\t")[2].split("_")[0]), Integer.valueOf(temp.split("\t")[2].split("_")[1]))) {
                    pos = Integer.valueOf(temp.split("\t")[2].split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {
                        start = gf.getGeneStart(i);
                        end = gf.getGeneEnd(i);
                        length = end - start;
                        int chunk = (pos - start) * 100 / length;
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    } else {
                        start = gf.getGeneEnd(i);
                        end = gf.getGeneStart(i);
                        length = start - end;
                        int chunk = (pos - end) * 100 / length;
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    }
                } else if (!gf.isWithinThisGene(i, Integer.valueOf(temp.split("\t")[2].split("_")[0]), Integer.valueOf(temp.split("\t")[2].split("_")[1]))) {
                    pos = Integer.valueOf(temp.split("\t")[2].split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {//1表示的是正链
                        start = gf.getGeneStart(i);
                        if (start >= pos) {
                            int chunk = (start - pos) / 1000;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = gf.getGeneEnd(i);
                            int chunk = (pos - end) / 1000;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    } else {
                        start = gf.getGeneEnd(i);
                        if (start <= pos) {
                            int chunk = (pos - start) / 1000;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = gf.getGeneStart(i);
                            int chunk = (end - pos) / 1000;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    }
                }
//                }
            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i;
                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
                if (efCUp[i] == 0) {
                    bwInterEf.write("Inter" + "\t" + chunk + "\t" + 0 + "\n");
                } else {
                    bwInterEf.write("Inter" + "\t" + chunk + "\t" + decFor.format((efCInter[i] / countInter[i]) * 1000000 / 1000000) + "\n");
                }
            }
            bwInter.flush();
            bwInter.close();
            bwInterEf.flush();
            bwInterEf.close();
            System.out.println(Inter);

            for (int i = 0; i < countUp.length; i++) {
                int chunk = i;
                bwUp.write("Up" + "\t" + chunk + "\t" + countUp[i] + "\n");
                if (efCUp[i] == 0) {
                    bwUpEf.write("Up" + "\t" + chunk + "\t" + 0 + "\n");
                } else {
                    bwUpEf.write("Up" + "\t" + chunk + "\t" + decFor.format((efCUp[i] / countUp[i]) * 1000000 / 1000000) + "\n");
                }
            }
            bwUp.flush();
            bwUp.close();
            bwUpEf.flush();
            bwUpEf.close();
            System.out.println(up);

            for (int i = 0; i < countDown.length; i++) {
                int chunk = i;
                bwDown.write("Down" + "\t" + chunk + "\t" + countDown[i] + "\n");
                if (efCDown[i] == 0) {
                    bwDownEf.write("Down" + "\t" + chunk + "\t" + 0 + "\n");
                } else {
                    bwDownEf.write("Down" + "\t" + chunk + "\t" + decFor.format((efCDown[i] / countDown[i]) * 1000000 / 1000000) + "\n");
                }
            }
            bwDown.flush();
            bwDown.close();
            bwDownEf.flush();
            bwDownEf.close();
            System.out.println(down);
        } catch (Exception ex) {
            System.out.println(geneName);
            System.out.println(pos);
            System.out.println(start);
//            ex.getStackTrace();
            ex.printStackTrace();
        }
    }

    public void intergenicPatternTransposon() {
        String infileDir = "/data1/home/xiaohan/Transposon/temp";
        String gff3 = "/data1/home/xiaohan/Transposon/chr";
        String outputdir = "/data1/home/xiaohan/Transposon/eQTLclassI";
        for (int m = 33; m < 42; m++) {
            int chr = m + 1;
            int size = 1000;
            int dis = 1000;
            int length = size * dis;
//        HashSet<String> nameSet = new HashSet();
//        String[] max = {"161314", "41550", "149688", "87337", "163527", "15883", "163284", "107539", "157986", "118253", "167585", "68613", "155416", "100464", "156750", "131745", "174312", "49935", "154526", "102576", "160736", "80592", "165370", "21126", "154742", "87533", "162763", "91419", "164112", "40868", "154323", "56209", "154150", "92592", "164392", "8422", "154761", "98233", "154984", "104129", "164596", "66656"};
//        HashMap<Integer, Integer> countMap = new HashMap<>();
//        for (int i = 0; i < 42; i++) {
//            int chr = i + 1;
//            countMap.put(chr, Integer.parseInt(max[i]));
//        }
//        for (int i = 1; i < 2; i++) {
//            String chr = String.valueOf(i);
//            nameSet.add(chr);
//        }
//        nameSet.stream().forEach(f -> {
            try {
                String temp = null;
                String[] temps = null;
                String temp1 = null;
                String[] temps1 = null;
                HashSet<String> DTCSet = new HashSet<>();
                HashSet<String> DTHSet = new HashSet<>();
                HashSet<String> DTMSet = new HashSet<>();
                HashSet<String> DTTSet = new HashSet<>();
                HashSet<String> DTXSet = new HashSet<>();
                BufferedReader br = IOUtils.getTextReader(new File(gff3, "chr" + chr + "_Transposon.gff3").getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    int strand = 0;
                    if (temps[6].equals("+")) {
                        strand = 1;
                    } else {
                        strand = 0;
                    }
                    String[] tems = temps[8].split(";");
//                    System.out.println(tems[3]);
                    if (tems[3].split("=")[1].split("_")[0].equals("RIX")) {
                        DTCSet.add(temps[3] + "_" + temps[4] + "_" + strand);
                    } else if (tems[3].split("=")[1].split("_")[0].equals("RLC")) {
                        DTHSet.add(temps[3] + "_" + temps[4] + "_" + strand);
                    } else if (tems[3].split("=")[1].split("_")[0].equals("RLG")) {
                        DTMSet.add(temps[3] + "_" + temps[4] + "_" + strand);
                    } else if (tems[3].split("=")[1].split("_")[0].equals("RLX")) {
                        DTTSet.add(temps[3] + "_" + temps[4] + "_" + strand);
                    } else if (tems[3].split("=")[1].split("_")[0].equals("XXX")) {
                        DTXSet.add(temps[3] + "_" + temps[4] + "_" + strand);
                    } else continue;
                }
                br.close();
                String[] DTC = DTCSet.toArray(new String[DTCSet.size()]);
                String[] DTH = DTHSet.toArray(new String[DTHSet.size()]);
                String[] DTM = DTMSet.toArray(new String[DTMSet.size()]);
                String[] DTT = DTTSet.toArray(new String[DTTSet.size()]);
                String[] DTX = DTXSet.toArray(new String[DTXSet.size()]);
                int[][] DTCup = new int[DTC.length][size];
                int[][] DTHup = new int[DTH.length][size];
                int[][] DTMup = new int[DTM.length][size];
                int[][] DTTup = new int[DTT.length][size];
                int[][] DTXup = new int[DTX.length][size];
                int[][] DTCdown = new int[DTC.length][size];
                int[][] DTHdown = new int[DTH.length][size];
                int[][] DTMdown = new int[DTM.length][size];
                int[][] DTTdown = new int[DTT.length][size];
                int[][] DTXdown = new int[DTX.length][size];
                int[][] DTCinter = new int[DTC.length][100];
                int[][] DTHinter = new int[DTH.length][100];
                int[][] DTMinter = new int[DTM.length][100];
                int[][] DTTinter = new int[DTT.length][100];
                int[][] DTXinter = new int[DTX.length][100];
                for (int i = 0; i < DTC.length; i++) {
                    for (int j = 0; j < size; j++) {
                        DTCup[i][j] = 0;
                        DTCdown[i][j] = 0;
                    }
                    for (int j = 0; j < 100; j++) {
                        DTCinter[i][j] = 0;
                    }
                }
                for (int i = 0; i < DTH.length; i++) {
                    for (int j = 0; j < size; j++) {
                        DTHup[i][j] = 0;
                        DTHdown[i][j] = 0;
                    }
                    for (int j = 0; j < 100; j++) {
                        DTHinter[i][j] = 0;
                    }
                }
                for (int i = 0; i < DTM.length; i++) {
                    for (int j = 0; j < size; j++) {
                        DTMup[i][j] = 0;
                        DTMdown[i][j] = 0;
                    }
                    for (int j = 0; j < 100; j++) {
                        DTMinter[i][j] = 0;
                    }
                }
                for (int i = 0; i < DTT.length; i++) {
                    for (int j = 0; j < size; j++) {
                        DTTup[i][j] = 0;
                        DTTdown[i][j] = 0;
                    }
                    for (int j = 0; j < 100; j++) {
                        DTTinter[i][j] = 0;
                    }
                }
                for (int i = 0; i < DTX.length; i++) {
                    for (int j = 0; j < size; j++) {
                        DTXup[i][j] = 0;
                        DTXdown[i][j] = 0;
                    }
                    for (int j = 0; j < 100; j++) {
                        DTXinter[i][j] = 0;
                    }
                }
                System.out.println("Finished initiating");
                System.out.println("DTC: " + DTC.length);
                System.out.println("DTH: " + DTH.length);
                System.out.println("DTM: " + DTM.length);
                System.out.println("DTT: " + DTT.length);
                System.out.println("DTX: " + DTX.length);
                Arrays.sort(DTC);
                Arrays.sort(DTH);
                Arrays.sort(DTM);
                Arrays.sort(DTT);
                Arrays.sort(DTX);
                BufferedReader br1 = IOUtils.getTextGzipReader(new File(infileDir, chr + ".cis_qtl_pairs_sig.txt.gz").getAbsolutePath());
                int countline = 0;
                while ((temp1 = br1.readLine()) != null) {
                    if (temp1.startsWith("Index")) continue;
                    countline++;
                    if (countline % 500 == 0) {
                        System.out.println(countline);
                    }
                    temps1 = temp1.split("\t");
                    int pos = Integer.parseInt(temps1[2].split("_")[1]);
                    for (int i = 0; i < DTC.length; i++) {
                        if (Integer.parseInt(DTC[i].split("_")[0]) == 1) {
                            int start = Integer.parseInt(DTC[i].split("_")[0]);
                            int end = Integer.parseInt(DTC[i].split("_")[1]);
                            if (pos < start && start - pos < length) {
                                int chunk = (start - pos) / dis;
                                DTCup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTCinter[i][chunk]++;
                            } else if (pos >= end && pos - end < length) {
                                int chunk = (pos - end) / dis;
                                DTCdown[i][chunk]++;
                            }
                        } else {
                            int end = Integer.parseInt(DTC[i].split("_")[0]);
                            int start = Integer.parseInt(DTC[i].split("_")[1]);
                            if (pos >= start && pos - start < length) {
                                int chunk = (pos - start) / dis;
                                DTCup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTCinter[i][chunk]++;
                            } else if (pos < end && end - pos < length) {
                                int chunk = (end - pos) / dis;
                                DTCdown[i][chunk]++;
                            }
                        }
                    }
                    for (int i = 0; i < DTH.length; i++) {
                        if (Integer.parseInt(DTH[i].split("_")[0]) == 1) {
                            int start = Integer.parseInt(DTH[i].split("_")[0]);
                            int end = Integer.parseInt(DTH[i].split("_")[1]);
                            if (pos < start && start - pos < length) {
                                int chunk = (start - pos) / dis;
                                DTHup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTHinter[i][chunk]++;
                            } else if (pos >= end && pos - end < length) {
                                int chunk = (pos - end) / dis;
                                DTHdown[i][chunk]++;
                            }
                        } else {
                            int end = Integer.parseInt(DTH[i].split("_")[0]);
                            int start = Integer.parseInt(DTH[i].split("_")[1]);
                            if (pos >= start && pos - start < length) {
                                int chunk = (pos - start) / dis;
                                DTHup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTHinter[i][chunk]++;
                            } else if (pos < end && end - pos < length) {
                                int chunk = (end - pos) / dis;
                                DTHdown[i][chunk]++;
                            }
                        }
                    }
                    for (int i = 0; i < DTM.length; i++) {
                        if (Integer.parseInt(DTM[i].split("_")[0]) == 1) {
                            int start = Integer.parseInt(DTM[i].split("_")[0]);
                            int end = Integer.parseInt(DTM[i].split("_")[1]);
                            if (pos < start && start - pos < length) {
                                int chunk = (start - pos) / dis;
                                DTMup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTMinter[i][chunk]++;
                            } else if (pos >= end && pos - end < length) {
                                int chunk = (pos - end) / dis;
                                DTMdown[i][chunk]++;
                            }
                        } else {
                            int end = Integer.parseInt(DTM[i].split("_")[0]);
                            int start = Integer.parseInt(DTM[i].split("_")[1]);
                            if (pos >= start && pos - start < length) {
                                int chunk = (pos - start) / dis;
                                DTMup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTMinter[i][chunk]++;
                            } else if (pos < end && end - pos < length) {
                                int chunk = (end - pos) / dis;
                                DTMdown[i][chunk]++;
                            }
                        }
                    }
                    for (int i = 0; i < DTT.length; i++) {
                        if (Integer.parseInt(DTT[i].split("_")[0]) == 1) {
                            int start = Integer.parseInt(DTT[i].split("_")[0]);
                            int end = Integer.parseInt(DTT[i].split("_")[1]);
                            if (pos < start && start - pos < length) {
                                int chunk = (start - pos) / dis;
                                DTTup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTTinter[i][chunk]++;
                            } else if (pos >= end && pos - end < length) {
                                int chunk = (pos - end) / dis;
                                DTTdown[i][chunk]++;
                            }
                        } else {
                            int end = Integer.parseInt(DTT[i].split("_")[0]);
                            int start = Integer.parseInt(DTT[i].split("_")[1]);
                            if (pos >= start && pos - start < length) {
                                int chunk = (pos - start) / dis;
                                DTTup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTTinter[i][chunk]++;
                            } else if (pos < end && end - pos < length) {
                                int chunk = (end - pos) / dis;
                                DTTdown[i][chunk]++;
                            }
                        }
                    }
                    for (int i = 0; i < DTX.length; i++) {
                        if (Integer.parseInt(DTX[i].split("_")[0]) == 1) {
                            int start = Integer.parseInt(DTX[i].split("_")[0]);
                            int end = Integer.parseInt(DTX[i].split("_")[1]);
                            if (pos < start && start - pos < length) {
                                int chunk = (start - pos) / dis;
                                DTXup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTXinter[i][chunk]++;
                            } else if (pos >= end && pos - end < length) {
                                int chunk = (pos - end) / dis;
                                DTXdown[i][chunk]++;
                            }
                        } else {
                            int end = Integer.parseInt(DTX[i].split("_")[0]);
                            int start = Integer.parseInt(DTX[i].split("_")[1]);
                            if (pos >= start && pos - start < length) {
                                int chunk = (pos - start) / dis;
                                DTXup[i][chunk]++;
                            } else if (pos >= start && pos < end) {
                                int chunk = (pos - start) * 100 / (end - start);
                                DTXinter[i][chunk]++;
                            } else if (pos < end && end - pos < length) {
                                int chunk = (end - pos) / dis;
                                DTXdown[i][chunk]++;
                            }
                        }
                    }
                }
                br1.close();
                System.out.println("Finished calculating");
                BufferedWriter bwDTC = IOUtils.getTextWriter(new File(outputdir, chr + "RIX.txt").getAbsolutePath());
                BufferedWriter bwDTH = IOUtils.getTextWriter(new File(outputdir, chr + "RLC.txt").getAbsolutePath());
                BufferedWriter bwDTM = IOUtils.getTextWriter(new File(outputdir, chr + "RLG.txt").getAbsolutePath());
                BufferedWriter bwDTT = IOUtils.getTextWriter(new File(outputdir, chr + "RLX.txt").getAbsolutePath());
                BufferedWriter bwDTX = IOUtils.getTextWriter(new File(outputdir, chr + "XXX.txt").getAbsolutePath());
                for (int i = 0; i < DTC.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    int start = Integer.parseInt(DTC[i].split("_")[0]);
                    int end = Integer.parseInt(DTC[i].split("_")[1]);
                    int translength = Math.abs(end - start);
                    sb.append(DTC[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(DTCup[i][j]).append("\t");
                    }
                    for (int j = 0; j < 100; j++) {
                        if (translength == 0) {
                            sb.append(0).append("\t");
                        } else {
                            int number = DTCinter[i][j] * 1000 / translength;
                            sb.append(number).append("\t");
                        }
                    }
                    for (int j = 0; j < size; j++) {
                        sb.append(DTCdown[i][j]).append("\t");
                    }
                    bwDTC.write(sb.toString());
                    bwDTC.newLine();
                }
                bwDTC.flush();
                bwDTC.close();
                for (int i = 0; i < DTH.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    int start = Integer.parseInt(DTH[i].split("_")[0]);
                    int end = Integer.parseInt(DTH[i].split("_")[1]);
                    int translength = Math.abs(end - start);
                    sb.append(DTH[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(DTHup[i][j]).append("\t");
                    }
                    for (int j = 0; j < 100; j++) {
                        if (translength == 0) {
                            sb.append(0).append("\t");
                        } else {
                            int number = DTHinter[i][j] * 1000 / translength;
                            sb.append(number).append("\t");
                        }
                    }
                    for (int j = 0; j < size; j++) {
                        sb.append(DTHdown[i][j]).append("\t");
                    }
                    bwDTH.write(sb.toString());
                    bwDTH.newLine();
                }
                bwDTH.flush();
                bwDTH.close();
                for (int i = 0; i < DTM.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    int start = Integer.parseInt(DTM[i].split("_")[0]);
                    int end = Integer.parseInt(DTM[i].split("_")[1]);
                    int translength = Math.abs(end - start);
                    sb.append(DTM[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(DTMup[i][j]).append("\t");
                    }
                    for (int j = 0; j < 100; j++) {
                        if (translength == 0) {
                            sb.append(0).append("\t");
                        } else {
                            int number = DTMinter[i][j] * 1000 / translength;
                            sb.append(number).append("\t");
                        }
                    }
                    for (int j = 0; j < size; j++) {
                        sb.append(DTMdown[i][j]).append("\t");
                    }
                    bwDTM.write(sb.toString());
                    bwDTM.newLine();
                }
                bwDTM.flush();
                bwDTM.close();
                for (int i = 0; i < DTT.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    int start = Integer.parseInt(DTT[i].split("_")[0]);
                    int end = Integer.parseInt(DTT[i].split("_")[1]);
                    int translength = Math.abs(end - start);
                    sb.append(DTT[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(DTTup[i][j]).append("\t");
                    }
                    for (int j = 0; j < 100; j++) {
                        if (translength == 0) {
                            sb.append(0).append("\t");
                        } else {
                            int number = DTTinter[i][j] * 1000 / translength;
                            sb.append(number).append("\t");
                        }
                    }
                    for (int j = 0; j < size; j++) {
                        sb.append(DTTdown[i][j]).append("\t");
                    }
                    bwDTT.write(sb.toString());
                    bwDTT.newLine();
                }
                bwDTT.flush();
                bwDTT.close();
                for (int i = 0; i < DTX.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    int start = Integer.parseInt(DTX[i].split("_")[0]);
                    int end = Integer.parseInt(DTX[i].split("_")[1]);
                    int translength = Math.abs(end - start);
                    sb.append(DTX[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(DTXup[i][j]).append("\t");
                    }
                    for (int j = 0; j < 100; j++) {
                        if (translength == 0) {
                            sb.append(0).append("\t");
                        } else {
                            int number = DTXinter[i][j] * 1000 / translength;
                            sb.append(number).append("\t");
                        }
                    }
                    for (int j = 0; j < size; j++) {
                        sb.append(DTXdown[i][j]).append("\t");
                    }
                    bwDTX.write(sb.toString());
                    bwDTX.newLine();
                }
                bwDTX.flush();
                bwDTX.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void intergenicPatternEnrichment() {
        int size = 100;
        int distance = 1000000;
        int bin = distance / size;
        String file = "1M";
        long startTime = System.currentTimeMillis();
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            String infile = "/data2/xiaohan/tensorQTL/" + "output" + "/" + chr + ".cis_qtl_pairs_sig.txt.gz";
            String outputDir = "/data2/xiaohan/tensorQTL/enrichment";
            String inputDir = "/data2/xiaohan/tensorQTL/tempvcf/";
            String vcfDir = "/data2/junxu/genotypeMaf005_87";
            BufferedReader br = IOUtils.getTextGzipReader(infile);
            BufferedReader br1 = IOUtils.getTextGzipReader(infile);
//            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, chr+"."+file+".enrichment.txt").getAbsolutePath());
            BufferedWriter bw1 = IOUtils.getTextWriter(new File(outputDir, chr + "." + file + ".cis_distance.txt").getAbsolutePath());
            BufferedWriter bw2 = IOUtils.getTextWriter(new File(outputDir, chr + "." + file + ".random_distance.txt").getAbsolutePath());
            String temp = null;
            String temp1 = null;
            String vcf = null;
            HashSet<String> geneSet = new HashSet<>();
            GeneFeature gf = new GeneFeature("/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3");
            try {
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Index")) continue;
                    geneSet.add(temp.split("\t")[1]);
                }
                br.close();
                String[] genelist = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genelist);
                int[][] random_distance = new int[genelist.length][size];
                int[][] cis_distance = new int[genelist.length][size];
                for (int i = 0; i < genelist.length; i++) {
                    for (int j = 0; j < size; j++) {
                        random_distance[i][j] = 0;
                        cis_distance[i][j] = 0;
                    }
                }
                HashMap<String, Integer> geneMap = new HashMap<>();
                for (int i = 0; i < genelist.length; i++) {
                    geneMap.put(genelist[i], i);
                }
                while ((temp1 = br1.readLine()) != null) {
                    if (temp1.startsWith("Index")) continue;
                    String geneName = temp1.split("\t")[1];
                    int pos = Integer.valueOf(temp1.split("\t")[2].split("_")[1]);
                    int index = gf.getGeneIndex(geneName);
                    int start = gf.getGeneStart(index);
                    int end = gf.getGeneEnd(index);
                    if (!gf.isWithinThisGene(index, chr, pos)) {
                        if (gf.getGeneStrand(index) == 1 && start >= pos) {//1表示的是正链
                            int chunk = (start - pos) / bin;
                            cis_distance[geneMap.get(geneName)][chunk]++;
                        } else if (gf.getGeneStrand(index) == 0 && end <= pos) {
                            int chunk = (pos - end) / bin;
                            cis_distance[geneMap.get(geneName)][chunk]++;
                        } else continue;
                    }
                }
                br1.close();
                for (int i = 0; i < genelist.length; i++) {
                    int index = gf.getGeneIndex(genelist[i]);
                    int start = gf.getGeneStart(index);
                    int end = gf.getGeneEnd(index);
                    int startsite = 0;
                    int endsite = 0;
                    String pos = null;
                    if (gf.getGeneStrand(index) == 1) {
                        startsite = start - distance + 1;
                        endsite = start;
                        if (startsite < 0) startsite = 0;
                        pos = chr + ":" + startsite + "-" + endsite;
                    } else {
                        startsite = end;
                        endsite = end + distance - 1;
                        pos = chr + ":" + startsite + "-" + endsite;
                    }
                    StringBuilder sb = new StringBuilder();
                    sb.append("tabix " + chr + ".87.B18.maf005.recode.vcf.gz " + pos + " > " + inputDir + "temp_enrichment.vcf ");
                    String command = sb.toString();
//                        System.out.println(command);
                    File dir = new File(new File(vcfDir).getAbsolutePath());
                    String[] cmdarry = {"/bin/bash", "-c", command};
                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                    p.waitFor();
                    BufferedReader brtemp = IOUtils.getTextReader(new File(inputDir, "temp_enrichment.vcf").getAbsolutePath());
                    while ((temp = brtemp.readLine()) != null) {
                        int snppos = Integer.parseInt(temp.split("\t")[1]);
                        if (snppos <= start) {
                            int chunk = (start - snppos) / bin;
                            if (chunk >= 100) {
                                System.out.println(chunk);
                                System.out.println(snppos);
                                System.out.println(start - snppos);
                                continue;
                            }
                            random_distance[geneMap.get(genelist[i])][chunk]++;
                        } else if (snppos >= end) {
                            int chunk = (snppos - end) / bin;
                            random_distance[geneMap.get(genelist[i])][chunk]++;
                        } else continue;
                    }
                    brtemp.close();
                    StringBuilder sb2 = new StringBuilder();
                    sb2.append("rm temp_enrichment.vcf ");
                    String command2 = sb2.toString();
                    File dir2 = new File(new File(inputDir).getAbsolutePath());
                    String[] cmdarry2 = {"/bin/bash", "-c", command2};
                    Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
                    p2.waitFor();
                    continue;
                }
                double[][] enrichment = new double[genelist.length][size];
                for (int i = 0; i < genelist.length; i++) {
                    for (int j = 0; j < size; j++) {
                        if (random_distance[i][j] == 0 && cis_distance[i][j] == 0) {
                            enrichment[i][j] = 0;
                        } else {
                            enrichment[i][j] = (double) cis_distance[i][j] / random_distance[i][j];
                        }
                    }
                }
                DecimalFormat decfor = new DecimalFormat("0.000");

                /*
                bw.write("Gene\t");
                for (int i = 0; i < size; i++) {
                    bw.write(i + "\t");
                }
                bw.newLine();
                for (int i = 0; i < genelist.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(genelist[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(decfor.format(enrichment[i][j])).append("\t");
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
                 */

                bw1.write("Gene\t");
                for (int i = 0; i < size; i++) {
                    bw1.write(i + "\t");
                }
                bw1.newLine();
                for (int i = 0; i < genelist.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(genelist[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(decfor.format(cis_distance[i][j])).append("\t");
                    }
                    bw1.write(sb.toString());
                    bw1.newLine();
                }
                bw1.flush();
                bw1.close();

                bw2.write("Gene\t");
                for (int i = 0; i < size; i++) {
                    bw2.write(i + "\t");
                }
                bw2.newLine();
                for (int i = 0; i < genelist.length; i++) {
                    StringBuilder sb = new StringBuilder();
                    sb.append(genelist[i]).append("\t");
                    for (int j = 0; j < size; j++) {
                        sb.append(decfor.format(random_distance[i][j])).append("\t");
                    }
                    bw2.write(sb.toString());
                    bw2.newLine();
                }
                bw2.flush();
                bw2.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
            long endTime = System.currentTimeMillis();
            System.out.println("程序运行时间： " + (endTime - startTime) + "ms");
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        new eQTL(args);
//        new eQTL();
    }
}
