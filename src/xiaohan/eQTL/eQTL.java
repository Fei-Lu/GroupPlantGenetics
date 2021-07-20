package xiaohan.eQTL;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import daxing.common.MD5;
import daxing.load.ancestralSite.Standardization;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import xiaohan.GBS.GBSsimulation;
import xiaohan.GBS.heterozygosity;
import xiaohan.eQTL.HapscanIdentifier.HapscannerParameters;
import xiaohan.eQTL.pipline.*;
import xiaohan.eQTL.simulation.simulationData;
import xiaohan.rareallele.GeneFeature;
import xiaohan.utils.IOUtils;
import xiaohan.rareallele.rareMutation;
import xiaohan.utils.SNPmappingInGene;

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
        effect size
         */
        this.effectsize(args);

        /*
        triads
         */
        this.getTraidsPattern();
        this.patternIdentify();

        /*
        VCF
         */
        this.vcf(args);

        /*
        pheno
         */
        this.pheno(args);
        this.covaraties(args);

        /*
        annotation
         */
        this.annotation(args);
        this.getExpressionWithaFC(args[0], args[1]);

        /*
        meta-tissue analysis
         */
        this.multiTissue(args);
//
        /*
        heterozygosity
         */
        this.heterozygosity(args);
        this.GBSsimulation(args);
//
        this.intergenicPattern();
        this.intergenicPatternEnrichment();
        this.intergenicPatternTransposon();
        this.getTransposonLength();
        this.getTransposonClass();

        this.ERCCRoc();
        this.ThreadPool(args[0],args[1]);
        this.subextractGeneRegion(args[0], args[1]);
        this.hetandMaf(args[0], args[1]);
        this.md5check();
        this.rareMutation(args);
        this.getsequencingDepth(args);
        this.getvariants(args);

    }


    public void getvariants(String[] args) {
        String inputDir = args[0];
        File[] fs = new File(new File(inputDir).getAbsolutePath()).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".txt");
        RowTable<String> rt1 = new RowTable<String>(args[1]);
        try {
            for (int i = 0; i < fs.length; i++) {
                RowTable<String> rt = new RowTable<String>(fs[i].getAbsolutePath());
                rt1.addColumn(fs[i].getName().replace("_counts_bygene.txt", ""), rt.getColumn(rt.getColumnIndex("n_variants")));
            }
            rt1.writeTextTable(args[2], IOFileFormat.Text);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void getsequencingDepth(String[] args) {
        String infile = args[0];
        String outfile = args[1];
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        BufferedReader br1 = IOUtils.getTextGzipReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        try {
            int[] count = null;
            String[] names = null;
            int[] depths = null;
            double[] avergeDepth = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    count = new int[temps.length - 9];
                    names = new String[temps.length - 9];
                    depths = new int[temps.length - 9];
                    avergeDepth = new double[temps.length - 9];
                    for (int i = 9; i < temps.length; i++) {
                        names[i - 9] = temps[i];
                        count[i - 9] = 0;
                        depths[i - 9] = 0;
                        avergeDepth[i - 9] = 0;
                    }
                    continue;
                }
                temps = temp.split("\t");
                for (int i = 9; i < temps.length; i++) {
                    if (!temps[i].startsWith("./.")) {
                        tems = temps[i].split(":");
                        depths[i - 9] += Integer.parseInt(tems[1].split(",")[0]) + Integer.parseInt(tems[1].split(",")[1]);
                        count[i - 9]++;
                    }
                }
            }
            br.close();
            for (int i = 0; i < avergeDepth.length; i++) {
                if (count[i] != 0) {
                    avergeDepth[i] = (double) depths[i] / count[i];
                } else {
                    avergeDepth[i] = -1;
                }

            }


            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < names.length; i++) {
                sb.append(names[i] + "\t");
            }
            bw.write(sb.toString().replaceAll("\\s+$", "")+"\n");

            sb.setLength(0);
            for (int i = 0; i < count.length; i++) {
                sb.append(count[i] + "\t");
            }
            bw.write(sb.toString().replaceAll("\\s+$", "")+"\n");

            sb.setLength(0);
            for (int i = 0; i < depths.length; i++) {
                sb.append(depths[i] + "\t");
            }
            bw.write(sb.toString().replaceAll("\\s+$", "")+"\n");

            DecimalFormat dec = new DecimalFormat("0.000");

            sb.setLength(0);
            for (int i = 0; i < avergeDepth.length; i++) {
                sb.append(dec.format(avergeDepth[i]) + "\t");
            }
            bw.write(sb.toString().replaceAll("\\s+$", "")+"\n");

            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void rareMutation(String[] args) {
        new rareMutation(args);
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

    public void ThreadPool(String infile, String thread) {
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        ExecutorService pool = Executors.newFixedThreadPool(Integer.parseInt(thread));
        File dir = new File(new File("/data1/home/xiaohan/jar").getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                String command = temp;
                Command com = new Command(command, dir);
                Future<Command> chrom = pool.submit(com);
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
                        double[] TPFPFN = simulationData.getPrecisionandRecall(expected, observed);
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

    public void heterozygosity(String[] args) throws IOException, InterruptedException {
        new heterozygosity(args);
    }

    public void GBSsimulation(String[] args) throws IOException, InterruptedException {
        new GBSsimulation(args);
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
