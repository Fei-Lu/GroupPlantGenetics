package xiaohan.eQTL;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import daxing.load.ancestralSite.Standardization;
import pgl.infra.table.RowTable;
import xiaohan.rareallele.GeneFeature;
import xiaohan.rareallele.IOUtils;
import xiaohan.rareallele.pheno;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class eQTL {

    int score = 0;
//    String infileDir = null;

    public eQTL(String[] args) throws IOException, InterruptedException {
         /*
//        simulation
//        */
//        this.simulationData(args);

        /*
//        Meta-Tissue Analysis pipeline
//         */
//        this.metasoft(args);

         /*
//        Hapscanner parameter
//         */
        this.hapscanner(args);

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

        /*
        meta-tissue analysis
         */
//        this.multiTissue(args);
//
//        this.hapscannercp(args);

//        this.intergenicPattern();
//        this.intergenicPatternEnrichment();
//        this.intergenicPatternhuman();
//        this.intergenicPatternTransposon();
//        this.getCisChrGene();
//        this.getNominalSig();
//        this.getcisSig();
//        this.getQTLcount();
//        this.getExprWithQTL();
//        this.getGeneUpSnpCount();
//        this.getGeneInterSnpCount();
//        this.getGeneDownSnpCount();
//        this.getRandomDistance();
//        this.getCisDistance();
//        this.getTransposonLength();
//        this.geteQTLpos();
//        this.getTransposonClass();
//        this.getoverlap();
//        this.snppos();
//        this.getsureofanno();
//        this.getsubGerp();
//        this.getNominalThreshold();
//        this.command();
    }

    public void pheno(String[] args){
        new pheno(args);
    }

    public void vcf(String[] args){
        new vcf(args);
    }

    public void effectsize(String[] args) {
        new effectsize(args);
    }

    public void multiTissue(String[] args) {
        new multiTissue(args);
    }

    public void getNominalThreshold() {
        String infileDir = "/data2/xiaohan/metasoft/v8/GTEx_Analysis_v8_eQTL";
        String infor = "/data2/xiaohan/metasoft/v8/samplesize.txt";
        File[] fs = new File(infileDir).listFiles();
        List<File> fList = new ArrayList(Arrays.asList());
        fs = IOUtils.listFilesEndsWith(fs, ".v8.egenes.txt.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
        }
        try {
            BufferedReader brinfo = IOUtils.getTextReader(infor);
            String temp1 = null;
            String[] temps1 = null;
            HashMap<String, Integer> TissueSize = new HashMap<>();
            while ((temp1 = brinfo.readLine()) != null) {
                if (temp1.startsWith("Tissue")) continue;
                temps1 = temp1.split("\t");
                System.out.println(temps1[0]);
                nameSet.add(temps1[0]);
                TissueSize.put(temps1[0], Integer.parseInt(temps1[1]));
            }
            brinfo.close();
            BufferedWriter bw = IOUtils.getTextWriter(new File(infileDir, "NominalT.txt").getAbsolutePath());
            bw.write("Group\tSampleSize\tNominalT\n");
            nameSet.stream().forEach(f -> {
                String temp = null;
                String[] temps = null;
                BufferedReader br = IOUtils.getTextGzipReader(new File(infileDir, f + ".v8.egenes.txt.gz").getAbsolutePath());
                String Tissue = String.valueOf(f);
                System.out.println(Tissue);
                int samplesize = TissueSize.get(Tissue);
                try {
                    while ((temp = br.readLine()) != null) {
                        if (temp.startsWith("gene")) continue;
                        temps = temp.split("\t");
                        double qval = Double.parseDouble(temps[28]);
                        double nomiT = 0.0000;
                        if (qval < 0.05) {
                            nomiT = Double.parseDouble(temps[29]);
                            bw.write(Tissue + "\t" + samplesize + "\t" + nomiT + "\n");
                        }
                    }
                    br.close();
                } catch (Exception e) {
                    e.printStackTrace();
                }
            });
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void command() {
        String[] sub = {"A", "B", "C", "D", "E", "F", "G", "H"};
        for (int i = 0; i < sub.length; i++) {
            for (int j = 1; j <= 9; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/14/sams/" + sub[i] + "0" + j + "_14_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/14sortedbam/" + sub[i] + "0" + j + "_14 && samtools index /data2/xiaohan/genotype/hapscanner/output/14sortedbam/" + sub[i] + "0" + j + "_14.bam");
            }
        }
        for (int i = 0; i < sub.length; i++) {
            for (int j = 10; j <= 12; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/14/sams/" + sub[i] + j + "_14_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/14sortedbam/" + sub[i] + j + "_14 && samtools index /data2/xiaohan/genotype/hapscanner/output/14sortedbam/" + sub[i] + j + "_14.bam");
            }
        }
        for (int i = 0; i < sub.length; i++) {
            for (int j = 1; j <= 9; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/29/sams/" + sub[i] + "0" + j + "_29_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/29sortedbam/" + sub[i] + "0" + j + "_29 && samtools index /data2/xiaohan/genotype/hapscanner/output/29sortedbam/" + sub[i] + "0" + j + "_29.bam");
            }
        }
        for (int i = 0; i < sub.length; i++) {
            for (int j = 10; j <= 12; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/29/sams/" + sub[i] + j + "_29_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/29sortedbam/" + sub[i] + j + "_29 && samtools index /data2/xiaohan/genotype/hapscanner/output/29sortedbam/" + sub[i] + j + "_29.bam");
            }
        }
        for (int i = 0; i < sub.length; i++) {
            for (int j = 1; j <= 9; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/46/sams/" + sub[i] + "0" + j + "_46_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/46sortedbam/" + sub[i] + "0" + j + "_46 && samtools index /data2/xiaohan/genotype/hapscanner/output/46sortedbam/" + sub[i] + "0" + j + "_46.bam");
            }
        }
        for (int i = 0; i < sub.length; i++) {
            for (int j = 10; j <= 12; j++) {
                System.out.println("samtools sort -n /data2/junxu/dataTest/test/46/sams/" + sub[i] + j + "_46_Aligned.out.bam /data2/xiaohan/genotype/hapscanner/output/46sortedbam/" + sub[i] + j + "_46 && samtools index /data2/xiaohan/genotype/hapscanner/output/46sortedbam/" + sub[i] + j + "_46.bam");
            }
        }
    }

    public void metasoft(String[] infile) {
        new MetaTissue(infile);
    }

    public void PlotPRplot(String infileDir) {
        StringBuilder sb = new StringBuilder();
        sb.append("R CMD BATCH /data1/home/xiaohan/simulationplot.R ");
        sb.append(new File(infileDir).getAbsolutePath());
        String command = sb.toString();
        System.out.println(command);
        try {
            File dir = new File(new File(infileDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished plot");
    }


    public void simulationData(String[] args) {
        new simulationData(args);
    }

    public void getsubGerp() {
        String infileDir = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/";
        String inputDir = "/data2/xiaohan/VCF_information/Gerpsnp";
        String outputDir = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/Gerpsnp";
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            BufferedReader br = IOUtils.getTextReader(new File(infileDir, chr + ".cis.sig.DE.log2.txt").getAbsolutePath());
            BufferedReader br1 = IOUtils.getTextGzipReader(new File(inputDir, chr + ".bed.gz").getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, chr + ".bed").getAbsolutePath());
            String temp = null;
            String[] temps = null;
            HashSet<String> posSet = new HashSet<>();
            try {
                System.out.println("Reading Gerp file chr" + chr + "...........");
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 5000 == 0) {
                        System.out.println(countline);
                    }
                    temps = temp.split("\t");
                    String pos = temps[1];
                    posSet.add(pos);
                }
                br.close();
                countline = 0;
                System.out.println("Reading Snp file chr" + chr + "...............");
                while ((temp = br1.readLine()) != null) {
                    countline++;
                    if (countline % 5000 == 0) {
                        System.out.println(countline);
                    }
                    temps = temp.split("\t");
                    String pos = temps[0] + "_" + temps[2];
                    double gerp = Double.parseDouble(temps[4]);
                    if (gerp < 0) continue;
                    if (posSet.contains(pos)) {
                        bw.write(temp + "\n");
                    }
                }
                br1.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void snppos() {
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            BufferedReader br = IOUtils.getTextGzipReader("/data2/junxu/genotypeMaf005_87/" + chr + ".87.B18.maf005.recode.vcf.gz");
            BufferedWriter bw = IOUtils.getTextWriter("/data2/xiaohan/VCF_information/snppos/" + chr + ".bed");
            String temp = null;
            String[] temps = null;
            try {
                System.out.println("Loading chr" + chr);
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    int end = Integer.parseInt(temps[1]);
                    int start = end - 1;
                    bw.write(temps[0] + "\t" + start + "\t" + end + "\t" + temps[3] + "\t" + temps[4] + "\n");
                }
                System.out.println(chr + "\t" + countline);
                System.out.println("Finish chr" + chr);
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getsureofanno() {
        String infile = "/data2/xiaohan/tensorQTL/001_exonSNP_anno.txt";
        String outfile = "/data2/xiaohan/tensorQTL/002_exonSNP_anno.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String mafdir = "/data2/xiaohan/VCF_information/snppos/";
        String tempdir = "/data2/xiaohan/tensorQTL/temp/";
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temp.startsWith("ID")) continue;
                if (temps[12].equals("NONCODING")) continue;
                String chr = temps[1];
                int end = Integer.parseInt(temps[2]);
                int start = end - 1;
                String ref = temps[3];
                String alt = temps[4];
                String pos = chr + ":" + start + "-" + end;
                StringBuilder sb = new StringBuilder();
                sb.append("tabix " + chr + ".bed.gz " + pos + " > " + tempdir + "temp.bed");
                String command = sb.toString();
                File dir = new File(new File(mafdir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
                BufferedReader brtemp = IOUtils.getTextReader(new File(tempdir, "temp.bed").getAbsolutePath());
                String temp1 = brtemp.readLine();
                if (temp1 != null) {
                    String[] temps1 = temp1.split("\t");
                    if (temps1[1].equals(String.valueOf(start))) {
                        if (temps1[3].equals(ref) && temps1[4].equals(alt)) {
                            bw.write(new StringBuilder(temps[1] + "\t" + temps[2] + "\t" + temps[3] + "\t" + temps[4] + "\t" + ref + "\t" + alt + "\t" + temps[12] + "\n").toString());
                        } else if (temps1[3].equals(alt) && temps1[4].equals(ref)) {
                            bw.write(new StringBuilder(temps[1] + "\t" + temps[2] + "\t" + temps[3] + "\t" + temps[4] + "\t" + ref + "\t" + alt + "\t" + temps[12] + "\n").toString());
                        } else {
                            System.out.println(chr + "\t" + pos + "\t" + ref + "\t" + alt + "\t" + temps1[3] + "\t" + temps1[4]);
                        }
                    }
                }
                brtemp.close();
                StringBuilder sb2 = new StringBuilder();
                sb2.append("rm temp.bed ");
                String command2 = sb2.toString();
//                        System.out.println(command2);
                File dir2 = new File(new File(tempdir).getAbsolutePath());
                String[] cmdarry2 = {"/bin/bash", "-c", command2};
                Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
                p2.waitFor();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }

    public void getoverlap() {
        String infile = "/data2/xiaohan/tensorQTL/annosite.txt";
        String infile2 = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect/all.cis.sig.DE.log2.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader br2 = IOUtils.getTextReader(infile2);
        String temp = null;
        String[] temps = null;
        HashSet<String> annosite = new HashSet<>();
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                for (int i = 0; i < temps.length; i++) {
                    annosite.add(temps[i]);
                }
            }
            br.close();
            while ((temp = br2.readLine()) != null) {
                temps = temp.split("\t");
                if (annosite.contains(temps[2])) {
                    System.out.println(temps[2]);
                }
            }
            br2.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
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

    public void gettransformedaFCfile() {
        String infile = "";
        String outfile = "";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        try {
            String temp = null;
            String[] temps = null;
            int logaFCindex = 0;
            int lowerindex = 0;
            int upperindex = 0;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("pheno") || temp.startsWith("pid") || temp.startsWith("Index")) {
                    bw.write(temp + "\n");
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals("log2_aFC")) logaFCindex = i;
                        if (temps[i].equals("log2_aFC_lower")) lowerindex = i;
                        if (temps[i].equals("log2_aFC_upper")) upperindex = i;
                    }
                    continue;
                }
                temps = temp.split("\t");
                if (temps[logaFCindex].equals("nan")) continue;
                if (temps[upperindex].equals("nan") || temps[lowerindex].equals("nan")) continue;
                double lower = Double.valueOf(temps[lowerindex]);
                double upper = Double.valueOf(temps[upperindex]);
                double ef = Math.abs(Double.valueOf(temps[logaFCindex]));
                if (ef > 2 / Math.log10(2)) continue;
                if (lower * upper < 0) continue;
                bw.write(temp + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    public void hapscannercp(String[] args){
        new HapScannercp(args[0]);
    }

    public void hapscanner(String[] args) throws IOException, InterruptedException {
        new HapscannerParameters(args);
    }

    public void geteQTLpos() {
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
//            File file = new File("/data1/home/xiaohan/tensorQTL/1M_log2/effectsize/" + chr + ".cis.sig.DE.log2.txt");
//            if (!file.exists()) continue;
            String infile = "/data1/home/xiaohan/tensorQTL/1M_log2/effectsize/" + chr + ".cis.sig.DE.log2.txt";
            String infor = "/data1/home/xiaohan/VCF_information/heter/7_87/site/chr" + chr + "heter.txt";
            String info = null;
            String[] infos = null;
            String outfile1 = "/data1/home/xiaohan/tensorQTL/1M_log2/heter/same" + chr + ".txt";
            String outfile2 = "/data1/home/xiaohan/tensorQTL/1M_log2/heter/dif" + chr + ".txt";
            String temp = null;
            String[] temps = null;
            HashMap<String, String> siteHeter = new HashMap<>();
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedReader br1 = IOUtils.getTextReader(infor);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfile2);
            try {
                while ((info = br1.readLine()) != null) {
                    infos = info.split("\t");
                    siteHeter.put(infos[0], infos[1]);
                }
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Index")) continue;
                    temps = temp.split("\t");
                    double slope = Double.parseDouble(temps[8]);
                    double aFC = Double.parseDouble(temps[11]);
                    if (slope * aFC > 0) {
                        bw1.write(temps[1] + "\t" + temps[2] + "\t" + temps[8] + "\t" + temps[11] + "\t" + siteHeter.get(temps[2]));
                        bw1.newLine();
                    } else if (slope * aFC < 0) {
                        bw2.write(temps[1] + "\t" + temps[2] + "\t" + temps[8] + "\t" + temps[11] + "\t" + siteHeter.get(temps[2]));
                        bw2.newLine();
                    } else continue;
                }
                br.close();
                bw1.flush();
                bw1.close();
                bw2.flush();
                bw2.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void calheterforSite() {
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            String infile = "/data2/junxu/genotypeMaf005_87/" + chr + ".87.B18.maf005.recode.vcf.gz";
            String outfile = "/data2/xiaohan/VCF_information/heter/7_87/site/chr" + chr + "heter.txt";
            BufferedReader br = IOUtils.getTextGzipReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            String temp = null;
            String[] temps = null;
            DecimalFormat decfor = new DecimalFormat("0.0000");
            try {
                while ((temp = br.readLine()) != null) {
                    int sum = 0;
                    int heter = 0;
                    double het = 0;
                    if (temp.startsWith("#")) continue;
                    temps = temp.split("\t");
                    for (int i = 9; i < temps.length; i++) {
                        String genotype = temps[i].split(":")[0];
//                    System.out.print(genotype);
                        if (genotype.equals("0/0") || genotype.equals("1/1")) {
                            sum += 1;
                        } else if (genotype.equals("0/1")) {
                            sum += 1;
                            heter += 1;
                        }
                    }
//                System.out.println(sum + "\t" + heter +"\t" + decfor.format(het));
                    het = (double) heter / sum;
//                het = heter/sum;
                    StringBuilder sb = new StringBuilder();
                    sb.append(temps[0] + "_" + temps[1] + "\t" + decfor.format(het));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
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

    public void getCisDistance() {
        HashSet<String> nameSet = new HashSet<>();
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            nameSet.add(String.valueOf(chr));
        }
        nameSet.parallelStream().forEach(chr -> {
            System.out.println(chr);
            String infile = "/data2/xiaohan/tensorQTL/5M/" + chr + ".pre.cis_qtl_pairs_sig.txt.gz";
            String outfile = "/data2/xiaohan/tensorQTL/snpcount/" + chr + ".5M.cis_distance.txt";
            BufferedReader br = IOUtils.getTextGzipReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            String temp = null;
            String[] temps = null;
            GeneFeature gf = new GeneFeature("/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3");
            try {
                int distance = 0;
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Index")) continue;
                    temps = temp.split("\t");
                    String geneName = temps[1];
                    String pos = temps[2].split("_")[1];
                    int position = Integer.parseInt(pos);
                    int index = gf.getGeneIndex(geneName);
                    int start = gf.getGeneStart(index);
                    int end = gf.getGeneEnd(index);
                    if (gf.getGeneStrand(gf.getGeneIndex(geneName)) == 1 && start >= position) {
                        distance = start - position;
                    } else if (gf.getGeneStrand(gf.getGeneIndex(geneName)) == 0 && end <= position) {
                        distance = position - end;
                    } else continue;
                    bw.write(geneName + "\t" + temp.split("\t")[2] + "\t" + distance);
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getRandomDistance() {
        int dis = 1000000;
        HashSet<String> nameSet = new HashSet<>();
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            nameSet.add(String.valueOf(chr));
        }
        nameSet.parallelStream().forEach(chr -> {
            String infile = "/data2/xiaohan/tensorQTL/output/" + chr + ".cis_qtl_pairs_sig.txt.gz";
            String outfile = chr + ".1M.random_distance.txt";
            String inputDir = "/data2/xiaohan/tensorQTL/tempvcf/";
            String outputDir = "/data2/xiaohan/tensorQTL/snpcount";
            String vcfdir = "/data2/junxu/genotypeMaf005_87";
            HashSet<String> geneSet = new HashSet<>();
            GeneFeature gf = new GeneFeature("/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3");
            try {
                String temp = null;
                BufferedReader br = IOUtils.getTextGzipReader(infile);
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, outfile).getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("Index")) continue;
                    String geneName = temp.split("\t")[1];
                    geneSet.add(geneName);
                }
                br.close();
                String[] genelist = geneSet.toArray(new String[geneSet.size()]);
                int startsite = 0;
                int endsite = 0;
                int start = 0;
                int end = 0;
                for (int i = 0; i < genelist.length; i++) {
                    String geneName = genelist[i];
                    start = gf.getGeneStart(gf.getGeneIndex(geneName));
                    end = gf.getGeneEnd(gf.getGeneIndex(geneName));
                    int strand = gf.getGeneStrand(gf.getGeneIndex(geneName));
                    if (strand == 1) {
                        startsite = start - dis + 1;
                        endsite = start;
                    } else {
                        startsite = end;
                        endsite = end + dis - 1;
                    }
                    if (startsite < 0) startsite = 0;
                    String pos = chr + ":" + startsite + "-" + endsite;
                    StringBuilder sb = new StringBuilder();
                    sb.append("tabix " + chr + ".87.B18.maf005.recode.vcf.gz " + pos + " > " + inputDir + "temp_" + chr + "1Mrandom.vcf ");
                    String command = sb.toString();
                    File dir = new File(new File(vcfdir).getAbsolutePath());
                    String[] cmdarry = {"/bin/bash", "-c", command};
                    Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                    p.waitFor();
                    BufferedReader brtemp = IOUtils.getTextReader(new File(inputDir, "temp_" + chr + "1Mrandom.vcf").getAbsolutePath());
                    String temp1 = null;
                    int countline = 0;
                    while ((temp1 = brtemp.readLine()) != null) {
                        countline++;
                        int distance = Math.abs(start - Integer.parseInt(temp1.split("\t")[1]));
                        bw.write(geneName + "\t" + temp1.split("\t")[2] + "\t" + distance);
                        bw.newLine();
                    }
                    brtemp.close();
                    StringBuilder sb2 = new StringBuilder();
                    sb2.append("rm temp_" + chr + "1Mrandom.vcf ");
                    String command2 = sb2.toString();
//                        System.out.println(command2);
                    File dir2 = new File(new File(inputDir).getAbsolutePath());
                    String[] cmdarry2 = {"/bin/bash", "-c", command2};
                    Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
                    p2.waitFor();
                    continue;
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getGeneInterSnpCount() {

    }

    public void getGeneDownSnpCount() {

    }

    public void getGeneUpSnpCount() throws IOException, InterruptedException {
        String inputDir = "/data2/xiaohan/tensorQTL/tempvcf/";
        String outputDir = "/data2/xiaohan/tensorQTL/snpcount";
        String vcfdir = "/data2/junxu/genotypeMaf005_87";
        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3");
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            String chr = String.valueOf(i + 1);
            nameSet.add(chr);
        }
        nameSet.parallelStream().forEach(f -> {
            try {
                String chr = f;
                String infile = "/data2/xiaohan/tensorQTL/bed/S7expression" + chr + ".bed";
                String[] genelist = utils.getGenelist(infile);
//            System.out.println(genelist[0]);
                int start = 0;
                int end = 0;
                String pos = null;
                int startsite = 0;
                int endsite = 0;
                for (int k = 1; k <= 5000; k++) {
                    int pos1 = k * 1000;
                    int pos2 = (k - 1) * 1000;
                    BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "chr" + chr + "_up_" + k + ".txt").getAbsolutePath());
                    for (int j = 0; j < genelist.length; j++) {
                        String geneName = genelist[j];
                        System.out.println(geneName);
                        if (gf.getGeneStrand(gf.getGeneIndex(geneName)) == 1) {
                            start = gf.getGeneStart(gf.getGeneIndex(geneName));
                            end = gf.getGeneEnd(gf.getGeneIndex(geneName));
                            startsite = start - pos1;
                            endsite = start - pos2;
                            pos = chr + ":" + startsite + "-" + endsite;
                            if (startsite < 0) {
                                startsite = 0;
                            }
                            StringBuilder sb = new StringBuilder();
                            sb.append("tabix " + chr + ".87.B18.maf005.recode.vcf.gz " + pos + " > " + inputDir + chr + "_temp_" + k + "_" + j + ".vcf ");
                            String command = sb.toString();
//                        System.out.println(command);
                            File dir = new File(new File(vcfdir).getAbsolutePath());
                            String[] cmdarry = {"/bin/bash", "-c", command};
                            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                            p.waitFor();
                            BufferedReader brtemp = IOUtils.getTextReader(new File(inputDir, chr + "_temp_" + k + "_" + j + ".vcf").getAbsolutePath());
                            String temp = null;
                            int countline = 0;
                            while ((temp = brtemp.readLine()) != null) {
                                countline++;
//                            System.out.println(countline);
                            }
                            brtemp.close();
                            bw.write(geneName + "\t" + countline);
                            bw.newLine();
                            StringBuilder sb2 = new StringBuilder();
                            sb2.append("rm " + chr + "_temp_" + k + "_" + j + ".vcf ");
                            String command2 = sb2.toString();
//                        System.out.println(command2);
                            File dir2 = new File(new File(inputDir).getAbsolutePath());
                            String[] cmdarry2 = {"/bin/bash", "-c", command2};
                            Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
                            p2.waitFor();
                            continue;
                        } else {
                            start = gf.getGeneStart(gf.getGeneIndex(geneName));
                            end = gf.getGeneEnd(gf.getGeneIndex(geneName));
                            startsite = start + pos2;
                            endsite = start + pos1;
                            pos = chr + ":" + startsite + "-" + endsite;
                            if (startsite < 0) {
                                startsite = 0;
                            }
                            StringBuilder sb = new StringBuilder();
                            sb.append("tabix " + chr + ".87.B18.maf005.recode.vcf.gz " + pos + " > " + inputDir + chr + "_temp_" + k + "_" + j + ".vcf ");
                            String command = sb.toString();
//                        System.out.println(command);
                            File dir = new File(new File(vcfdir).getAbsolutePath());
                            String[] cmdarry = {"/bin/bash", "-c", command};
                            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                            p.waitFor();
                            BufferedReader brtemp = IOUtils.getTextReader(new File(inputDir, chr + "_temp_" + k + "_" + j + ".vcf").getAbsolutePath());
                            String temp = null;
                            int countline = 0;
                            while ((temp = brtemp.readLine()) != null) {
                                countline++;
                            }
                            brtemp.close();
                            bw.write(geneName + "\t" + countline);
                            bw.newLine();
                            StringBuilder sb2 = new StringBuilder();
                            sb2.append("rm " + chr + "_temp_" + k + "_" + j + ".vcf ");
                            String command2 = sb2.toString();
//                        System.out.println(command2);
                            File dir2 = new File(new File(inputDir).getAbsolutePath());
                            String[] cmdarry2 = {"/bin/bash", "-c", command2};
                            Process p2 = Runtime.getRuntime().exec(cmdarry2, null, dir2);
                            p2.waitFor();
                            continue;
                        }
                    }
                    bw.flush();
                    bw.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getExprWithQTL() {
        String infileDir = "/data2/xiaohan/tensorQTL/1Mqtl/";
        String temp = null;
        String[] temps = null;
        String sample = "posgeneRegion";
        for (int m = 0; m < 87; m++) {
            int indi = m + 2;
            try {
                BufferedReader br = IOUtils.getTextGzipReader(new File(infileDir + sample + "/", "all.snp." + sample + ".txt.gz").getAbsolutePath());
                HashMap<String, Integer> geneCountMap = new HashMap<>();
                HashSet<String> geneSet = new HashSet<>();
                String index = null;
                temp = br.readLine();
                temps = temp.split("\t");
                if (temp.startsWith("Gene")) {
                    index = temps[indi];
                }
                System.out.println(index);
                temp = br.readLine();
                temps = temp.split("\t");
                String geneName = temps[0];
                System.out.println(geneName);
                geneSet.add(geneName);
                int count = Integer.parseInt(temps[indi]);
//                int countline = 0;
                while ((temp = br.readLine()) != null) {
//                    countline++;
//                    if (countline % 50000 == 0) {
//                        System.out.println(countline);
//                    }
                    temps = temp.split("\t");
                    if (temps[0].equals(geneName)) {
                        count += Integer.parseInt(temps[indi]);
//                        System.out.println(count);
                    } else {
                        geneCountMap.put(geneName, count);
                        geneName = temps[0];
                        geneSet.add(geneName);
                        count = Integer.parseInt(temps[indi]);
                    }
                }
                geneCountMap.put(geneName, count);
                br.close();
                String triad = null;
                String[] triads = null;
                BufferedReader br1 = IOUtils.getTextReader("/data2/xiaohan/tensorQTL/1Mqtl/Science/" + index + "-region.txt");
//                System.out.println("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expression_pattern/Science/" + index + "-region.txt");
                BufferedWriter bw = IOUtils.getTextWriter("/data2/xiaohan/tensorQTL/1Mqtl/geneRegion_output/" + index + "-region_count_" + sample + ".txt");
                bw.write("geneA\texprA\tratioA\tgeneB\texprB\tratioB\tgeneD\texprD\tratioD\tgroup\tcountA\tcountB\tcountD\tcountratioA\tcountratioB\tcountratioD\tcountGroup");
                bw.newLine();
                int countA = 0;
                int countB = 0;
                int countD = 0;
                while ((triad = br1.readLine()) != null) {
                    triads = triad.split("\t");
                    if (!triads[2].equals("NaN") && !triads[5].equals("NaN") && !triads[8].equals("NaN")) {
                        String geneA = triads[0];
                        if (geneSet.contains(geneA)) {
                            countA = geneCountMap.get(geneA);
                        } else {
                            countA = 0;
                        }
                        String geneB = triads[3];
                        if (geneSet.contains(geneB)) {
                            countB = geneCountMap.get(geneB);
                        } else {
                            countB = 0;
                        }
                        String geneD = triads[6];
                        if (geneSet.contains(geneD)) {
                            countD = geneCountMap.get(geneD);
                        } else {
                            countD = 0;
                        }
                        StringBuilder sb = new StringBuilder();
                        if (countA + countB + countD == 0) {
                            String region = "M000";
                            sb.append(triad).append("\t").append(countA).append("\t").append(countB).append("\t").append(countD).append("\t");
                            sb.append("0").append("\t").append("0").append("\t").append("0").append("\t");
                            sb.append(region);
                        } else {
                            double countratioA = (double) countA / (countA + countB + countD) / 2;
                            double countratioB = (double) countB / (countA + countB + countD) / 3;
                            double countratioD = (double) countD / (countA + countB + countD);
                            double[] ratiodABD = {countratioA, countratioB, countratioD};
                            String region = Standardization.getNearestPointIndex(ratiodABD).getRegion();
                            sb.append(triad).append("\t").append(countA).append("\t").append(countB).append("\t").append(countD).append("\t");
                            sb.append(countratioA).append("\t").append(countratioB).append("\t").append(countratioD).append("\t");
                            sb.append(region);
                        }
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                br1.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getQTLcount() {
        String inputDir = "/data2/junxu/genotypeMaf005_87";
        String inforDir = "/data2/xiaohan/tensorQTL/output";
        String outputDir = "/data2/xiaohan/tensorQTL/1Mqtl/neggeneRegion";
        File[] fs = new File(inputDir).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, "maf005.recode.vcf.gz");
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) {
                continue;
            }
            String Name = fs[i].getName().split("\\.")[0];
            nameSet.add(Name);
//            System.out.println(Name);
        }
        nameSet.stream().forEach(f -> {
            BufferedReader br = IOUtils.getTextGzipReader(new File(inputDir, f + ".87.B18.maf005.recode.vcf.gz").getAbsolutePath());
            BufferedReader brinfo = IOUtils.getTextGzipReader(new File(inforDir, f + ".cis_qtl_pairs_sig.txt.gz").getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, f + ".snp.neg.txt").getAbsolutePath());
            String info = null;
            String[] infos = null;
            String temp = null;
            String[] temps = null;
            HashSet<String> SNPset = new HashSet<>();
            Multimap<String, String> geneposMap = ArrayListMultimap.create();
            GeneFeature gf = new GeneFeature("/data1/home/xiaohan/rareallele/SiPASpipeline/reference/wheat_v1.1_Lulab.gff3");
            try {
                while ((info = brinfo.readLine()) != null) {
                    if (info.startsWith("Index")) continue;
                    infos = info.split("\t");
                    String snp = infos[2];
                    String gene = infos[1];
                    if (gf.isWithinThisGene(gf.getGeneIndex(gene), Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                        geneposMap.put(snp, gene);
                        double ef = Double.parseDouble(infos[8]);
                        if (ef > 0) {
                            SNPset.add(snp);
                        }
                    }
                }
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("##")) continue;
                    if (temp.startsWith("#C")) {
                        temps = temp.split("\t");
                        bw.write("Gene\tPOS\t");
                        for (int i = 9; i < temps.length; i++) {
                            bw.write(temps[i] + "\t");
                        }
                        bw.newLine();
                        continue;
                    }
                    temps = temp.split("\t");
                    if (SNPset.contains(temps[2])) {
                        Collection<String> genelist = geneposMap.get(temps[2]);
                        String[] genes = genelist.toArray(new String[genelist.size()]);
                        for (int j = 0; j < genes.length; j++) {
                            StringBuilder sb = new StringBuilder();
                            sb.append(genes[j]).append("\t");
                            sb.append(temps[2]).append("\t");
                            for (int i = 9; i < temps.length; i++) {
                                if (temps[i].split(":")[0].equals("0/0")) {
                                    sb.append("0\t");
                                } else if (temps[i].split(":")[0].equals("0/1")) {
                                    sb.append("1\t");
                                } else if (temps[i].split(":")[0].equals("1/1")) {
                                    sb.append("2\t");
                                } else {
                                    sb.append("0\t");
                                }
                            }
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                br.close();
                brinfo.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
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

    public void getcisSig() {
        for (int m = 0; m < 42; m++) {
            int chr = m + 1;
            String infile = "/data1/home/xiaohan/tensorQTL/1M_log2/" + chr + ".cis_qtl.txt.gz";
            String outfile = "/data1/home/xiaohan/tensorQTL/1M_log2/" + chr + ".top.cis_qtl.txt";
            BufferedReader br = IOUtils.getTextGzipReader(infile);
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            String temp = null;
            String[] temps = null;
            try {
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("phenotype_id")) {
                        temps[0] = "pid";
                        temps[6] = "sid";
                        StringBuilder sb = new StringBuilder();
                        sb.append(temps[0] + "\t" + temps[6] + "\t" + temps[7] + "\t" + temps[10] + "\t" + temps[13] + "\t" + temps[14]);
                        bw.write(sb.toString());
                        bw.newLine();
                        continue;
                    }
                    double qvalue = Double.parseDouble(temps[temps.length - 2]);
                    if (qvalue < 0.05) {
                        StringBuilder sb = new StringBuilder();
                        sb.append(temps[0] + "\t" + temps[6] + "\t" + temps[7] + "\t" + temps[10] + "\t" + temps[13] + "\t" + temps[14]);
                        bw.write(sb.toString());
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
    }

    public void getNominalSig() throws IOException {
        HashSet<String> nameSet = new HashSet<>();
        String geneName = null;
        String cis = null;
        String nominal = null;
        String[] cistemp = null;
        String[] nominaltemp = null;
        String index = null;
        double threshold = 0.00;
        double qvalue = 0.00;
        double pval = 0.00;
        String infileDir = "/data2/xiaohan/tensorQTL/Filter";
        String Dir = "/data2/xiaohan/tensorQTL/Filter";
        BufferedWriter bwerror = IOUtils.getTextWriter(new File(Dir, "error.txt").getAbsolutePath());
        try {
            for (int j = 0; j < 42; j++) {
                int f = j + 1;
                BufferedReader brcis = IOUtils.getTextGzipReader(new File(infileDir, f + ".cis_qtl.txt.gz").getAbsolutePath());
                BufferedReader brnominal = IOUtils.getTextGzipReader(new File(infileDir, f + ".cis_qtl_pairs." + f + ".txt.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(Dir, f + ".cis_qtl_pairs_sig.txt").getAbsolutePath());
                HashMap<String, Double> geneT = new HashMap<>();
                HashSet<String> eGene = new HashSet<>();
                while ((cis = brcis.readLine()) != null) {
                    if (cis.startsWith("p")) continue;
                    cistemp = cis.split("\t");
                    geneName = cistemp[0];
                    threshold = Double.parseDouble(cistemp[cistemp.length - 1]);
                    qvalue = Double.parseDouble(cistemp[cistemp.length - 2]);
                    if (qvalue < 0.05) {
                        eGene.add(geneName);
                        geneT.put(geneName, threshold);
                    }
                }
                brcis.close();
                while ((nominal = brnominal.readLine()) != null) {
                    nominaltemp = nominal.split("\t");
                    if (nominaltemp[1].startsWith("p")) {
                        bw.write("Index\tphenotype_id\tvariant_id\t");
                        for (int i = 3; i < nominaltemp.length; i++) {
                            bw.write(nominaltemp[i] + "\t");
                        }
                        bw.write("threshold");
                        bw.newLine();
                        continue;
                    }
                    geneName = nominaltemp[1];
                    if (!eGene.contains(nominaltemp[1])) continue;
                    threshold = geneT.get(nominaltemp[1]);
                    if (nominaltemp.length > 7) {
                        pval = Double.parseDouble(nominaltemp[7]);
                        index = nominaltemp[0];
                        if (pval < threshold) {
                            bw.write(nominal + "\t" + threshold);
                            bw.newLine();
                        }
                    }
                }
                geneT.clear();
                eGene.clear();
                brnominal.close();
                bw.flush();
                bw.close();

                StringBuilder sb1 = new StringBuilder();
                sb1.append("bgzip " + f + ".cis_qtl_pairs_sig.txt");
                String command = sb1.toString();
                File dir = new File(new File(Dir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
            System.out.println(geneName);
            System.out.println(pval);
            System.out.println(index);
            bwerror.write(nominal);
            bwerror.flush();
            bwerror.close();
        }
    }

    public void getCisChrGene() {
        RowTable rt = new RowTable("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable/DEnorm7_87chr1-42_donor02.txt");
        xujun.analysis.rnaseq.GeneFeature gf = new xujun.analysis.rnaseq.GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
//        String geneName=gf.getGeneName(19374);
        try {
            BufferedWriter[] bw = new BufferedWriter[42];
            for (int i = 0; i < bw.length; i++) {
                int chr = i + 1;
                bw[i] = pgl.infra.utils.IOUtils.getTextWriter("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/" + chr + ".txt");
                bw[i].write("Chr\tgene");
                bw[i].newLine();
            }
            for (int i = 0; i < rt.getRowNumber(); i++) {
                if (rt.getCell(i, 0).equals("Chr")) {
                    continue;
                }
                int index = gf.getGeneIndex(rt.getCell(i, 3).toString());
                int chro = gf.getGeneChromosome(index);
                System.out.println(rt.getCellAsString(i, 0) + "\t" + rt.getCellAsString(i, 3));
                bw[chro - 1].write(rt.getCellAsString(i, 0) + "\t" + rt.getCellAsString(i, 3));
                bw[chro - 1].newLine();
            }
            for (int i = 0; i < bw.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception ex) {
            ex.getStackTrace();
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

    public void intergenicPatternhuman() {
        String infileinfo = "/Users/yxh/Documents/eQTL/GTEx/Whole_Blood.v8.egenes.txt";
        String inputFileS = "/Users/yxh/Documents/eQTL/GTEx/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt";
        String outputFileUp = "/Users/yxh/Documents/eQTL/GTEx/up.distribution.txt";
//        String outputFileDown = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".down.distribution.txt";
//        String outputFileInter = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".inter.distribution.txt";
        String outputFileUpEf = "/Users/yxh/Documents/eQTL/GTEx/up.ef.txt";
//        String outputFileDownEf = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".down.ef.txt";
//        String outputFileInterEf = "/Users/yxh/Documents/eQTL/data_explain/" + mode + ".inter.ef.txt";
//        GeneFeature gf = new GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
//        int a=gf.getGeneStart(gf.getGeneIndex("TraesCS1A02G004400"));
//        int b=gf.getGeneEnd(gf.getGeneIndex("TraesCS1A02G004400"));
        int[] countUp = new int[1000];
        double[] efCUp = new double[1000];
        int up = 0;
//        int[] countDown = new int[100];
//        double[] efCDown = new double[100];
//        int down = 0;
//        int[] countInter = new int[1000];
//        double[] efCInter = new double[1000];
//        int Inter = 0;
        String temp = null;
        String info = null;
        int pos = 0;
        int start = 0;
        int end = 0;
        int length = 0;
        String geneName = null;
        double ef = 0;
        HashMap<String, String> geneStartMap = new HashMap<>();
        HashMap<String, String> geneEndMap = new HashMap<>();
        HashMap<String, String> geneStrandMap = new HashMap<>();
        try {
            BufferedReader brinfo = IOUtils.getTextReader(infileinfo);
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bwUp = IOUtils.getTextWriter(outputFileUp);
//            BufferedWriter bwDown = IOUtils.getTextWriter(outputFileDown);
//            BufferedWriter bwInter = IOUtils.getTextWriter(outputFileInter);
            BufferedWriter bwUpEf = IOUtils.getTextWriter(outputFileUpEf);
//            BufferedWriter bwDownEf = IOUtils.getTextWriter(outputFileDownEf);
//            BufferedWriter bwInterEf = IOUtils.getTextWriter(outputFileInterEf);
            while ((info = brinfo.readLine()) != null) {
                if (info.startsWith("gene")) continue;
                geneStartMap.put(info.split("\t")[0], info.split("\t")[3]);
                geneEndMap.put(info.split("\t")[0], info.split("\t")[4]);
                if (info.split("\t")[5].equals("+")) {
                    geneStrandMap.put(info.split("\t")[0], "1");
                } else {
                    geneStrandMap.put(info.split("\t")[0], "0");
                }
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("variant_id")) continue;
                if (temp.split("\t")[7].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Math.abs(Double.valueOf(temp.split("\t")[7]));
                }
                pos = Integer.valueOf(temp.split("\t")[0].split("_")[1]);
                if (geneStrandMap.get(temp.split("\t")[1]).equals("1")) {//1表示的是正链
                    start = Integer.parseInt(geneStartMap.get(temp.split("\t")[1]));
                    if (start >= pos) {
                        int chunk = (start - pos) / 1000;
                        countUp[chunk]++;
                        up++;
                        efCUp[chunk] += ef;
//                        } else {
//                            end = gf.getGeneEnd(i);
//                            int chunk = (pos - end) / 1000;
//                            countDown[chunk]++;
//                            down++;
//                            efCDown[chunk] += ef;
                    }
                } else {
                    start = Integer.parseInt(geneEndMap.get(temp.split("\t")[1]));
                    if (start <= pos) {
                        int chunk = (pos - start) / 1000;
                        countUp[chunk]++;
                        up++;
                        efCUp[chunk] += ef;
//                        } else {
//                            end = gf.getGeneStart(i);
//                            int chunk = (end - pos) / 1000;
//                            countDown[chunk]++;
//                            down++;
//                            efCDown[chunk] += ef;
//                        }
                    }
//                }
                }
            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
//            for (int i = 0; i < countInter.length; i++) {
//                int chunk = 10 * (i + 1) + 1000;
//                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
//                if (efCUp[i] == 0) {
//                    bwInterEf.write("Inter" + "\t" + chunk + "\t" + 0 + "\n");
//                } else {
//                    bwInterEf.write("Inter" + "\t" + chunk + "\t" + decFor.format((efCInter[i] / countInter[i]) * 1000000 / 1000000) + "\n");
//                }
//            }
//            bwInter.flush();
//            bwInter.close();
//            bwInterEf.flush();
//            bwInterEf.close();
//            System.out.println(Inter);

            for (int i = 0; i < countUp.length; i++) {
                int chunk = 1000 - i;
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

//            for (int i = 0; i < countDown.length; i++) {
//                int chunk = i + 1 + 2000;
//                bwDown.write("Down" + "\t" + chunk + "\t" + countDown[i] + "\n");
//                if (efCDown[i] == 0) {
//                    bwDownEf.write("Down" + "\t" + chunk + "\t" + 0 + "\n");
//                } else {
//                    bwDownEf.write("Down" + "\t" + chunk + "\t" + decFor.format((efCDown[i] / countDown[i]) * 1000000 / 1000000) + "\n");
//                }
//            }
//            bwDown.flush();
//            bwDown.close();
//            bwDownEf.flush();
//            bwDownEf.close();
//            System.out.println(down);
        } catch (Exception ex) {
            System.out.println(geneName);
            System.out.println(pos);
            System.out.println(start);
//            ex.getStackTrace();
            ex.printStackTrace();
        }
    }

    public static void main(String[] args) throws IOException, InterruptedException {
        new eQTL(args);
//        new eQTL();
    }
}
