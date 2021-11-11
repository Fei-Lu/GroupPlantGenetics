package xiaohan.rareallele;

import pgl.infra.utils.PStringUtils;
import xiaohan.eQTL.GeneFeature;
import xiaohan.utils.IOUtils;
import xiaohan.utils.chrUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

/**
 * new scripts to analysis rare alleles
 *
 * @ author: yxh
 * @ created: 2021-10-24 : 9:56 AM
 */
public class NewStart {

    String input = null;
    String output = null;

    public NewStart(String[] args) {
//        this.mergeExpressionTable(args);
//        this.getbedFile(args);
//        this.splitBychr(args);
//        this.outliersSNP(args);
//        this.annotation(args);
//        this.summary(args);
        this.Version360(args);
//        this.getFileSplitInABD(args);
    }

    public void getFileSplitInABD(String[] args) {
        BufferedReader br;
        BufferedWriter[] bw = new BufferedWriter[3];
        br = IOUtils.getInFile(new File(args[0]).getAbsolutePath());
        for (int i = 0; i < bw.length; i++) {
            bw[i] = IOUtils.getOutFile(new File(args[i + 1]).getAbsolutePath());
        }
        int index = Integer.parseInt(args[args.length - 1]);
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (!temps[index].startsWith("Tr")) {
                    for (int i = 0; i < bw.length; i++) {
                        bw[i].write(temp + "\n");
                    }
                } else {
                    String CHR = temps[index].substring(8,9);
                    if(CHR.equals("A")){
                        bw[0].write(temp+"\n");
                    }else if (CHR.equals("B")){
                        bw[1].write(temp+"\n");
                    }else {
                        bw[2].write(temp+"\n");
                    }
                }
            }
            br.close();
            for (int i = 0; i < bw.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void Version360(String[] args) {
        String inputfile = new File(args[0]).getAbsolutePath();
        String outputfile = new File(args[1]).getAbsolutePath();
        BufferedReader br = IOUtils.getTextGzipReader(inputfile);
        BufferedWriter bw = IOUtils.getTextGzipWriter(outputfile);
        String temp = null;
        String[] temps = null;
        StringBuilder sb = new StringBuilder();
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) {
                    bw.write(temp + "\n");
                    continue;
                }
                temps = PStringUtils.fastSplit(temp).toArray(new String[0]);
                sb.setLength(0);
                for (int i = 0; i < temps.length; i++) {
                    if (i > 9 && (i - 8) % 24 == 0) {
                        if (temp.startsWith("#C")) {
                            System.out.println(temps[353]);
                            int sample = ((i - 8) / 24) * 25;
                            sb.append(temps[i].replace("B18-", "") + "\t" + "E" + PStringUtils.getNDigitNumber(3, sample) + "\t");
                        } else {
                            sb.append(temps[i] + "\t" + temps[353] + "\t");
                        }
                    } else {
                        if (temp.startsWith("#C")) {
                            sb.append(temps[i].replace("B18-", "") + "\t");
                        } else {
                            sb.append(temps[i] + "\t");
                        }
                    }
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void summary(String[] args) {
        BufferedReader br = IOUtils.getTextGzipReader(new File(args[0]).getAbsolutePath());
        BufferedReader br1 = IOUtils.getTextReader(new File(args[1]).getAbsolutePath());
        BufferedReader br2 = IOUtils.getTextReader(new File(args[2]).getAbsolutePath());
        BufferedWriter bw = IOUtils.getTextWriter(new File(args[3]).getAbsolutePath());
        double depth1 = 0;
        double depth2 = 0;
        String maf = null;
        String missingrate = null;
        String site = null;
        try {
            String temp = null;
            String[] temps = null;
            String temp1 = null;
            String[] temps1 = null;
            String temp2 = null;
            String[] temps2 = null;
            for (int i = 0; i < 15; i++) {
                temp = br.readLine();
            }
            temp1 = br1.readLine();
            temp2 = br2.readLine();
            StringBuilder sb = new StringBuilder();
            bw.write("Site\tdepthtotal\tdepthmean\tmaf\tmissingrate\n");
            while ((temp = br.readLine()) != null) {
                temp1 = br1.readLine();
                temp2 = br2.readLine();
                temps = temp.split("\t");
                temps1 = temp1.split("\t");
                temps2 = temp2.split("\t");
                site = temps[2];
                depth1 = depthSite(temps)[0];
                depth2 = depthSite(temps)[1];
                maf = MAF(temps1);
                missingrate = missingrate(temps2);
                sb.setLength(0);
                sb.append(site + "\t" + depth1 + "\t" + depth2 + "\t" + maf + "\t" + missingrate);
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static String missingrate(String[] temps) {
        return temps[5];
    }

    private static String MAF(String[] temps) {
        String maf = null;
        String[] tems = new String[2];
        tems[0] = temps[4].split(":")[1];
        tems[1] = temps[5].split(":")[1];
        if (tems[0].equals("-nan") || tems[1].equals("-nan")) {
            maf = "NA";
        } else {
            maf = String.valueOf(Math.min(Double.parseDouble(tems[0]), Double.parseDouble(tems[1])));
        }
        return maf;
    }

    private static double[] depthSite(String[] temps) {
        double depth = 0;
        double[] depths = new double[2];
        String[] tems = null;
        for (int i = 9; i < temps.length; i++) {
            if (!temps[i].startsWith("./.")) {
                tems = temps[i].split(":")[1].split(",");
                depth += Integer.parseInt(tems[0]);
                depth += Integer.parseInt(tems[1]);
            }
        }
        depths[0] = depth;
        depths[1] = depth / (temps.length - 9);
        return depths;
    }

    public void annotation(String[] args) {
        String input = args[0];
        String output = args[1];
        String[] category = {"three_prime_UTR", "exon", "CDS", "five_prime_UTR"};
        BufferedReader br = IOUtils.getTextGzipReader(new File(input).getAbsolutePath());
        BufferedWriter bw = IOUtils.getTextWriter(new File(output).getAbsolutePath());
        String temp = null;
        String[] temps = null;
        try {
            int countline = 0;
            temp = br.readLine();
            temps = temp.split("\t");
            int chr = Integer.parseInt(temps[0]);
            int start = Integer.parseInt(temps[1]);
            int end = Integer.parseInt(temps[2]);
            String anno = temps[3];
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 5000 == 0) {
//                    System.out.println(countline);
                }
                temps = temp.split("\t");
                if (temps[3].equals(anno)) {
                    if (Integer.parseInt(temps[1]) == end) {
                        end = Integer.parseInt(temps[2]);
                    } else {
                        bw.write("chr" + chr + "\t" + start + "\t" + end + "\t" + anno + "\n");
                        chr = Integer.parseInt(temps[0]);
                        start = Integer.parseInt(temps[1]);
                        end = Integer.parseInt(temps[2]);
                    }
                } else {
                    bw.write("chr" + chr + "\t" + start + "\t" + end + "\t" + anno + "\n");
                    chr = Integer.parseInt(temps[0]);
                    start = Integer.parseInt(temps[1]);
                    end = Integer.parseInt(temps[2]);
                    anno = temps[3];
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void outliersSNP(String[] args) {
        String input = args[0];
        String output = args[1];
        BufferedReader brinfo = null;
        BufferedReader br = null;
        BufferedWriter bw1 = null;
        BufferedWriter bw2 = null;
        BufferedWriter bw3 = null;
        String temp = null;
        String[] temps = null;
        String symbol = null;
        String site = null;
        String sample = null;
        String gene = null;
        int sampleindex = -1;
        int end = -1;
        int start = -1;
        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                brinfo = IOUtils.getTextReader(new File(input, chr + "_outliers_picked.txt").getAbsolutePath());
                br = IOUtils.getTextReader(new File(output, "chr" + PStringUtils.getNDigitNumber(3, chr) + "/all_rare_variants_chr" + PStringUtils.getNDigitNumber(3, chr) + "_SNPs.txt").getAbsolutePath());
                bw1 = IOUtils.getTextWriter(new File(output, "chr" + PStringUtils.getNDigitNumber(3, chr) + "/Underoutliers_chr" + PStringUtils.getNDigitNumber(3, chr) + "_SNPs.txt").getAbsolutePath());
                bw2 = IOUtils.getTextWriter(new File(output, "chr" + PStringUtils.getNDigitNumber(3, chr) + "/Overoutliers_chr" + PStringUtils.getNDigitNumber(3, chr) + "_SNPs.txt").getAbsolutePath());
                bw3 = IOUtils.getTextWriter(new File(output, "chr" + PStringUtils.getNDigitNumber(3, chr) + "/Nonoutliers_chr" + PStringUtils.getNDigitNumber(3, chr) + "_SNPs.txt").getAbsolutePath());
                HashSet<String> OveroutlierSet = new HashSet<>();
                HashSet<String> UnderoutlierSet = new HashSet<>();
                HashSet<String> NonoutlierSet = new HashSet<>();
                HashSet<String> siteSet = new HashSet<>();
                HashMap<String, String> zscoreMap = new HashMap<>();
                while ((temp = brinfo.readLine()) != null) {
                    if (temp.startsWith("GENE")) continue;
                    temps = temp.split("\t");
                    symbol = temps[0] + "_" + temps[1];
                    gene = temps[0];
                    sample = temps[1];
                    sampleindex = Integer.parseInt(temps[1].substring(1, 4));
                    if (Double.parseDouble(temps[3]) < 0) {
                        UnderoutlierSet.add(symbol);
                        zscoreMap.put(symbol, temps[3]);
                    }
                    if (Double.parseDouble(temps[3]) > 0) {
                        OveroutlierSet.add(symbol);
                        zscoreMap.put(symbol, temps[3]);
                    }
                    for (int j = 1; j <= 360; j++) {
                        if (sampleindex == j) continue;
                        sample = "E" + PStringUtils.getNDigitNumber(3, j);
                        symbol = gene + "_" + sample;
//                        System.out.println(sample);
                        NonoutlierSet.add(symbol);
                    }
                }
                brinfo.close();
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    symbol = temps[1] + "_" + temps[0];
                    end = Integer.parseInt(temps[3]);
                    start = end - 1;
                    site = temps[2] + "\t" + start + "\t" + end + "\t" + symbol + "\t" + zscoreMap.get(symbol);
                    if (UnderoutlierSet.contains(symbol)) {
                        bw1.write(site + "\n");
                    }
                    if (OveroutlierSet.contains(symbol)) {
                        bw2.write(site + "\n");
                    }
                    if (NonoutlierSet.contains(symbol)) {
                        bw3.write(site + "\n");
                    }
                }
                bw1.flush();
                bw1.close();
                bw2.flush();
                bw2.close();
                bw3.flush();
                bw3.close();
            }


        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void splitBychr(String[] args) {
        String input = args[0];
        String output = args[1];
        String suffix = args[2];
        GeneFeature gf = new GeneFeature("/data1/publicData/wheat/annotation/gene/v1.1/wheat_v1.1_Lulab.gff3");
        BufferedReader br = IOUtils.getTextReader(new File(input).getAbsolutePath());
        BufferedWriter[] bw = new BufferedWriter[42];
        for (int i = 0; i < bw.length; i++) {
            int chr = i + 1;
            bw[i] = IOUtils.getTextWriter(new File(output, chr + suffix).getAbsolutePath());
        }
        String[] temps = null;
        String temp = null;
        String gene = null;
        String subgenome = null;
        int chr = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("G")) {
                    for (int i = 0; i < bw.length; i++) {
                        bw[i].write(temp + "\n");
                    }
                    continue;
                }
                temps = temp.split("\t");
                gene = temps[0];
                chr = gf.getChromosomeOfGene(gf.getGeneIndex(gene));
                bw[chr - 1].write(temp + "\n");
            }
            for (int i = 0; i < bw.length; i++) {
                bw[i].flush();
                bw[i].close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getbedFile(String[] args) {
        String gff = args[0];
        int window = Integer.parseInt(args[1]);
        String outfile = args[2];
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        GeneFeature gf = new GeneFeature(gff);
        String[] genelist = new String[gf.getGeneNumber()];
//        System.out.println(gf.getGeneNumber());
        int start = 0;
        int end = 0;
        int chr = 0;
        try {
            for (int i = 0; i < genelist.length; i++) {
                System.out.println(i);
                genelist[i] = gf.getGeneName(i);
                chr = gf.getChromosomeOfGene(i);
                if (gf.getGeneStrand(i) == 1) {
                    start = gf.getGeneStart(i) - window - 1;
                    end = start + 1;
//                start = gf.getGeneStart(i);
//                end = gf.getGeneEnd(i);
                } else {
                    start = gf.getGeneEnd(i) + window - 1;
                    end = start + 1;
//                    start = gf.getGeneEnd(i);
//                    end = gf.getGeneStart(i);
                }
                if (start >= 0 && end >= 0) {
                    bw.write("chr" + chr + "\t" + start + "\t" + end + "\t" + genelist[i] + "\n");
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void mergeExpressionTable(String[] args) {
        this.input = args[0];
        this.output = args[1];
        BufferedWriter bw = IOUtils.getTextWriter(output);
        String temp = null;
        try {
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                BufferedReader br = IOUtils.getTextGzipReader(new File(input, chr + ".bed.gz").getAbsolutePath());
                while ((temp = br.readLine()) != null) {
                    if (temp.startsWith("#")) {
                        if (chr == 1) {
                            bw.write(temp);
                            bw.newLine();
                            continue;
                        } else continue;
                    }
                    bw.write(temp);
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {
        new NewStart(args);
    }
}