package xiaohan.eQTL;

import smile.stat.Stat;
import xiaohan.rareallele.GeneFeature;
import xiaohan.rareallele.IOUtils;
import xuebo.analysis.annotation.BinarySearch;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;

public class annotation {

    String annotationfile = "/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3";
    String annotationfile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3";

    public annotation(String... args) {
        this.geneGerp();
//        this.genePhyloP();
    }

    public void geneGerp() {
        System.out.println("This program is running ********************************************************************");
        String inputdir = "/data2/xiaohan/annotation/GerpOrigin/chr/Conservation";
        String outputdir = "/data2/xiaohan/annotation/GerpOrigin/chr/Conservation/gene";
        System.out.println(inputdir);
        HashSet<String> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            chrSet.add(String.valueOf(chr));
        }
        GeneFeature gf = new GeneFeature(annotationfile);
        System.out.println("Finished genefeature ***********************************************************************");
        chrSet.parallelStream().forEach(f -> {
            try {
                HashSet<String> geneSet = new HashSet<>();
                for (int i = 0; i < gf.getGeneNumber(); i++) {
                    if (gf.getGeneChromosome(i) == Integer.parseInt(f)) {
                        geneSet.add(gf.getGeneName(i));
                    }
                }
                String[] genelist = geneSet.toArray(new String[geneSet.size()]);
                Arrays.sort(genelist);

                int[][] geneRange = new int[genelist.length][3];

                for (int i = 0; i < genelist.length; i++) {
                    geneRange[i][0] = gf.getGeneChromosome(gf.getGeneIndex(genelist[i]));
                    geneRange[i][1] = gf.getGeneStart(gf.getGeneIndex(genelist[i]));
                    geneRange[i][2] = gf.getGeneEnd(gf.getGeneIndex(genelist[i]));
                }

                int[] genelength = new int[genelist.length];
                double[] geneGerpsum = new double[genelist.length];
                double[] geneGerpbylength = new double[genelist.length];

                for (int i = 0; i < genelist.length; i++) {
                    genelength[i] = gf.getGeneEnd(gf.getGeneIndex(genelist[i])) - gf.getGeneStart(gf.getGeneIndex(genelist[i]));
                    geneGerpsum[i] = 0;
                }

                BufferedReader br = IOUtils.getTextGzipReader(new File(inputdir, "chr" + f + ".bed.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, "chr" + f + "geneGerpbylength.txt").getAbsolutePath());

                String temp = null;
                String[] temps = null;
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 500000 == 0) {
//                        System.out.println("Dealing line :" + countline);
                    }
                    temps = temp.split("\t");
//                    if (Double.parseDouble(temps[4]) <= 0) continue;
                    int pos = Integer.parseInt(temps[2]);
                    double gerp = Double.parseDouble(temps[4]);
//                    System.out.println(countline+"\t"+gerp);

                    int[] index = SNPmappingInGene.binarySearch(geneRange, pos);
                    for (int i = 0; i < index.length; i++) {
//                        System.out.println(index[i]);
                    }

                    if (index[0]!=-1) {
                        for (int i = 0; i < index.length; i++) {
                            geneGerpsum[index[i]] = geneGerpsum[index[i]] + gerp;
                        }
                    }
                }
                br.close();

                DecimalFormat decfor = new DecimalFormat("0.000000");

                for (int i = 0; i < genelist.length; i++) {
                    geneGerpbylength[i] = geneGerpsum[i]/genelength[i];
                }

                for (int i = 0; i < genelist.length; i++) {
                    bw.write(genelist[i] + "\t" + decfor.format(geneGerpbylength[i] * 1000000 / 1000000) + "\n");
                }
                bw.flush();
                bw.close();

            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void genePhyloP() {

    }

    public static void main(String[] args) {
        new annotation(args);
    }
}
