package xiaohan.eQTL.pipline;

import com.koloboke.collect.map.hash.HashIntIntMaps;
import pgl.infra.range.Range;
import xiaohan.rareallele.GeneFeature;
import xiaohan.utils.IOUtils;
import xiaohan.utils.SNPmappingInGene;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.text.DecimalFormat;
import java.util.*;

public class annotation {

    String annotationfile = "/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3";
    String annotationfile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3";

    public annotation(String[] args) {
//        this.geneGerp();
//        this.genePhyloP();
//        this.eQTLgerp(args[0], args[1]);
//        this.eQTLphyloP(args[0], args[1]);
//        this.getGeneannotation(args[0], args[1], args[2]);
        this.getWholeGenomeAnnotationMap();
    }

    public void getWholeGenomeAnnotationMap() {
//        String outputdir = "/Users/yxh/Documents/important/Lulab/exon_anno";
        String outputdir = "/data2/xiaohan/annotation/FCE";
        GeneFeature gf = new GeneFeature(annotationfile);
        System.out.println("Finished gathering gff3");
        HashSet<String> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            chrSet.add(String.valueOf(chr));
        }
//        BufferedReader brinfo = IOUtils.getTextGzipReader("/Users/yxh/Documents/important/Lulab/exon_anno/001_exonSNP_anno.txt.gz");
        BufferedReader brinfo = IOUtils.getTextGzipReader("/data2/xiaohan/annotation/exon_anno/001_exonSNP_anno.txt.gz");
        HashMap<String, String> annotationMap = new HashMap<>();
        HashSet<String> snpposition = new HashSet<>();
        String info = null;
        String[] infos = null;
        try {
            while ((info = brinfo.readLine()) != null) {
                if (info.startsWith("ID")) continue;
                infos = info.split("\t");
                snpposition.add(infos[0].split("-")[0] + "_" + infos[0].split("-")[1]);
                annotationMap.put(infos[0].split("-")[0] + "_" + infos[0].split("-")[1], infos[12]);
            }
            brinfo.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished gathering annotations");

        chrSet.stream().forEach(f -> {
            try {
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, "chr" + f + ".bed").getAbsolutePath());
                int chr = Integer.parseInt(f);

                HashSet<String> geneSet = new HashSet<>();
                for (int i = 0; i < gf.getGeneNumber(); i++) {
                    if (gf.getGeneChromosome(i) == chr) {
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
                System.out.println("Finished gathering geneRange");


                int chrStart = geneRange[0][1];
                int chrEnd = geneRange[geneRange.length - 1][2] + 1;
//                System.out.println(chrStart);
//                System.out.println(chrEnd);

                StringBuilder sb = new StringBuilder();
                for (int i = chrStart; i < chrEnd; i++) {
                    sb.setLength(0);
                    String annotation = null;
                    int startsite = i - 1;
                    int endsite = i;
                    String site = f + "_" + String.valueOf(i);
                    if (annotationMap.containsKey(site)) {
                        annotation = annotationMap.get(site);
                        sb.append(annotation);
                    } else {
                        int pos = i;
                        int[] index = SNPmappingInGene.binarySearch(geneRange, pos);
                        int finalindex = index[0];
                        if (index[0] == -1) {
                            annotation = "Intergenic";
                            sb.append(annotation);
                        } else {
                            if (index.length > 1) {
                                System.out.println(site);
                                for (int j = 0; j < index.length; j++) {
                                    System.out.println(genelist[index[j]]);
                                    if (finalindex > index[j]) {
                                        finalindex = index[j];
                                    }
                                }
                            }
                            Map<Integer, Integer> snpannotationMap = HashIntIntMaps.newMutableMap();
                            int geneindex = gf.getGeneIndex(genelist[finalindex]);
                            int t = gf.getLongestTranscriptIndex(geneindex);
                            if (gf.get5UTRList(geneindex, t).size() != 0) {
                                List<Range> fU = gf.get5UTRList(geneindex, t);
                                for (int a = 0; a < fU.size(); a++) {
                                    int start = fU.get(a).getRangeStart();
                                    int end = fU.get(a).getRangeEnd();
                                    for (int k = start; k <= end; k++) {
                                        snpannotationMap.put(k, 1);
                                    }
                                }
                            }
                            if (gf.getCDSList(geneindex, t).size() != 0) {
                                List<Range> CDS = gf.getCDSList(geneindex, t);
                                for (int a = 0; a < CDS.size(); a++) {
                                    int start = CDS.get(a).getRangeStart();
                                    int end = CDS.get(a).getRangeEnd();
                                    for (int k = start; k <= end; k++) {
                                        snpannotationMap.put(k, 2);
                                    }
                                }
                            }
                            if (gf.get3UTRList(geneindex, t).size() != 0) {
                                List<Range> tU = gf.get3UTRList(geneindex, t);
                                for (int a = 0; a < tU.size(); a++) {
                                    int start = tU.get(a).getRangeStart();
                                    int end = tU.get(a).getRangeEnd();
                                    for (int k = start; k <= end; k++) {
                                        snpannotationMap.put(k, 3);
                                    }
                                }
                            }
                            if (!snpannotationMap.containsKey(pos)) {
                                annotation = "Intron";
                            } else if (snpannotationMap.get(pos) == 1) {
                                annotation = "5-UTR";
                            } else if (snpannotationMap.get(pos) == 2) {
                                annotation = "CDS";
                            } else if (snpannotationMap.get(pos) == 3) {
                                annotation = "3-UTR";
                            }
                        }
                        if (sb.toString().replaceAll("\\s+$", "").contains("\t")) {
                            System.out.println(sb.toString());
                            if (sb.toString().contains("CDS")) {
                                annotation = "CDS";
                                continue;
                            } else if (sb.toString().contains("Intron")) {
                                annotation = "Intron";
                                continue;
                            } else if (sb.toString().contains("3-UTR")) {
                                annotation = "3-UTR";
                                continue;
                            } else if (sb.toString().contains("3-UTR")) {
                                annotation = "5-UTR";
                                continue;
                            }
                        }
                        sb.setLength(0);
                        sb.append(annotation);
                    }
                    bw.write(chr + "\t" + startsite + "\t" + endsite + "\t" + sb.toString() + "\n");
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public void getGeneannotation(String infile, String outfile, String chr) {
        BufferedReader br = null;
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(infile);
        } else {
            br = IOUtils.getTextReader(infile);
        }
        try {
            BufferedReader brinfo = IOUtils.getTextReader("/data2/xiaohan/tensorQTL/001_exonSNP_anno.txt");
            HashMap<String, String> annotationMap = new HashMap<>();
            HashSet<String> snpposition = new HashSet<>();
            String info = null;
            String[] infos = null;
            while ((info = brinfo.readLine()) != null) {
                if (info.startsWith("ID")) continue;
                infos = info.split("\t");
                snpposition.add(infos[0].split("-")[0] + "_" + infos[0].split("-")[1]);
                annotationMap.put(infos[0].split("-")[0] + "_" + infos[0].split("-")[1], infos[12]);
            }
            brinfo.close();

            GeneFeature gf = new GeneFeature(annotationfile);
            HashSet<String> geneSet = new HashSet<>();
            for (int i = 0; i < gf.getGeneNumber(); i++) {
                if (gf.getGeneChromosome(i) == Integer.parseInt(chr)) {
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

            String temp = null;
            String[] temps = null;
            int vindex = 0;
            int pindex = 0;
            StringBuilder sb = new StringBuilder();
            StringBuilder sb1 = new StringBuilder();

            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("p") || temp.startsWith("In")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals("variant_id")) {
                            vindex = i;
                        }
                        if (temps[i].equals("phenotype_id")) {
                            pindex = i;
                        }
                    }
                    sb.setLength(0);
                    sb.append("phenotype_id\tvariant_id\tannotation\n");
                    bw.write(sb.toString());
                    continue;
                }
                temps = temp.split("\t");

                sb1.setLength(0);
                String annotation = null;
                if (snpposition.contains(temps[vindex])) {
                    System.out.println(temps[vindex]);
                    System.out.println(annotationMap.get(temps[vindex]));
                    sb1.append(annotationMap.get(temps[vindex]));
                } else if (!snpposition.contains(temps[vindex])) {
                    int pos = Integer.parseInt(temps[vindex].split("_")[1]);
                    int[] index = SNPmappingInGene.binarySearch(geneRange, pos);
                    int finalindex = 0;
                    if (index.length != 1) {
                        for (int i = 0; i < index.length; i++) {
                            if (gf.isWithinThisGene(gf.getGeneIndex(genelist[index[i]]), Integer.parseInt(chr), pos)) {
                                finalindex = index[i];
                            } else {
                                finalindex = -1;
                            }
                        }
                    } else {
                        finalindex = index[0];
                    }

                    if (finalindex == -1) {
                        sb1.append("Intergenic");
                    } else {
//                        System.out.println(index[i]+"\t"+genelist[index[i]]+"\t");
                        int geneindex = gf.getGeneIndex(genelist[finalindex]);
                        int t = gf.getLongestTranscriptIndex(geneindex);
                        if (gf.get5UTRList(geneindex, t).size() != 0) {
                            List<Range> fU = gf.get5UTRList(geneindex, t);
                            for (int a = 0; a < fU.size(); a++) {
                                int start = fU.get(a).getRangeStart();
                                int end = fU.get(a).getRangeEnd();
                                if (pos >= start && pos <= end) {
                                    annotation = "5-UTR";
                                }
                            }
                            if (annotation == "5-UTR") {
                                sb1.append("5-UTR");
                            }
                        }
                        if (gf.getCDSList(geneindex, t).size() != 0) {
                            List<Range> CDS = gf.getCDSList(geneindex, t);
                            for (int a = 0; a < CDS.size(); a++) {
                                int start = CDS.get(a).getRangeStart();
                                int end = CDS.get(a).getRangeEnd();
                                if (pos >= start && pos <= end) {
                                    annotation = "3-UTR";
                                }
                            }
                            if (annotation == "3-UTR") {
                                sb1.append("3-UTR");
                            }
                        }
                        if (gf.get3UTRList(geneindex, t).size() != 0) {
                            List<Range> tU = gf.get3UTRList(geneindex, t);
                            for (int a = 0; a < tU.size(); a++) {
                                int start = tU.get(a).getRangeStart();
                                int end = tU.get(a).getRangeEnd();
                                if (pos >= start && pos <= end) {
                                    annotation = "CDS";
                                }
                            }
                            if (annotation == "CDS") {
                                sb1.append("CDS");
                            }
                        }
                        if (annotation == null) {
                            sb1.append("CDS");
                        }
                    }
                }
                sb.setLength(0);
                sb.append(temps[pindex] + "\t" + temps[vindex] + "\t" + sb1.toString());
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }

    public void eQTLphyloP(String infile, String outfile) {
        BufferedReader br = null;
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(infile);
        } else {
            br = IOUtils.getTextReader(infile);
        }
        try {
            String temp = null;
            String[] temps = null;
            int index = 0;
            int pindex = 0;
            DecimalFormat decfor = new DecimalFormat("0.000000");
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("p") || temp.startsWith("In")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals("variant_id")) {
                            index = i;
                        }
                        if (temps[i].equals("phenotype_id")) {
                            pindex = i;
                        }
                    }
                    sb.setLength(0);
                    sb.append("phenotype_id\tvariant_id\tPhyloP\n");
                    bw.write(sb.toString());
                    continue;
                }
                temps = temp.split("\t");
                int chr = Integer.parseInt(temps[index].split("_")[0]);
                int pos = Integer.parseInt(temps[index].split("_")[1]);
                double PhyloP = getPhyloP(chr, pos);
                sb.setLength(0);
                sb.append(temps[pindex] + "\t" + temps[index] + "\t" + decfor.format(PhyloP));
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

    public void eQTLgerp(String infile, String outfile) {
        BufferedReader br = null;
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(infile);
        } else {
            br = IOUtils.getTextReader(infile);
        }
        try {
            String temp = null;
            String[] temps = null;
            int index = 0;
            int pindex = 0;
            DecimalFormat decfor = new DecimalFormat("0.000000");
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("p") || temp.startsWith("In")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (temps[i].equals("variant_id")) {
                            index = i;
                        }
                        if (temps[i].equals("phenotype_id")) {
                            pindex = i;
                        }
                    }
                    System.out.print(index);
                    sb.setLength(0);
                    sb.append("phenotype_id\tvariant_id\tgerp\n");
                    bw.write(sb.toString());
                    continue;
                }
                temps = temp.split("\t");
                int chr = Integer.parseInt(temps[index].split("_")[0]);
                int pos = Integer.parseInt(temps[index].split("_")[1]);
                double gerp = getGerp(chr, pos);
                sb.setLength(0);
                sb.append(temps[pindex] + "\t" + temps[index] + "\t" + decfor.format(gerp));
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

    public static double getGerp(int chr, int pos) {
        String infileDir = "/data2/xiaohan/annotation/GerpOrigin/chr";
        String infile = new File(infileDir, "chr" + chr + ".bed.gz").getAbsolutePath();
        StringBuilder sb = new StringBuilder();
        sb.append("tabix ");
        sb.append(infile + " ");
        sb.append(chr);
        sb.append(":");
        sb.append(pos);
        sb.append("-");
        sb.append(pos);
        double gerp = 0.000;
        try {
            Process p = Runtime.getRuntime().exec(sb.toString());
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] temps = temp.split("\t");
                gerp = Double.parseDouble(temps[4]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return gerp;
    }

    public static double getPhyloP(int chr, int pos) {
        String infileDir = "/data2/xiaohan/annotation/PhyloP";
        String infile = new File(infileDir, "chr" + chr + ".bed.gz").getAbsolutePath();
        StringBuilder sb = new StringBuilder();
        sb.append("tabix ");
        sb.append(infile + " ");
        sb.append(chr);
        sb.append(":");
        sb.append(pos);
        sb.append("-");
        sb.append(pos);
        double PhyloP = 0.000;
        try {
            Process p = Runtime.getRuntime().exec(sb.toString());
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                String[] temps = temp.split("\t");
                PhyloP = Double.parseDouble(temps[3]);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return PhyloP;
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

                    if (index[0] != -1) {
                        for (int i = 0; i < index.length; i++) {
                            geneGerpsum[index[i]] = geneGerpsum[index[i]] + gerp;
                        }
                    }
                }
                br.close();

                DecimalFormat decfor = new DecimalFormat("0.000000");

                for (int i = 0; i < genelist.length; i++) {
                    geneGerpbylength[i] = geneGerpsum[i] / genelength[i];
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
