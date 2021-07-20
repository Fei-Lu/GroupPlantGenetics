package xiaohan.eQTL.pipline;

import pgl.infra.range.Range;
import xiaohan.rareallele.GeneFeature;
import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.text.DecimalFormat;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class effectsize {

    public effectsize(String[] args){

    }

    public void substractEffectsize() {
        String infile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/all.cis.sig.DE.log2.txt";
        String outfile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/all.cis.sig.DE.log2.clear.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temps[0].equals("pid")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                double ef = 0;
                if (temps[6].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Math.abs(Double.parseDouble(temps[6]));
                }
                if (ef >= (2 / Math.log10(2))) continue;
                if (temps[7].equals("nan") || temps[8].equals("nan")) continue;
                if (Double.valueOf(temps[7]) * Double.parseDouble(temps[8]) < 0) continue;
                bw.write(temp + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getDistanceEffectSizeTopupdate() {
        String mode = null;
//        mode = "all";
        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < ; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
        String f = "all";
//        nameSet.parallelStream().forEach(f -> {
        String Dir = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/";
//        String homogenefile = "/data2/xiaohan/reference/TheABD.txt";
        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = Dir + "distribution/" + f + ".up.distribution.txt";
        String outputFileDown = Dir + "distribution/" + f + ".down.distribution.txt";
        String outputFileInter = Dir + "distribution/" + f + ".inter.distribution.txt";
        String outputFileUpEf = Dir + "distribution/" + f + ".up.ef.txt";
        String outputFileDownEf = Dir + "distribution/" + f + ".down.ef.txt";
        String outputFileInterEf = Dir + "distribution/" + f + ".inter.ef.txt";
        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[1000];
        int[] countDown = new int[1000];
        int size = 100;
        int size1 = 1000 / size;
        int[] countInter = new int[size];
//        double[] efCInter = new double[size];
        for (int i = 0; i < 100; i++) {
            countDown[i] = 0;
            countUp[i] = 0;
        }
        for (int i = 0; i < size; i++) {
            countInter[i] = 0;
//            efCInter[i] = 0;
        }
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
                if (temp.startsWith("pid")) continue;
                String[] temps = temp.split("\t");
                geneName = temps[0];
                String snp = temps[1];
                if (temps[6].equals("nan")) {
//                    ef = 0;
                    continue;
                } else {
                    ef = Math.abs(Double.valueOf(temps[6]));
                }
                if (ef > 2 / Math.log10(2)) continue;
                if (temps[7].equals("nan") || temps[8].equals("nan")) continue;
                if (Double.valueOf(temps[7]) * Double.parseDouble(temps[8]) < 0) continue;
                int i = gf.getGeneIndex(geneName);
                if (gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {
                        start = gf.getGeneStart(i);
                        end = gf.getGeneEnd(i);
                        length = end - start;
                        int chunk = (pos - start - 1) * size / length;
                        double chunk1 = (pos - start - 1) * size / length;
                        chunk1 = 1000 + chunk1 * 10;
                        bwInterEf.write(chunk1 + "\t" + ef + "\n");
                    } else {
                        start = gf.getGeneEnd(i);
                        end = gf.getGeneStart(i);
                        length = start - end;
                        int chunk = (start - pos - 1) * size / length;
                        double chunk1 = (start - pos - 1) * size / length;
                        countInter[chunk]++;
                        chunk1 = 1000 + chunk1 * 10;
                        bwInterEf.write(chunk1 + "\t" + ef + "\n");
                    }
                } else if (!gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {//1表示的是正链
                        start = gf.getGeneStart(i);
                        if (start >= pos) {
                            int chunk = (start - pos - 1) / 100;
                            double chunk1 = (start - pos - 1) / 100;
                            if (chunk <= 999) {
                                chunk1 = 1000 - chunk1;
                                bwUpEf.write(chunk1 + "\t" + ef + "\n");
                                countUp[chunk]++;
                            }

                        } else {
                            end = gf.getGeneEnd(i);
                            int chunk = (pos - end - 1) / 100;
                            double chunk1 = (pos - end - 1) / 100;
                            if (chunk <= 999) {
                                chunk1 = chunk1 + 2001;
                                bwDownEf.write(chunk1 + "\t" + ef + "\n");
                                countDown[chunk]++;
                            }
                        }
                    } else {
                        start = gf.getGeneEnd(i);
                        if (start <= pos) {
                            int chunk = (pos - start - 1) / 100;
                            double chunk1 = (pos - start - 1) / 100;
                            if (chunk <= 999) {
                                chunk1 = 1000 - chunk1;
                                bwUpEf.write(chunk1 + "\t" + ef + "\n");
                                countUp[chunk]++;
                            }
                        } else {
                            end = gf.getGeneStart(i);
                            int chunk = (end - pos - 1) / 100;
                            double chunk1 = (end - pos - 1) / 100;
                            if (chunk <= 999) {
                                chunk1 = chunk1 + 2001;
                                bwDownEf.write(chunk1 + "\t" + ef + "\n");
                                countDown[chunk]++;
                            }
                        }
                    }
                }
            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i * size1 + 1001;
                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
            }
            bwInter.flush();
            bwInter.close();
            bwInterEf.flush();
            bwInterEf.close();

            for (int i = 0; i < countUp.length; i++) {
                int chunk = 1000 - i;
                bwUp.write("Up" + "\t" + chunk + "\t" + countUp[i] + "\n");
            }
            bwUp.flush();
            bwUp.close();
            bwUpEf.flush();
            bwUpEf.close();

            for (int i = 0; i < countDown.length; i++) {
                int chunk = i + 2001;
                bwDown.write("Down" + "\t" + chunk + "\t" + countDown[i] + "\n");
            }
            bwDown.flush();
            bwDown.close();
            bwDownEf.flush();
            bwDownEf.close();
        } catch (Exception ex) {
            System.out.println(geneName);
            System.out.println(pos);
            System.out.println(start);
            ex.printStackTrace();
        }
//        });
    }

    public void getDistanceEffectSizeTop() {
        String mode = null;
//        mode = "all";
        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < ; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
        String infor = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/genelist/donor02GeneName.txt";
        String f = "all";
//        nameSet.parallelStream().forEach(f -> {
        String Dir = "/data2/xiaohan/tensorQTL/updated_1Mlog2/homoeffect_top/";
//        String homogenefile = "/data2/xiaohan/reference/TheABD.txt";
        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = Dir + "distribution/" + f + ".up.distribution.txt";
        String outputFileDown = Dir + "distribution/" + f + ".down.distribution.txt";
        String outputFileInter = Dir + "distribution/" + f + ".inter.distribution.txt";
        String outputFileUpEf = Dir + "distribution/" + f + ".up.ef.txt";
        String outputFileDownEf = Dir + "distribution/" + f + ".down.ef.txt";
        String outputFileInterEf = Dir + "distribution/" + f + ".inter.ef.txt";
        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[100];
        double[] efCUp = new double[100];
        int up = 0;
        int[] countDown = new int[100];
        double[] efCDown = new double[100];
        int down = 0;
        int size = 10;
        int size1 = 100 / size;
        int[] countInter = new int[size];
        double[] efCInter = new double[size];
        for (int i = 0; i < 100; i++) {
            countDown[i] = 0;
            countUp[i] = 0;
            efCDown[i] = 0;
            efCUp[i] = 0;
        }
        for (int i = 0; i < size; i++) {
            countInter[i] = 0;
            efCInter[i] = 0;
        }
        int Inter = 0;
        String temp = null;
        int pos = 0;
        int start = 0;
        int end = 0;
        int length = 0;
        String geneName = null;
        double ef = 0;
        try {
            HashSet<String> geneSet = new HashSet<>();
            BufferedReader brinfo = IOUtils.getTextReader(infor);
            while ((temp = brinfo.readLine()) != null) {
                geneSet.add(temp);
            }
            BufferedReader br = IOUtils.getTextReader(inputFileS);
            BufferedWriter bwUp = IOUtils.getTextWriter(outputFileUp);
            BufferedWriter bwDown = IOUtils.getTextWriter(outputFileDown);
            BufferedWriter bwInter = IOUtils.getTextWriter(outputFileInter);
            BufferedWriter bwUpEf = IOUtils.getTextWriter(outputFileUpEf);
            BufferedWriter bwDownEf = IOUtils.getTextWriter(outputFileDownEf);
            BufferedWriter bwInterEf = IOUtils.getTextWriter(outputFileInterEf);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("pid")) continue;
                String[] temps = temp.split("\t");
                geneName = temps[0];
                if (!geneSet.contains(geneName)) continue;
                String snp = temps[6];
//                if (String.valueOf(geneName.charAt(8)).equals("D")) {
//                if(getHomoGene.ishomoGene(homogenefile,geneName).equals("A")){
                if (temps[19].equals("nan")) {
//                    ef = 0;
                    continue;
                } else {
                    ef = Math.abs(Double.valueOf(temps[19]));
                }
                if (ef > 2 / Math.log10(2)) continue;
                if (temps[20].equals("nan") || temps[21].equals("nan")) continue;
                if (Double.valueOf(temps[20]) * Double.parseDouble(temps[21]) < 0) continue;
                int i = gf.getGeneIndex(geneName);
                if (gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {
                        start = gf.getGeneStart(i);
                        end = gf.getGeneEnd(i);
                        length = end - start;
                        int chunk = (pos - start - 1) * size / length;
//                        int chunk = (pos - start - 1) / 100;
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    } else {
                        start = gf.getGeneEnd(i);
                        end = gf.getGeneStart(i);
                        length = start - end;
                        int chunk = (start - pos - 1) * size / length;
//                        int chunk = (start - pos - 1) / 100 ;
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    }
                } else if (!gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {//1表示的是正链
                        start = gf.getGeneStart(i);
                        if (start >= pos) {
                            int chunk = (start - pos - 1) / 1000;
                            if (chunk <= 99) {
                                countUp[chunk]++;
                                up++;
                                efCUp[chunk] += ef;
                            }
                        } else {
                            end = gf.getGeneEnd(i);
                            int chunk = (pos - end - 1) / 1000;
                            if (chunk <= 99) {
                                countDown[chunk]++;
                                down++;
                                efCDown[chunk] += ef;
                            }
                        }
                    } else {
                        start = gf.getGeneEnd(i);
                        if (start <= pos) {
                            int chunk = (pos - start - 1) / 1000;
                            if (chunk <= 99) {
                                countUp[chunk]++;
                                up++;
                                efCUp[chunk] += ef;
                            }
                        } else {
                            end = gf.getGeneStart(i);
                            int chunk = (end - pos - 1) / 100;
                            if (chunk <= 99) {
                                countDown[chunk]++;
                                down++;
                                efCDown[chunk] += ef;
                            }
                        }
                    }
                }
            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i * size1 + 101;
                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
                if (efCInter[i] == 0) {
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
                int chunk = 100 - i;
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
                int chunk = i + 201;
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
            ex.printStackTrace();
        }
    }


    public void getDistanceEffectSizeShuf() {
//        String mode = null;
//        mode = "all";
//        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < 42; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
//        nameSet.parallelStream().forEach(f -> {
        String f = "all";
        String Dir = "/data1/home/xiaohan/tensorQTL/1M_log2/effectsize/";
        String homogenefile = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/TheABD.txt";
        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = Dir + "distribution/" + f + ".up.distribution.txt";
        String outputFileDown = Dir + "distribution/" + f + ".down.distribution.txt";
        String outputFileInter = Dir + "distribution/" + f + ".inter.distribution.txt";
        String outputFileUpEf = Dir + "distribution/" + f + ".up.ef.txt";
        String outputFileDownEf = Dir + "distribution/" + f + ".down.ef.txt";
        String outputFileInterEf = Dir + "distribution/" + f + ".inter.ef.txt";
        GeneFeature gf = new GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[1000];
        double[] efCUp = new double[1000];
        int up = 0;
        int[] countDown = new int[1000];
        double[] efCDown = new double[1000];
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
                if (temp.split("\t")[11].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Math.abs(Double.valueOf(temp.split("\t")[11]));
                }
                int i = gf.getGeneIndex(geneName);
                if (gf.isWithinThisGene(i, Integer.valueOf(temp.split("\t")[2].split("_")[0]), Integer.valueOf(temp.split("\t")[2].split("_")[1]))) {
                    pos = Integer.valueOf(temp.split("\t")[2].split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {
                        start = gf.getGeneStart(i);
                        end = gf.getGeneEnd(i);
                        length = end - start;
                        int chunk = (pos - start) * 100 / length;
                        int number = 1000 / (length / 100);
                        countInter[chunk] += number;
                        Inter++;
                        efCInter[chunk] += ef * number;
                    } else {
                        start = gf.getGeneEnd(i);
                        end = gf.getGeneStart(i);
                        length = start - end;
                        int chunk = (pos - end) * 100 / length;
                        int number = 1000 / (length / 100);
                        countInter[chunk] += number;
                        Inter++;
                        efCInter[chunk] += ef * number;
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
            ex.printStackTrace();
        }
    }

    public void getDistanceEffectSizehumannominal() {
//        String mode = null;
//        mode = "all";
//        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < ; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
        String info = "/Users/yxh/Documents/eQTL/GTEx/Whole_Blood.v8.egenes.txt";
        String infile = "/Users/yxh/Documents/eQTL/GTEx/GTEx_Analysis_v8_eQTL/Whole_Blood.v8.signif_variant_gene_pairs.txt";
//        String f = "all";
//        nameSet.parallelStream().forEach(f -> {
//        String Dir = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/";
//        String homogenefile = "/data1/home/xiaohan/reference/TheABD.txt";
//        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/distribution_up.txt";
        String outputFileDown = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/distribution_down.txt";
        String outputFileInter = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/distribution_inter.txt";
        String outputFileUpEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/ef_up.txt";
        String outputFileDownEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/ef_down.txt";
        String outputFileInterEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance_nominal/ef_inter.txt";
//        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[1000];
        double[] efCUp = new double[1000];
        int up = 0;
        int[] countDown = new int[1000];
        double[] efCDown = new double[1000];
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
            BufferedReader brinfo = IOUtils.getTextReader(info);
            String tempinfo = null;
            String[] tempinfos = null;
            HashMap<String, Integer> geneStrandMap = new HashMap<>();
            HashMap<String, String> geneEndMap = new HashMap<>();
            HashMap<String, String> geneStartMap = new HashMap<>();
            while ((tempinfo = brinfo.readLine()) != null) {
                tempinfos = tempinfo.split("\t");
                geneStartMap.put(tempinfos[0], tempinfos[3]);
                geneEndMap.put(tempinfos[0], tempinfos[4]);
                if (tempinfos[5].equals("+")) {
                    geneStrandMap.put(tempinfos[0], 1);
                } else if (tempinfos[5].equals("-")) {
                    geneStrandMap.put(tempinfos[0], 0);
                }
            }
            brinfo.close();
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedWriter bwUp = IOUtils.getTextWriter(outputFileUp);
            BufferedWriter bwDown = IOUtils.getTextWriter(outputFileDown);
            BufferedWriter bwInter = IOUtils.getTextWriter(outputFileInter);
            BufferedWriter bwUpEf = IOUtils.getTextWriter(outputFileUpEf);
            BufferedWriter bwDownEf = IOUtils.getTextWriter(outputFileDownEf);
            BufferedWriter bwInterEf = IOUtils.getTextWriter(outputFileInterEf);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("variant")) continue;
                String[] temps = temp.split("\t");
                geneName = temps[1];
                int snp = Integer.parseInt(temps[0].split("_")[1]);
//                if (String.valueOf(geneName.charAt(8)).equals("D")) {
//                if(getHomoGene.ishomoGene(homogenefile,geneName).equals("A")){
                if (temps[7].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Math.abs(Double.valueOf(temps[7]));
//                    ef = Double.valueOf(temps[30]);
                }
                if (ef > 2 / Math.log10(2)) continue;
//                if (temps[31].equals("nan") || temps[32].equals("nan")) continue;
//                if (Double.valueOf(temps[31]) * Double.parseDouble(temps[32]) < 0) continue;
//                int i = gf.getGeneIndex(geneName);
                int strand = geneStrandMap.get(geneName);
                if (Integer.parseInt(geneStartMap.get(geneName)) < snp && snp < Integer.parseInt(geneEndMap.get(geneName))) {
                    pos = snp;
                    if (strand == 1) {
                        start = Integer.parseInt(geneStartMap.get(geneName));
                        end = Integer.parseInt(geneEndMap.get(geneName));
                        length = end - start;
                        int chunk = (pos - start - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    } else {
                        start = Integer.parseInt(geneEndMap.get(geneName));
                        end = Integer.parseInt(geneStartMap.get(geneName));
                        length = start - end;
                        int chunk = (start - pos - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    }
                } else if (Integer.parseInt(temps[3]) > snp || snp > Integer.parseInt(temps[4])) {
                    pos = snp;
                    if (strand == 1) {//1表示的是正链
                        start = Integer.parseInt(geneStartMap.get(geneName));
                        if (start >= pos) {
                            int chunk = (start - pos - 1) / 1000;
                            if (chunk > 99) continue;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = Integer.parseInt(geneEndMap.get(geneName));
                            int chunk = (pos - end - 1) / 1000;
                            if (chunk > 99) continue;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    } else {
                        start = Integer.parseInt(geneStartMap.get(geneName));
                        if (start <= pos) {
                            int chunk = (pos - start - 1) / 1000;
                            if (chunk > 99) continue;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = Integer.parseInt(geneEndMap.get(geneName));
                            int chunk = (end - pos - 1) / 1000;
                            if (chunk > 99) continue;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    }
                }
//                }
            }
//            for (int i = 0; i < efCInter.length; i++) {
//                efCInter[i] = 0;
//            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i * 10 + 1001;
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

            for (int i = 0; i < countDown.length; i++) {
                int chunk = i + 2001;
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
            ex.printStackTrace();
            ex.printStackTrace();
        }
//        });
    }

    public void getDistanceEffectSizehuman() {
//        String mode = null;
//        mode = "all";
//        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < ; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
        String infile = "/Users/yxh/Documents/eQTL/GTEx/Whole_Blood.v8.egenes.txt";
//        String f = "all";
//        nameSet.parallelStream().forEach(f -> {
        String Dir = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/";
//        String homogenefile = "/data1/home/xiaohan/reference/TheABD.txt";
//        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/distribution_up.txt";
        String outputFileDown = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/distribution_down.txt";
        String outputFileInter = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/distribution_inter.txt";
        String outputFileUpEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/ef_up.txt";
        String outputFileDownEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/ef_down.txt";
        String outputFileInterEf = "/Users/yxh/Documents/eQTL/GTEx/GTEx_distance/ef_inter.txt";
//        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[100];
        double[] efCUp = new double[100];
        int up = 0;
        int[] countDown = new int[100];
        double[] efCDown = new double[100];
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
            BufferedReader br = IOUtils.getTextReader(infile);
            BufferedWriter bwUp = IOUtils.getTextWriter(outputFileUp);
            BufferedWriter bwDown = IOUtils.getTextWriter(outputFileDown);
            BufferedWriter bwInter = IOUtils.getTextWriter(outputFileInter);
            BufferedWriter bwUpEf = IOUtils.getTextWriter(outputFileUpEf);
            BufferedWriter bwDownEf = IOUtils.getTextWriter(outputFileDownEf);
            BufferedWriter bwInterEf = IOUtils.getTextWriter(outputFileInterEf);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("gene")) continue;
                String[] temps = temp.split("\t");
                geneName = temps[0];
                int snp = Integer.parseInt(temps[14]);
//                if (String.valueOf(geneName.charAt(8)).equals("D")) {
//                if(getHomoGene.ishomoGene(homogenefile,geneName).equals("A")){
                if (temps[30].equals("nan")) {
                    continue;
                } else {
                    ef = Math.abs(Double.valueOf(temps[30]));
//                    ef = Double.valueOf(temps[30]);
                }
                if (ef > 2 / Math.log10(2)) continue;
                if (temps[31].equals("nan") || temps[32].equals("nan")) continue;
                if (Double.valueOf(temps[31]) * Double.parseDouble(temps[32]) < 0) continue;
//                int i = gf.getGeneIndex(geneName);
                int strand = 0;
                if (temps[5].equals("+")) {
                    strand = 1;
                }
                if (Integer.parseInt(temps[3]) < snp && snp < Integer.parseInt(temps[4])) {
                    pos = snp;
                    if (strand == 1) {
                        start = Integer.parseInt(temps[3]);
                        end = Integer.parseInt(temps[4]);
                        length = end - start;
                        int chunk = (pos - start - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    } else {
                        start = Integer.parseInt(temps[4]);
                        end = Integer.parseInt(temps[3]);
                        length = start - end;
                        int chunk = (start - pos - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    }
                } else if (Integer.parseInt(temps[3]) > snp || snp > Integer.parseInt(temps[4])) {
                    pos = snp;
                    if (strand == 1) {//1表示的是正链
                        start = Integer.parseInt(temps[3]);
                        if (start >= pos) {
                            int chunk = (start - pos - 1) / 1000;
                            if (chunk <= 99) {
                                countUp[chunk]++;
                                up++;
                                efCUp[chunk] += ef;
                            }
                        } else {
                            end = Integer.parseInt(temps[4]);
                            int chunk = (pos - end - 1) / 1000;
                            if (chunk <= 99) {
                                countDown[chunk]++;
                                down++;
                                efCDown[chunk] += ef;
                            }
                        }
                    } else {
                        start = Integer.parseInt(temps[4]);
                        if (start <= pos) {
                            int chunk = (pos - start - 1) / 1000;
                            if (chunk <= 99) {
                                countUp[chunk]++;
                                up++;
                                efCUp[chunk] += ef;
                            }
                        } else {
                            end = Integer.parseInt(temps[3]);
                            int chunk = (end - pos - 1) / 1000;
                            if (chunk <= 99) {
                                countDown[chunk]++;
                                down++;
                                efCDown[chunk] += ef;
                            }
                        }
                    }
                }
//                }
            }
//            for (int i = 0; i < efCInter.length; i++) {
//                efCInter[i] = 0;
//            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i + 101;
                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
                if (efCInter[i] == 0) {
                    bwInterEf.write("Inter" + "\t" + chunk + "\t" + "0.0" + "\n");
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
                int chunk = 100 - i;
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
                int chunk = i + 201;
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
            ex.printStackTrace();
            ex.printStackTrace();
        }
//        });
    }

    public void getDistanceEffectSizeNominal() {
        String mode = null;
//        mode = "all";
        HashSet<String> nameSet = new HashSet<>();
//        for (int m = 0; m < ; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
        String f = "all";
//        nameSet.parallelStream().forEach(f -> {
        String Dir = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect/";
        String homogenefile = "/data2/xiaohan/reference/TheABD.txt";
        String inputFileS = Dir + f + ".cis.sig.DE.log2.txt";
        String outputFileUp = Dir + "distribution/" + f + ".up.distribution.txt";
        String outputFileDown = Dir + "distribution/" + f + ".down.distribution.txt";
        String outputFileInter = Dir + "distribution/" + f + ".inter.distribution.txt";
        String outputFileUpEf = Dir + "distribution/" + f + ".up.ef.txt";
        String outputFileDownEf = Dir + "distribution/" + f + ".down.ef.txt";
        String outputFileInterEf = Dir + "distribution/" + f + ".inter.ef.txt";
        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        int[] countUp = new int[1000];
        double[] efCUp = new double[1000];
        int up = 0;
        int[] countDown = new int[1000];
        double[] efCDown = new double[1000];
        int down = 0;
        int[] countInter = new int[100];
        double[] efCInter = new double[100];
        for (int i = 0; i < 1000; i++) {
            countDown[i] = 0;
            countUp[i] = 0;
            efCDown[i] = 0;
            efCUp[i] = 0;
        }
        for (int i = 0; i < 100; i++) {
            countInter[i] = 0;
            efCInter[i] = 0;
        }
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
                String[] temps = temp.split("\t");
                geneName = temps[1];
                String snp = temps[2];
//                if (String.valueOf(geneName.charAt(8)).equals("D")) {
//                if(getHomoGene.ishomoGene(homogenefile,geneName).equals("A")){
                if (temps[11].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Math.abs(Double.valueOf(temps[11]));
                }
                if (ef > 2 / Math.log10(2)) continue;
                if (temps[12].equals("nan") || temps[13].equals("nan")) continue;
                if (Double.valueOf(temps[12]) * Double.parseDouble(temps[13]) < 0) continue;
                int i = gf.getGeneIndex(geneName);
                if (gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {
                        start = gf.getGeneStart(i);
                        end = gf.getGeneEnd(i);
                        length = end - start;
                        int chunk = (pos - start - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    } else {
                        start = gf.getGeneEnd(i);
                        end = gf.getGeneStart(i);
                        length = start - end;
                        int chunk = (start - pos - 1) * 100 / length;
//                            int number = 1000 / (length / 100);
                        countInter[chunk]++;
                        Inter++;
                        efCInter[chunk] += ef;
                    }
                } else if (!gf.isWithinThisGene(i, Integer.valueOf(snp.split("_")[0]), Integer.valueOf(snp.split("_")[1]))) {
                    pos = Integer.valueOf(snp.split("_")[1]);
                    if (gf.getGeneStrand(i) == 1) {//1表示的是正链
                        start = gf.getGeneStart(i);
                        if (start >= pos) {
                            int chunk = (start - pos - 1) / 1000;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = gf.getGeneEnd(i);
                            int chunk = (pos - end - 1) / 1000;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    } else {
                        start = gf.getGeneEnd(i);
                        if (start <= pos) {
                            int chunk = (pos - start - 1) / 1000;
                            countUp[chunk]++;
                            up++;
                            efCUp[chunk] += ef;
                        } else {
                            end = gf.getGeneStart(i);
                            int chunk = (end - pos - 1) / 1000;
                            countDown[chunk]++;
                            down++;
                            efCDown[chunk] += ef;
                        }
                    }
                }
//                }
            }
//            for (int i = 0; i < efCInter.length; i++) {
//                efCInter[i] = 0;
//            }
            br.close();
            DecimalFormat decFor = new DecimalFormat("0.000000");
            for (int i = 0; i < countInter.length; i++) {
                int chunk = i * 10 + 1001;
                bwInter.write("Inter" + "\t" + chunk + "\t" + countInter[i] + "\n");
                if (efCInter[i] == 0) {
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

            for (int i = 0; i < countDown.length; i++) {
                int chunk = i + 2001;
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
            ex.printStackTrace();
        }
//        });
    }

    public void getMafWithEffectSize() {
        String f = "all";
//        HashSet<String> nameSet = new HashSet();
//        for (int m = 0; m < 42; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
//        nameSet.parallelStream().forEach(f -> {
        String infile = "/data1/home/xiaohan/tensorQTL/1M_log2/homoeffect/" + f + ".cis.sig.DE.log2.txt";
        String outfile = "/data1/home/xiaohan/tensorQTL/1M_log2/homoeffect/selection/" + f + ".selection.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temps[11].startsWith("log")) {
                    bw.write("Index\tGene\tID\tEff\tMaf\tGroup");
                    bw.newLine();
                    continue;
                }
                double ef = 0;
                if (temps[11].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Double.parseDouble(temps[11]);
                }
                String group = null;
                if (Math.abs(ef) <= 1) {
                    group = "Small";
                } else if (Math.abs(ef) > 1 && Math.abs(ef) <= 2) {
                    group = "Medium";
                } else if (Math.abs(ef) > 2) {
                    group = "Large";
                }
                StringBuilder sb = new StringBuilder();
                sb.append(temps[0] + "\t" + temps[1] + "\t" + temps[2] + "\t" + ef + "\t" + temps[4] + "\t" + group);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
//        });
    }

    public void getGerpEffectSizeTop() {
        String f = "all";
        String infile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/" + f + ".cis.sig.DE.log2.txt";
        String outfile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/Gerp/" + f + ".gerp.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        HashSet<String> annoSet = new HashSet<>();
        HashMap<String, String> annoMap = new HashMap<>();
        BufferedReader brinfo = IOUtils.getTextGzipReader("/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/Gerpsnp/all.sorted.bed.gz");
        String info = null;
        String[] infos = null;
        String temp = null;
        String[] temps = null;
        try {
            System.out.println("Reading ....................");
            int countline = 0;
            while ((info = brinfo.readLine()) != null) {
                countline++;
                if (countline % 500 == 0) System.out.println(countline);
                infos = info.split("\t");
                if (info.startsWith("ID")) continue;
                String name = infos[0] + "_" + infos[2];
                annoSet.add(name);
                annoMap.put(name, infos[4]);
//                System.out.println(name + "\t" + infos[12] + "\n");
            }
            brinfo.close();
            countline = 0;
            System.out.println("Writing ....................");
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                countline++;
                if (countline % 500 == 0) System.out.println(countline);
                if (temps[0].equals("pid")) {
                    bw.write(temps[0] + "\t" + temps[1] + "\t" + "Eff" + "\t" + "Maf" + "\t" + "Group" + "\t" + "Gerp");
                    bw.newLine();
                    continue;
                }
                String snp = temps[1];
                double ef = 0;
                if (temps[6].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Double.parseDouble(temps[6]);
                }
                if (ef >= (2 / Math.log10(2))) continue;
                if (temps[7].equals("nan") || temps[8].equals("nan")) continue;
                if (Double.valueOf(temps[7]) * Double.parseDouble(temps[8]) < 0) continue;
                String group = null;
                if (annoMap.get(snp) == null) {
                    group = "NA";
                    StringBuilder sb = new StringBuilder();
                    sb.append(temps[0] + "\t" + temps[1] + "\t" + ef + "\t" + temps[4] + "\t" + group + "\t" + "NA");
                    bw.write(sb.toString());
                    bw.newLine();
                } else if (Double.parseDouble(annoMap.get(snp)) <= 1) {
                    group = "Gerp0-1";
                    StringBuilder sb = new StringBuilder();
                    sb.append(temps[0] + "\t" + temps[1] + "\t" + ef + "\t" + temps[4] + "\t" + group + "\t" + Double.parseDouble(annoMap.get(snp)));
                    bw.write(sb.toString());
                    bw.newLine();
                } else if (Double.parseDouble(annoMap.get(snp)) <= 2 && Double.parseDouble(annoMap.get(snp)) > 1) {
                    group = "Gerp1-2";
                    StringBuilder sb = new StringBuilder();
                    sb.append(temps[0] + "\t" + temps[1] + "\t" + ef + "\t" + temps[4] + "\t" + group + "\t" + Double.parseDouble(annoMap.get(snp)));
                    bw.write(sb.toString());
                    bw.newLine();
                } else if (Double.parseDouble(annoMap.get(snp)) > 2) {
                    group = "Gerp>2";
                    StringBuilder sb = new StringBuilder();
                    sb.append(temps[0] + "\t" + temps[1] + "\t" + ef + "\t" + temps[4] + "\t" + group + "\t" + Double.parseDouble(annoMap.get(snp)));
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


    public void getRegionEffectSizeNominal() {
        String annotationfile = "/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3";
        GeneFeature gf = new GeneFeature(annotationfile);
//        HashSet<String> nameSet = new HashSet();
//        for (int m = 0; m < 42; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
//        nameSet.parallelStream().forEach(f -> {
        String f = "all";
        String infile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect/" + f + ".cis.sig.DE.log2.txt";
        String outfile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect/region/" + f + ".region.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        HashSet<String> annoSet = new HashSet<>();
        HashMap<String, String> annoMap = new HashMap<>();
        BufferedReader brinfo = IOUtils.getTextReader("/data2/xiaohan/tensorQTL/002_exonSNP_anno.txt");
        String info = null;
        String[] infos = null;
        String temp = null;
        String[] temps = null;
        try {
            while ((info = brinfo.readLine()) != null) {
                infos = info.split("\t");
                if (info.startsWith("ID")) continue;
                String name = infos[0] + "_" + infos[1];
                annoSet.add(name);
                annoMap.put(name, infos[6]);
//                System.out.println(name + "\t" + infos[12] + "\n");
            }
            brinfo.close();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temps[0].equals("Index")) {
                    bw.write(temps[1] + "\t" + temps[2] + "\t" + "Eff" + "\t" + "Maf" + "\t" + "Group");
                    bw.newLine();
                    continue;
                }
                String ID = temps[2];
                int pos = Integer.parseInt(temps[2].split("_")[1]);
                int chr = Integer.parseInt(temps[2].split("_")[0]);
                String geneName = temps[1];
                int index = gf.getGeneIndex(geneName);
                int startsite = gf.getGeneStart(index);
                int endsite = gf.getGeneEnd(index);
                int strand = gf.getGeneStrand(index);
                double ef = 0;
                if (temps[11].equals("nan")) {
                    ef = 0;
                } else {
                    ef = Double.parseDouble(temps[11]);
                }
                if (ef >= (2 / Math.log10(2))) continue;
                if (temps[12].equals("nan") || temps[13].equals("nan")) continue;
                if (Double.valueOf(temps[12]) * Double.parseDouble(temps[13]) < 0) continue;
                String group = null;
//                String position = "not";
                int j = gf.getLongestTranscriptIndex(index);
                if (pos < startsite || pos > endsite) {
                    if (strand == 1) {
                        if (pos <= startsite) {
                            group = "Upstream";
                        } else {
                            group = "Downstream";
                        }
                    } else {
                        if (pos >= endsite) {
                            group = "Upstream";
                        } else {
                            group = "Downstream";
                        }
                    }
                } else if (pos >= startsite && pos <= endsite) {
                    group = "Intron";
                    if (gf.get5UTRList(index, j).size() != 0) {
                        List<Range> fU = gf.get5UTRList(index, j);
                        for (int a = 0; a < fU.size(); a++) {
                            int start = fU.get(a).getRangeStart();
                            int end = fU.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "FiveUTR";
                            }
                            continue;
                        }
                    }
                    if (gf.getCDSList(index, j).size() != 0) {
                        List<Range> CDS = gf.getCDSList(index, j);
                        for (int a = 0; a < CDS.size(); a++) {
                            int start = CDS.get(a).getRangeStart();
                            int end = CDS.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "CDS";
                                if (annoSet.contains(ID)) {
                                    group = annoMap.get(ID);
                                }
                                continue;
                            }
                        }
                    }
                    if (gf.get3UTRList(index, j).size() != 0) {
//                        System.out.println("3'");
                        List<Range> tU = gf.get3UTRList(index, j);
                        for (int a = 0; a < tU.size(); a++) {
                            int start = tU.get(a).getRangeStart();
                            int end = tU.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "ThreeUTR";
                            }
                            continue;
                        }
                    }
                }
//                if (group == null) {
//                    System.out.println(chr + "\t" + strand + "\t" + startsite + "\t" + endsite + "\t" + pos + "\t" + position);
//                }
                StringBuilder sb = new StringBuilder();
                sb.append(temps[1] + "\t" + temps[2] + "\t" + ef + "\t" + temps[4] + "\t" + group);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
//        });
    }

    public void getRegionEffectSizeTop() {
        String annotationfile = "/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3";
        GeneFeature gf = new GeneFeature(annotationfile);
//        HashSet<String> nameSet = new HashSet();
//        for (int m = 0; m < 42; m++) {
//            int chr = m + 1;
//            nameSet.add(String.valueOf(chr));
//        }
//        nameSet.parallelStream().forEach(f -> {
        String f = "all";
        String infile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/" + f + ".cis.sig.DE.log2.txt";
        String outfile = "/data2/xiaohan/tensorQTL/1M_log2/homoeffect_top/region/" + f + ".region.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        HashSet<String> annoSet = new HashSet<>();
        HashMap<String, String> annoMap = new HashMap<>();
        BufferedReader brinfo = IOUtils.getTextReader("/data2/xiaohan/tensorQTL/002_exonSNP_anno.txt");
        String info = null;
        String[] infos = null;
        String temp = null;
        String[] temps = null;
        try {
            while ((info = brinfo.readLine()) != null) {
                infos = info.split("\t");
                if (info.startsWith("ID")) continue;
                String name = infos[0] + "_" + infos[1];
                annoSet.add(name);
                annoMap.put(name, infos[6]);
//                System.out.println(name + "\t" + infos[12] + "\n");
            }
            brinfo.close();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temps[0].equals("pid")) {
//                    + "\t" + temps[2]
                    bw.write(temps[0] + "\t" + temps[1] + "\t" + "Eff" + "\t" + "Maf" + "\t" + "Group");
                    bw.newLine();
                    continue;
                }
                int pos = Integer.parseInt(temps[1].split("_")[1]);
                int chr = Integer.parseInt(temps[1].split("_")[0]);
                String geneName = temps[0];
                String ID = temps[1];
                int index = gf.getGeneIndex(geneName);
                int startsite = gf.getGeneStart(index);
                int endsite = gf.getGeneEnd(index);
                int strand = gf.getGeneStrand(index);
                double ef = 0;
                if (temps[6].equals("nan")) {
//                    ef = 0;
                    continue;
                } else {
                    ef = Double.parseDouble(temps[6]);
                }
                if (ef >= (2 / Math.log10(2))) continue;
                if (temps[7].equals("nan") || temps[8].equals("nan")) continue;
                if (Double.valueOf(temps[7]) * Double.parseDouble(temps[8]) < 0) continue;
                String group = null;
//                String position = "not";
                int j = gf.getLongestTranscriptIndex(index);
                if (pos < startsite || pos > endsite) {
                    if (strand == 1) {
                        if (pos <= startsite) {
                            group = "Upstream";
                        } else {
                            group = "Downstream";
                        }
                    } else {
                        if (pos >= endsite) {
                            group = "Upstream";
                        } else {
                            group = "Downstream";
                        }
                    }
                } else if (pos >= startsite && pos <= endsite) {
                    group = "Intron";
                    if (gf.get5UTRList(index, j).size() != 0) {
                        List<Range> fU = gf.get5UTRList(index, j);
                        for (int a = 0; a < fU.size(); a++) {
                            int start = fU.get(a).getRangeStart();
                            int end = fU.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "FiveUTR";
                            }
                            continue;
                        }
                    }
                    if (gf.getCDSList(index, j).size() != 0) {
                        List<Range> CDS = gf.getCDSList(index, j);
                        for (int a = 0; a < CDS.size(); a++) {
                            int start = CDS.get(a).getRangeStart();
                            int end = CDS.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "CDS";
                                if (annoSet.contains(ID)) {
                                    group = annoMap.get(ID);
                                }
                                continue;
                            }
                        }
                    }
                    if (gf.get3UTRList(index, j).size() != 0) {
//                        System.out.println("3'");
                        List<Range> tU = gf.get3UTRList(index, j);
                        for (int a = 0; a < tU.size(); a++) {
                            int start = tU.get(a).getRangeStart();
                            int end = tU.get(a).getRangeEnd();
                            if (pos >= start && pos <= end) {
                                group = "ThreeUTR";
                            }
                            continue;
                        }
                    }
                }
//                if (group == null) {
//                    System.out.println(chr + "\t" + strand + "\t" + startsite + "\t" + endsite + "\t" + pos + "\t" + position);
//                }
                StringBuilder sb = new StringBuilder();
//                + "\t" + temps[2]
                sb.append(temps[0] + "\t" + temps[1] + "\t" + ef + "\t" + temps[4] + "\t" + group);
                bw.write(sb.toString());
                bw.newLine();
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
//        });
    }



    public static void main (String[] args){
        new effectsize(args);
    }
}
