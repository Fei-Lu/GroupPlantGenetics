package xiaohan.eQTL.simulation;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

public class simulationData {
    public simulationData(String[] args) {
        this.generateFastq(args);
        this.generateReads(args[0]);
        this.getTrueSet();
        this.getReadslength(args[0]);
        this.getalignment(args[0]);
        this.getCount(args[0]);
        this.getRealCount(args[0]);
        this.getSummary(args[0]);
        this.getShuf(args[0]);
        this.getPrecisionRecall(args[0]);
        this.getPrecisionPlot(args[0]);
        this.getPrecisionandRecall(args[0],args[1]);
    }

    public static double[] getPrecisionandRecall(String tru, String obs) {
        double TP = 0;
        double FN = 0;
        double FP = 0;
        double observedvalue = Double.parseDouble(obs);
        double truevalue = Double.parseDouble(tru);
        if (observedvalue >= truevalue) {
            TP += truevalue;
            FP += observedvalue - truevalue;
        }
        if (observedvalue <= truevalue) {
            TP += observedvalue;
            FN += truevalue - observedvalue;
        }
        double[] TPFPFN = {TP,FP,FN};
        return TPFPFN;
    }

    private String getphred(int score) {
        int[] scores = new int[41];
        for (int i = 0; i < 41; i++) {
            scores[i] = i;
        }
        String[] phreds = {"!", "\"", "#", "$", "%", "&", " ", "(", ")", "*", "+", ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I"};
        HashMap<Integer, String> scoreQMap = new HashMap<>();
        for (int i = 0; i < 41; i++) {
            scoreQMap.put(scores[i], phreds[i]);
        }
        String phred = scoreQMap.get(score);
        return phred;
    }

    public void getPrecisionPlot(String infileDir) {
        String[] dir = {"SE", "PE"};
        for (int m = 0; m < dir.length; m++) {
            String infile = new File(infileDir, dir[m] + "boot/simulation.txt").getAbsolutePath();
            String outfile = new File(infileDir, dir[m] + "boot/simulationforplot.txt").getAbsolutePath();
            BufferedReader br = IOUtils.getTextReader(infile);
            String temp = null;
            String[] temps = null;
            try {
                int countline = 0;
                double[][] precision = new double[101][100];
                double[][] recall = new double[101][100];
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    if (temp.startsWith("FileName")) continue;
                    precision[countline / 100][countline % 100] = Double.parseDouble(temps[4]);
                    recall[countline / 100][countline % 100] = Double.parseDouble(temps[5]);
                    countline++;
                }
                br.close();
                double[] precisionsum = new double[101];
                double[] precisionave = new double[101];
                double[] precisionsd = new double[101];
                double[] recallsum = new double[101];
                double[] recallave = new double[101];
                double[] recallsd = new double[101];
                for (int i = 0; i < 101; i++) {
                    for (int j = 0; j < 100; j++) {
                        precisionsum[i] += precision[i][j];
                        recallsum[i] += recall[i][j];
                    }
                    precisionave[i] = (double) precisionsum[i] / 100;
                    recallave[i] = (double) recallsum[i] / 100;
                    for (int j = 0; j < 100; j++) {
                        precisionsd[i] += (precision[i][j] - precisionave[i]) * (precision[i][j] - precisionave[i]);
                        recallsd[i] += (recall[i][j] - recallave[i]) * (recall[i][j] - recallave[i]);
                    }
                    precisionsd[i] = Math.sqrt(precisionsd[i] / 99);
                    recallsd[i] = Math.sqrt(recallsd[i] / 99);
                }
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                DecimalFormat decfor = new DecimalFormat("0.00000000");
                bw.write("size\tprecision\tsd1\trecall\tsd2\n");
                for (int i = 0; i < 101; i++) {
                    int index = 50 + i;
                    bw.write(index + "\t" + decfor.format(precisionave[i]) + "\t" + decfor.format(precisionsd[i]) + "\t" + decfor.format(recallave[i]) + "\t" + decfor.format(recallsd[i]) + "\n");
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getSummary(String arg) {
        for (int i = 50; i < 151; i++) {
            String command = null;
            StringBuilder sb4 = new StringBuilder();
            sb4.append("cp SE/SE3724-" + i + "_Count.txt SE/SE3724-" + i + "_forcalu.txt\n");
            sb4.append("sed -i -e '/__/d' SE/SE3724-" + i + "_forcalu.txt\n");
            sb4.append("cut -f 2 SE/SE3724-" + i + "_forcalu.txt > SE/SE3724-" + i + "_forcalu_2.txt\n");
            sb4.append("paste TrueSet.txt SE/SE3724-" + i + "_forcalu_2.txt > SE/SE3724-" + i + ".simu\n");
            sb4.append("cat SE/SE3724-" + i + ".simu | awk -F '\\t' '$2!=0||$3!=0{print $0}' > SE/SE3724-" + i + ".s\n");
            command = sb4.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished summary");
        for (int i = 50; i < 151; i++) {
            String command = null;
            StringBuilder sb4 = new StringBuilder();
            sb4.append("cp PE/SE3724-" + i + "_Count.txt PE/SE3724-" + i + "_forcalu.txt\n");
            sb4.append("sed -i -e '/__/d' PE/SE3724-" + i + "_forcalu.txt\n");
            sb4.append("cut -f 2 PE/SE3724-" + i + "_forcalu.txt > PE/SE3724-" + i + "_forcalu_2.txt\n");
            sb4.append("paste TrueSet.txt PE/SE3724-" + i + "_forcalu_2.txt > PE/SE3724-" + i + ".simu\n");
            sb4.append("cat PE/SE3724-" + i + ".simu | awk -F '\\t' '$2!=0||$3!=0{print $0}' > PE/SE3724-" + i + ".s\n");
            command = sb4.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished summary");
    }

    public void getPrecisionRecall(String infileDir) {
        String[] dir = {"SE", "PE"};
        for (int m = 0; m < dir.length; m++) {
            String outfile = new File(infileDir, dir[m] + "boot/simulation.txt").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            try {
                bw.write("FileName\tTP\tFP\tFN\tPrecision\tRecall\n");
                for (int i = 50; i < 151; i++) {
                    for (int j = 1; j < 101; j++) {
                        int number = i;
                        String index = String.valueOf(number);
                        BufferedReader br = IOUtils.getTextReader(new File(infileDir, dir[m] + "boot/SE3724-" + i + "-" + j).getAbsolutePath());
//                    BufferedReader br = IOUtils.getTextReader("/data1/home/xiaohan/SiPASsimu/SE3724-" + i +".s");
                        String temp = null;
                        String[] temps = null;
                        int TP = 0;
                        int FN = 0;
                        int FP = 0;
                        double precision = 0.00000000;
                        double recall = 0.00000000;
                        try {
                            while ((temp = br.readLine()) != null) {
                                temps = temp.split("\t");
                                int truevalue = Integer.parseInt(temps[1]);
                                int observedvalue = Integer.parseInt(temps[2]);
                                if (truevalue > 0 && observedvalue >= truevalue) {
                                    TP += truevalue;
                                    if (observedvalue > truevalue) {
                                        FP += observedvalue - truevalue;
                                    }
                                }
                                if (observedvalue > 0 && observedvalue <= truevalue) {
                                    TP += observedvalue;
                                    if (observedvalue < truevalue) {
                                        FN += truevalue - observedvalue;
                                    }
                                }
                                if (truevalue > 0 && observedvalue == 0) {
                                    FN += truevalue;
                                }
                                if (truevalue == 0 && observedvalue > 0) {
                                    FP += observedvalue;
                                }
                            }
                            br.close();
                            DecimalFormat decfor = new DecimalFormat("0.00000000");
                            precision = (double) TP / (TP + FP);
                            recall = (double) TP / (TP + FN);
                            StringBuffer sb = new StringBuffer();
//                        sb.append(i +"\t" + TP + "\t" + FP + "\t" + FN + "\t" + decfor.format(precision) + "\t" + decfor.format(recall) + "\n");
                            sb.append(i + "-" + j + "\t" + TP + "\t" + FP + "\t" + FN + "\t" + decfor.format(precision) + "\t" + decfor.format(recall) + "\n");
                            bw.write(sb.toString());
                        } catch (Exception e) {
                            e.printStackTrace();
                        }
                    }
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getShuf(String arg) {
        String command = null;
        StringBuilder sb5 = new StringBuilder();
        sb5.append("mkdir SEboot\n" +
                "for m in $(seq 1 100)\n" +
                "do\n" +
                "\tfor c in $(seq 50 150)\n" +
                "\tdo shuf -n 5000 SE/SE3724-${c}.s > SEboot/SE3724-${c}-${m}\n" +
                "done\n" +
                "done");
        command = sb5.toString();
        System.out.println(command);
        try {
            File dir = new File(new File(arg).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished shufling");
        command = null;
        StringBuilder sb6 = new StringBuilder();
        sb6.append("mkdir PEboot\n" +
                "for m in $(seq 1 100)\n" +
                "do\n" +
                "\tfor c in $(seq 50 150)\n" +
                "\tdo shuf -n 5000 PE/SE3724-${c}.s > PEboot/SE3724-${c}-${m}\n" +
                "done\n" +
                "done");
        command = sb6.toString();
        System.out.println(command);
        try {
            File dir = new File(new File(arg).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("Finished shufling");
    }

    public void getRealCount(String arg) {
        String infile = new File(arg, "SE3724-150_R2.fq.gz").getAbsolutePath();
        String infile1 = new File(arg, "SE/SE3724-150_Count.txt").getAbsolutePath();
        String outfile = new File(arg, "TrueSet.txt").getAbsolutePath();
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        try {
            HashSet<String> geneSet = new HashSet<>();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("@")) {
                    String geneName = temp.split("@")[1];
                    geneSet.add(geneName);
                }
            }
            br.close();
            String[] genelist = geneSet.toArray(new String[geneSet.size()]);
            int[] number = new int[genelist.length];
            for (int i = 0; i < genelist.length; i++) {
                number[i] = 0;
            }
            Arrays.sort(genelist);
            HashMap<String, Integer> geneIndexmap = new HashMap<>();
            for (int i = 0; i < genelist.length; i++) {
                geneIndexmap.put(genelist[i], i);
            }
            BufferedReader br1 = IOUtils.getTextGzipReader(infile);
            while ((temp = br1.readLine()) != null) {
                if (temp.startsWith("@")) {
                    String geneName = temp.split("@")[1];
//                    System.out.println(geneName);
                    number[geneIndexmap.get(geneName)]++;
                }
            }
            br1.close();
            HashMap<String, Integer> geneNumbermap = new HashMap<>();
            for (int i = 0; i < genelist.length; i++) {
                geneNumbermap.put(genelist[i], number[i]);
            }
            BufferedReader br2 = IOUtils.getTextReader(infile1);
            while ((temp = br2.readLine()) != null) {
                String geneName = temp.split("\t")[0];
                if (geneName.startsWith("T")) {
                    if (!geneSet.contains(geneName)) {
                        bw.write(geneName + "\t0\n");
                    } else {
                        bw.write(geneName + "\t" + geneNumbermap.get(geneName) + "\n");
                    }
                }
            }
            br2.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getCount(String arg) {
        for (int i = 50; i < 151; i++) {
            String command = null;
            StringBuilder sb3 = new StringBuilder();
            sb3.append("htseq-count -f bam -m intersection-nonempty -s no SE/SE3724-" + i + "_Aligned.out.bam /data1/home/junxu/wheat_v1.1_Lulab.gtf >  SE/SE3724-" + i + "_Count.txt\n");
            sb3.append("htseq-count -f bam -m intersection-nonempty -s no PE/SE3724-" + i + "_Aligned.out.bam /data1/home/junxu/wheat_v1.1_Lulab.gtf >  PE/SE3724-" + i + "_Count.txt\n");
            command = sb3.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished counting");
    }

    public void getalignment(String arg) {
        String command = null;
        StringBuilder sb5 = new StringBuilder();
        sb5.append("mkdir PE && mkdir SE");
        command = sb5.toString();
        System.out.println(command);
        try {
            File dir = new File(new File(arg).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 50; i < 151; i++) {
            command = null;
            StringBuilder sb2 = new StringBuilder();
            sb2.append("/data1/home/junxu/software/STAR-2.6.1c/bin/Linux_x86_64/STAR --runThreadN 32 --genomeDir /data1/home/junxu/starLib1.1 --genomeLoad LoadAndKeep --readFilesCommand zcat --readFilesIn SE3724-" + i + "_R2.fq.gz --outFileNamePrefix SE/SE3724-" + i + "_ --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM Unsorted --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \n");
            sb2.append("/data1/home/junxu/software/STAR-2.6.1c/bin/Linux_x86_64/STAR --runThreadN 32 --genomeDir /data1/home/junxu/starLib1.1 --genomeLoad LoadAndKeep --readFilesCommand zcat --readFilesIn SE3724-" + i + "_R1.fq.gz SE3724-150_R2.fq.gz --outFileNamePrefix PE/SE3724-" + i + "_ --outFilterMultimapNmax 10 --outFilterMismatchNoverLmax 0.1 --outFilterIntronMotifs RemoveNoncanonicalUnannotated --outSAMtype BAM Unsorted --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \n");
            command = sb2.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished alignment");
    }

    public void getReadslength(String arg) {
        for (int i = 50; i < 150; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append("zcat SE3724-150_R2.fq.gz | cut -c1-" + i + " > SE3724-" + i + "_R2.fq && bgzip SE3724-" + i + "_R2.fq");
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished readslength");
        for (int i = 50; i < 150; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append("zcat SE3724-150_R1.fq.gz | cut -c1-" + i + " > SE3724-" + i + "_R1.fq && bgzip SE3724-" + i + "_R1.fq");
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        System.out.println("Finished readslength");
    }

    public void getTrueSet() {
        String infile = "/Users/yxh/Documents/eQTL/SiPAS/simu/0219/gene.txt";
        String infile1 = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/SE3724-150_Count.txt";
        String outfile = "/Users/yxh/Documents/eQTL/SiPAS/simu/0219/JunTrue.txt";
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        try {
            HashSet<String> geneSet = new HashSet<>();
            HashMap<String, Integer> geneNumbermap = new HashMap<>();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                String geneName = temps[1];
                geneSet.add(geneName);
                geneNumbermap.put(temps[1], Integer.parseInt(temps[0]));
            }
            br.close();
            String[] genelist = geneSet.toArray(new String[geneSet.size()]);
            int[] number = new int[genelist.length];
            BufferedReader br2 = IOUtils.getTextReader(infile1);
            while ((temp = br2.readLine()) != null) {
                String geneName = temp.split("\t")[0];
                if (geneNumbermap.get(geneName) == null) {
                    bw.write(geneName + "\t0\n");
                } else {
                    bw.write(geneName + "\t" + geneNumbermap.get(geneName) + "\n");
                }
            }
            br2.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void generateReads(String arg) {
        String infile = new File(arg, "exons.fq").getAbsolutePath();
        String outfile1 = new File(arg, "SE3724-150_R1.fq").getAbsolutePath();
        String outfile2 = new File(arg, "SE3724-150_R2.fq").getAbsolutePath();
        String outfile = new File(arg, "SE3724-350.fq").getAbsolutePath();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        BufferedWriter bw1 = IOUtils.getTextWriter(outfile1);
        BufferedWriter bw2 = IOUtils.getTextWriter(outfile2);
        HashMap<String, String> geneReadsMap = new HashMap();
        HashMap<Integer, String> IndexGeneMap = new HashMap<>();
        try {
            String temp = null;
            String[] temps = null;
            int countline = 0;
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                String geneName = temps[0].split("\\.")[0];
                String Reads = temps[1];
                geneReadsMap.put(geneName, Reads);
                IndexGeneMap.put(countline, geneName);
                countline++;
            }
            br.close();
            for (int i = 0; i < 100000; i++) {
                Random r = new Random();
                int num = r.nextInt(countline);
                System.out.println("This is reading gene : " + num);
                String geneName = IndexGeneMap.get(num);
                String Reads1 = null;
                String Reads2 = null;
                //reads1 加随机突变率1%
                Reads1 = geneReadsMap.get(geneName).substring(0, 150);
                StringBuilder sbreads1 = new StringBuilder();
                for (int j = 0; j < Reads1.length(); j++) {
                    String index = Reads1.substring(j, j + 1);
                    String index1 = null;
                    num = r.nextInt(100);
                    if (num == 1) {
                        int num1 = r.nextInt(100);
                        if (num1 < 25) {
                            index1 = "A";
                        } else if (num1 >= 25 && num1 < 50) {
                            index1 = "T";
                        } else if (num1 >= 50 && num1 < 75) {
                            index1 = "C";
                        } else {
                            index1 = "G";
                        }
                        sbreads1.append(index1);
                    } else {
                        sbreads1.append(index);
                    }
                }
                Reads1 = sbreads1.toString();
                //reads2 加随机突变率 1% 反向互补
                Reads2 = geneReadsMap.get(geneName).substring(200, 350);
                StringBuilder sbreads2 = new StringBuilder();
                for (int j = 0; j < Reads2.length(); j++) {
                    String index = Reads2.substring(j, j + 1);
                    String index1 = null;
                    num = r.nextInt(100);
                    if (num == 1) {
                        int num1 = r.nextInt(100);
                        if (num1 < 25) {
                            index1 = "A";
                        } else if (num1 >= 25 && num1 < 50) {
                            index1 = "T";
                        } else if (num1 >= 50 && num1 < 75) {
                            index1 = "C";
                        } else {
                            index1 = "G";
                        }
                        sbreads2.append(index1);
                    } else {
                        sbreads2.append(index);
                    }
                }
                Reads2 = sbreads2.toString();
                Reads2 = new StringBuffer(Reads2).reverse().toString();
                StringBuffer sb = new StringBuffer();
                for (int j = 0; j < Reads2.length(); j++) {
                    String index = Reads2.substring(j, j + 1);
                    String index1 = null;
                    if (index.equals("A")) {
                        index1 = "T";
                    }
                    if (index.equals("G")) {
                        index1 = "C";
                    }
                    if (index.equals("T")) {
                        index1 = "A";
                    }
                    if (index.equals("C")) {
                        index1 = "G";
                    }
                    sb.append(index1);
                }
                Reads2 = sb.toString();
                String Reads350 = geneReadsMap.get(geneName);
                bw.write("@" + geneName + "\n");
                bw.write(Reads350 + "\n");

                bw1.write("@" + geneName + "\n");
                bw1.write(Reads1 + "\n");
                bw1.write("+\n");
                StringBuffer sb1 = new StringBuffer();
                for (int j = 0; j < 150; j++) {
                    int score = (j - 150) * (j - 150) * 14 / 22201 + 24;
//                    System.out.println(score);
                    String Q = getphred(score);
                    sb1.append(Q);
                }
                bw1.write(sb1.toString() + "\n");

                bw2.write("@" + geneName + "\n");
                bw2.write(Reads2 + "\n");
                bw2.write("+\n");
                StringBuffer sb2 = new StringBuffer();
                for (int j = 0; j < 150; j++) {
                    int score = -(j - 1) * (j - 1) * 8 / 22201 + 37;
//                    System.out.println(score);
                    String Q = getphred(score);
                    sb2.append(Q);
                }
                bw2.write(sb2.toString() + "\n");
            }

            bw.flush();
            bw.close();
            bw1.flush();
            bw1.close();
            bw2.flush();
            bw2.close();

            StringBuilder sb = new StringBuilder();
            sb.append("bgzip SE3724-150_R1.fq\n");
            sb.append("bgzip SE3724-150_R2.fq\n");
            sb.append("bgzip SE3724-350.fq\n");
            String command = sb.toString();
            System.out.println(command);
            try {
                File dir = new File(new File(arg).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
            System.out.println("Finished generate Reads");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void generateFastq(String[] args) {
        String infile = new File(args[1]).getAbsolutePath();
        String outfile = new File(args[0], "exons.fq").getAbsolutePath();
        HashSet<String> geneNameSet = new HashSet<>();
        HashMap<String, Integer> geneLengthMap = new HashMap<>();
        HashMap<String, Integer> geneEndNumber = new HashMap<>();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedReader brinfo = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String geneName = null;
//        GeneFeature gf = new GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
        try {
            int countline = 0;
            int geneLength = 0;
            int start = 0;
            int end = 0;
            int endNumber = 0;
            temp = brinfo.readLine();
            geneName = temp.split(" ")[0].split(">")[1];
//            System.out.println(geneName);
            while ((temp = brinfo.readLine()) != null) {
                countline++;
                if (temp.startsWith(">")) {
                    geneNameSet.add(geneName);
                    geneLengthMap.put(geneName, (countline - 2) * 70 + endNumber);
                    geneName = temp.split(" ")[0].split(">")[1];
                    countline = 0;
                }
                endNumber = temp.length();
            }
            geneNameSet.add(geneName);
            geneLengthMap.put(geneName, (countline - 1) * 70 + endNumber);
            brinfo.close();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    geneName = temp.split(" ")[0].split(">")[1];
                    geneLength = geneLengthMap.get(geneName);
                }
                if (geneLength < 350) {
                    int line = geneLength / 70;
                    int left = geneLength % 70;
                    if (left == 0) {
                        for (int i = 0; i < line; i++) {
                            br.readLine();
                        }
                    } else {
                        for (int i = 0; i <= line; i++) {
                            br.readLine();
                        }
                    }
                    continue;
                } else if (geneLength == 350) {
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 5; i++) {
                        temp = br.readLine();
                        sb.append(temp);
                    }
                    bw.write(geneName + "\t" + sb.toString() + "\n");
                } else {
                    Random r = new Random();
                    int totalline = 0;
                    if (geneLength % 70 == 0) {
                        totalline = geneLength / 70;
                    } else {
                        totalline = geneLength / 70 + 1;
                    }
                    int startsite = r.nextInt(geneLength - 350 + 1);
                    int line = startsite / 70;
                    int site1 = startsite % 70;
                    StringBuilder sb = new StringBuilder();
                    if (site1 != 0) {
                        for (int i = 0; i < line; i++) {
                            br.readLine();
                        }
                        temp = br.readLine();
                        sb.append(temp.substring(site1, temp.length()));
                        for (int i = 0; i < 4; i++) {
                            temp = br.readLine();
                            sb.append(temp);
                        }
                        temp = br.readLine();
//                        int lineend = 70 - site1 - 1;
                        sb.append(temp.substring(0, site1));
                        int leftline = totalline - line - 6;
                        for (int i = 0; i < leftline; i++) {
                            br.readLine();
                        }
                    } else {
                        for (int i = 0; i < line; i++) {
                            br.readLine();
                        }
                        for (int i = 0; i < 5; i++) {
                            temp = br.readLine();
                            sb.append(temp);
                        }
                        int leftline = totalline - line - 5;
                        for (int i = 0; i < leftline; i++) {
                            br.readLine();
                        }
                    }
                    String seq = sb.toString();
                    if (!seq.contains("N") && seq.length() == 350) {
                        bw.write(geneName + "\t" + sb.toString() + "\n");
                    }
                }
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println(geneName);
        }
    }

    public static void main(String[] args) {
        new simulationData(args);
    }
}
