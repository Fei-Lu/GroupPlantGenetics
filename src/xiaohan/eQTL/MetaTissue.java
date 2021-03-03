package xiaohan.eQTL;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;

public class MetaTissue {
    String javaDir = null;
    String shellDir = null;
    String MetasoftDir = null;
    String StudySigfile = null;
    String StudyNominalfile = null;
    String StudyOutputDir = null;
    String[] subDirS = {"Extracted", "Input", "Output", "Summary"};

    public MetaTissue(String[] args) {
        this.parseparameters(args);
        this.extractsig();
        this.extractall();
        this.prepareinput();
        this.metasoft();
        this.postprocess();
    }

    public void postprocess() {
        String outputfile = new File(StudyOutputDir, "output.txt").getAbsolutePath();
        String inputfile = new File(StudyOutputDir, "allextract.txt").getAbsolutePath();
        String temp = null;
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("python ");
            if (shellDir.endsWith("/")) {
                sb.append(shellDir);
            } else sb.append(shellDir).append("/");
            sb.append("metasoft_postprocess.py ");
            sb.append(outputfile + " ");
            sb.append(inputfile + " ");
            sb.append("meta_summary -o ");
            sb.append(new File(StudyOutputDir, subDirS[3]).getAbsolutePath());
            String command = sb.toString();
            System.out.println("-----------------------------------Metasoft summary--------------------------------------");
            System.out.println(command);
            File dir = new File(new File(StudyOutputDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void metasoft() {
        String temp = null;
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("python ");
            if (shellDir.endsWith("/")) {
                sb.append(shellDir);
            } else sb.append(shellDir).append("/");
            sb.append("run_metasoft.py ");
            if (MetasoftDir.endsWith("/")) {
                sb.append(MetasoftDir);
            } else sb.append(MetasoftDir).append("/");
            sb.append("Metasoft.jar ");
            sb.append(new File(StudyOutputDir, subDirS[1]).getAbsolutePath());
            sb.append("/metasoft.metasoft_input.txt.gz all ");
            sb.append("--pvalue_table ");
            sb.append(new File(MetasoftDir, "HanEskinPvalueTable.txt").getAbsolutePath());
            sb.append(" --seed 100 -o ");
            sb.append(new File(StudyOutputDir, subDirS[2]).getAbsolutePath());
            String command = sb.toString();
            System.out.println("-----------------------------------Run metasoft--------------------------------------");
            System.out.println(command);
            File dir = new File(new File(StudyOutputDir, subDirS[2]).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(StudyOutputDir, "output.txt").getAbsolutePath()));
            bw.write(new File(StudyOutputDir, subDirS[2]).getAbsolutePath());
            bw.write("/all.metasoft.txt.gz");
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void prepareinput() {
        String inputfile = new File(StudyOutputDir).getAbsolutePath() + "/allextract.txt";
        String temp = null;
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("python ");
            if (shellDir.endsWith("/")) {
                sb.append(shellDir);
            } else sb.append(shellDir).append("/");
            sb.append("metasoft_prepare_input.py ");
            sb.append(inputfile).append(" ");
            sb.append("metasoft -o ");
            sb.append(new File(StudyOutputDir, subDirS[1]).getAbsolutePath());
            sb.append(" --write_full ");
//            sb.append(" > logsigcombine.txt &");
            String command = sb.toString();
            System.out.println("-----------------------------------Prepare Meta input--------------------------------------");
            System.out.println(command);
            File dir = new File(new File(StudyOutputDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void extractall() {
        HashSet<Integer> numCores = new HashSet<>();
        String outfile = "allextract.txt";
        String temp = null;
        String[] temps = null;
        List<String> nominalList = new ArrayList<>();
        List<String> prefixList = new ArrayList<>();
        try {
            BufferedReader br1 = new BufferedReader(new FileReader(StudyNominalfile));
            BufferedWriter bw = new BufferedWriter(new FileWriter(new File(StudyOutputDir, outfile).getAbsolutePath()));
            int countline = 0;
            while ((temp = br1.readLine()) != null) {
                nominalList.add(temp);
                temps = temp.split("/");
//                System.out.println(temp);
//                System.out.println(temps.length);
                prefixList.add(temps[temps.length - 1].split("\\.")[0]);
//                System.out.println(temps[temps.length - 1].split("\\.")[0]);
                numCores.add(countline);
                bw.write(new File(StudyOutputDir, subDirS[0]).getAbsolutePath());
                bw.write("/"+prefixList.get(countline) + ".extracted_pairs.txt.gz\n");
                countline++;
            }
            br1.close();
            bw.flush();
            bw.close();
            if (nominalList.size() == prefixList.size()) {
                numCores.parallelStream().forEach(f -> {
                    try {
                        StringBuilder sb = new StringBuilder();
                        sb.append("python ");
                        if (shellDir.endsWith("/")) {
                            sb.append(shellDir);
                        } else sb.append(shellDir).append("/");
                        sb.append("extract_pairs.py ");
                        sb.append(nominalList.get(f)).append(" ");
                        if (StudyOutputDir.endsWith("/")) {
                            sb.append(StudyOutputDir);
                        } else sb.append(StudyOutputDir).append("/");
                        sb.append("meta-sig.combined_signifpairs.txt.gz ");
                        sb.append(prefixList.get(f) + " -o ");
                        if (StudyOutputDir.endsWith("/")) {
                            sb.append(StudyOutputDir);
                        } else sb.append(StudyOutputDir).append("/");
                        sb.append(subDirS[0]);
                        String command = sb.toString();
                        System.out.println("-----------------------------------Extract nominal pairs--------------------------------------");
                        System.out.println(command);
                        File dir = new File(new File(StudyOutputDir).getAbsolutePath());
                        String[] cmdarry = {"/bin/bash", "-c", command};
                        Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                        BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                        String temp1 = null;
                        while ((temp1 = br.readLine()) != null) {
                            System.out.println(temp1);
                        }
                        p.waitFor();
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                });
            } else {
                System.out.println("Can't find prefix of tissue..................................");
            }
        } catch (
                Exception e) {
            e.printStackTrace();
        }

    }

    public void extractsig() {
        String temp = null;
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("python ");
            if (shellDir.endsWith("/")) {
                sb.append(shellDir);
            } else sb.append(shellDir).append("/");
            sb.append("combine_signif_pairs.py ");
            sb.append(StudySigfile).append(" ");
            sb.append("meta-sig -o ");
            sb.append(StudyOutputDir);
            String command = sb.toString();
            System.out.println("-----------------------------------combind sig pairs--------------------------------------");
            System.out.println(command);
            File dir = new File(new File(StudyOutputDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            while ((temp = br.readLine()) != null) {
                System.out.println(temp);
            }
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void parseparameters(String[] infileS) {
        System.out.println("-----------------------------------Parsing parameters--------------------------------------");
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = new BufferedReader(new FileReader(infileS[0]));
            String temp = null;
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) continue;
                if (temp.isEmpty()) continue;
                pLineList.add(temp);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        if(pLineList.get(0).endsWith("/")) {
            this.shellDir = pLineList.get(0) + "shell";
        }
        else this.shellDir = pLineList.get(0) + "/shell";
        if(pLineList.get(0).endsWith("/")) {
            this.MetasoftDir = pLineList.get(0) + "metasoft";
        }else this.MetasoftDir = pLineList.get(0) + "/metasoft";
        this.StudySigfile = pLineList.get(1);
        this.StudyNominalfile = pLineList.get(2);
        this.StudyOutputDir = pLineList.get(3);
        System.out.println("Loading full......");
        System.out.println("-----------------------------------Building outputDirs--------------------------------------");
        for (int i = 0; i < this.subDirS.length; i++) {
            new File(this.StudyOutputDir, subDirS[i]).mkdir();
        }
        System.out.println("Building all......");
    }

    public static void main(String[] args) {
        new MetaTissue(args);
    }
}
