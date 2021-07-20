package xiaohan.eQTL.pipline;

import xiaohan.utils.MathUtils;
import xiaohan.utils.RowTable;
import xiaohan.rareallele.GeneFeature;
import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class pheno {

    ///data2/junxu/dataTest/WEGA/S3/leaf
    String TissueDir = null;
    ///data2/xiaohan/genotype/hapscanner/output/S3/leaf
    String replacesample = null;
    String plate = null;
    String outputDir = null;
    String JarDir = "/data1/home/xiaohan/jar";


    //args1 : TissueDir ;    args2: replace sample
    public pheno(String[] args) {
        this.parseparameter(args);
        this.getMergedCountTable();
        this.DESeqnormalization();
        this.fixation();
        this.filtersampe();
//        this.changeName();
        this.countExpDonor02();
        this.SplitPhenoBychr();
        this.sortandbgzip();
    }

    public void parseparameter(String[] args) {
        System.out.println("*********Parsing paramter ******************************************************************");
//        if (!args.equals(2)) {
//            System.out.println("args1 : countTable tissue dir");
//            System.out.println("args2 : phenolist tissue dir ");
//        }
        this.plate = args[0];
        this.replacesample = args[1];
        HashSet<String> tissues = new HashSet<>();
        this.TissueDir = args[2];
//        this.TissueDir[1] = args[3];
        this.outputDir = new File("/data2/xiaohan/pheno/", plate).getAbsolutePath();
        File out = new File(outputDir);
        out.mkdir();
    }

    public void getMergedCountTable() {
        System.out.println("*********Merging Tables ********************************************************************");
        //get merged countTable
//        this.TissueDir = args[0];
        ArrayList<String> files = new ArrayList<>();
        String[] dirs = new File(TissueDir).list();
        if (plate.equals("S4leaf")) {
            for (int i = 0; i < dirs.length; i++) {
                String name = dirs[i].toString();
                String name1 = name.substring(name.length() - 2, name.length());
                System.out.println(name1);
                if (Integer.parseInt(name1) > 70) {
                    if (!dirs[i].endsWith("DS_Store") && !dirs[i].startsWith("countTable")) {
                        String file = new File(TissueDir, dirs[i] + "/countTable/" + dirs[i] + "_countResult.txt").getAbsolutePath();
                        files.add(file);
                    }
                } else continue;
            }
        } else {
            for (int i = 0; i < dirs.length; i++) {
                if (!dirs[i].endsWith("DS_Store") && !dirs[i].startsWith("countTable")) {
                    String file = new File(TissueDir, dirs[i] + "/countTable/" + dirs[i] + "_countResult.txt").getAbsolutePath();
                    files.add(file);
                }
            }
        }
        RowTable<String>[] ts = new RowTable[files.size()];
        for (int i = 0; i < files.size(); i++) {
            ts[i] = new RowTable<String>(files.get(i));
        }
        BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, plate + "countResult.txt").getAbsolutePath());
        try {
            StringBuilder sb = new StringBuilder();
            sb.append("Gene\t");
            for (int i = 0; i < ts.length; i++) {
                String[] header = ts[i].getHeader().toArray(new String[0]);
                for (int j = 1; j < header.length; j++) {
                    sb.append(header[j] + "\t");
                }
            }
            bw.write(sb.toString().replaceAll("\\s+$", ""));
            bw.newLine();
            for (int k = 0; k < ts[0].getRowNumber(); k++) {
                int rowNumber = k;
                String geneNumber = ts[0].getCell(k, 0);
                StringBuilder sb1 = new StringBuilder();
                sb1.append(geneNumber + "\t");
                for (int i = 0; i < ts.length; i++) {
                    RowTable t = ts[i];
                    for (int j = 1; j < t.getColumnNumber(); j++) {
                        sb1.append(t.getCell(rowNumber, j)).append("\t");
                    }
                }
                bw.write(sb1.toString().replaceAll("\\s+$", ""));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void DESeqnormalization() {
        System.out.println("*********Deseq normalization ***************************************************************");
        String infile = new File(this.outputDir, plate + "countResult.txt").getAbsolutePath();
        String tempfile = new File(this.outputDir, plate + "expressiontemp.txt").getAbsolutePath();
        try {
            File f = new File(tempfile);
            StringBuilder sb = new StringBuilder();
            sb.append("Rscript /data2/xiaohan/tensorQTL/scripts/DEseq_Normalization.R " + infile + " " + tempfile);
            String command = sb.toString();
            System.out.println(command);
            File dir = new File(new File(JarDir).getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void fixation() {
        System.out.println("*********Fixation **************************************************************************");
        String tempfile = new File(this.outputDir, plate + "expressiontemp.txt").getAbsolutePath();
        String exprfile = new File(this.outputDir, plate + "expression.txt").getAbsolutePath();
        try {
            BufferedReader br = IOUtils.getTextReader(tempfile);
            BufferedWriter bw = IOUtils.getTextWriter(exprfile);
            String temp = null;
            String[] temps = null;
            while ((temp = br.readLine()) != null) {
                if (!temp.startsWith("T")) {
                    bw.write("Gene\t");
                    bw.write(temp + "\n");
                    continue;
                }
                bw.write(temp + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
//            File f = new File(tempfile);
//            f.delete();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void filtersampe() {
        System.out.println("*********Filtering samples *****************************************************************");
        String file = new File(replacesample).getAbsolutePath();
        HashMap<String, String> RNAtoDNAmap = new HashMap<>();
        HashSet<String> RNASet = new HashSet<>();
        BufferedReader br = IOUtils.getTextReader(new File(file).getAbsolutePath());
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                RNASet.add(temps[0]);
                RNAtoDNAmap.put(temps[0], temps[1]);
            }
            br.close();

            BufferedReader br1 = IOUtils.getTextReader(new File(outputDir, plate + "expression.txt").getAbsolutePath());
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, plate + "expression_hapscanner.txt").getAbsolutePath());

            HashSet<Integer> indexSet = new HashSet<>();
            indexSet.add(0);
            StringBuilder sb = new StringBuilder();
            temp = br1.readLine();
            temps = temp.split("\t");
            if (temp.startsWith("Gene")) {
                sb.append("Gene\t");
                for (int i = 0; i < temps.length; i++) {
                    String name = null;
                    if (plate.equals("S6grain") && temps[i].equals("E039")) {
                        continue;
                    } else if (plate.equals("S6grain") && temps[i].equals("E338")) {
                        continue;
                    } else if (plate.equals("S7leaf") && temps[i].equals("E098")) {
                        continue;
                    } else if (plate.equals("S7leaf") && temps[i].equals("E100")) {
                        continue;
                    } else if (plate.equals("S7leaf") && temps[i].equals("E107")) {
                        continue;
                    } else if (plate.equals("S7leaf") && temps[i].equals("E113")) {
                        continue;
                    } else if (temps[i].startsWith("B18")) {
                        name = "RNA" + temps[i].replace("B18.", "").replace(".", "-");
                    } else {
                        name = "RNA" + temps[i].replace(".", "-");
                    }
                    if (RNASet.contains(name)) {
                        indexSet.add(i);
                        if (temps[i].endsWith(".2")) {
                            sb.append(RNAtoDNAmap.get(name) + "\t");
                        } else {
                            sb.append(RNAtoDNAmap.get(name) + "\t");
                        }
                    }
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            while ((temp = br1.readLine()) != null) {
                sb.setLength(0);
                temps = temp.split("\t");
                for (int i = 0; i < temps.length; i++) {
                    if (indexSet.contains(i)) {
                        sb.append(temps[i] + "\t");
                    }
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
//            RowTable<String> t = new RowTable<>(new File(outputDir, plate + "expression.txt").getAbsolutePath());
//            HashSet<String> newheader = new HashSet<>();
//            for (int i = 1; i < t.getHeader().size(); i++) {
//                String name = "RNA" + t.getColumnName(i).replace(".", "-");
//                System.out.println(name);
//                if (!RNASet.contains(name)) {
//                    System.out.println("Remove " + name);
//                    t.removeColumn(i);
//                }
//            }


//            for (int i = 0; i < t.getColumnNumber(); i++) {
//                System.out.println(t.getColumnName(i));
//            }

//            t.writeTextTable(new File(outputDir, plate + "expression_hapscannertemp.txt").getAbsolutePath(),IOFileFormat.Text);

            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void changeName() {
        String infile = new File(outputDir, plate + "expression_hapscannertemp.txt").getAbsolutePath();
        String outfile = new File(outputDir, plate + "expression_hapscanner.txt").getAbsolutePath();
        String infor = new File(replacesample).getAbsolutePath();
        String temp = null;
        String[] temps = null;
        HashMap<String, String> RNAtoDNAmap = new HashMap<>();
        BufferedReader br = IOUtils.getTextReader(new File(infor).getAbsolutePath());
        BufferedReader br1 = IOUtils.getTextReader(new File(infile).getAbsolutePath());
        BufferedWriter bw = IOUtils.getTextWriter(new File(outfile).getAbsolutePath());
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                System.out.println(temps[0]);
                RNAtoDNAmap.put(temps[0], temps[1]);
            }
            br.close();

            while ((temp = br1.readLine()) != null) {
                if (temp.startsWith("Gene")) {
                    temps = temp.split("\t");
                    StringBuilder sb = new StringBuilder();
                    sb.append("Gene\t");
                    for (int i = 1; i < temps.length; i++) {
                        String name = "RNA" + temps[i].replace(".", "-");
                        System.out.println(name);
                        sb.append(RNAtoDNAmap.get(name) + "\t");
                        System.out.println(RNAtoDNAmap.get(name));
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                    continue;
                }
                bw.write(temp);
                bw.newLine();
            }
            br1.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void countExpDonor02() {
        System.out.println("*********phenotype making ******************************************************************");
        String infileS = new File(outputDir, plate + "expression_hapscanner.txt").getAbsolutePath();
        String outfileS = new File(outputDir, plate + "expression_hapscanner_donor02.txt").getAbsolutePath();
        String outfileS1 = new File(outputDir, plate + "expression_hapscanner_donor02_zscore.txt").getAbsolutePath();
        GeneFeature gf = new GeneFeature("/data1/home/xiaohan/reference/wheat_v1.1_Lulab.gff3");
        try {
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Gene")) {
                    bw.write("Chr\tstart\tend\tID\t");
                    bw.write(temp.replace("Gene\t", ""));
                    bw.newLine();
                    bw1.write("Chr\tstart\tend\tID\t");
                    bw1.write(temp.replace("Gene\t", ""));
                    bw1.newLine();
                    continue;
                }
                temps = temp.split("\t");
                int threshold = (int) ((temps.length - 1) * 0.2);
//                System.out.println(threshold);
                int Expcount = 0;
                for (int i = 1; i < temps.length; i++) {
                    if (Double.parseDouble(temps[i]) != 0) {
                        Expcount++;
                    }
                }
                if (Expcount > threshold) {
                    String geneName = temps[0];
                    int index = gf.getGeneIndex(geneName);
                    int chr = gf.getGeneChromosome(index);
                    int start = -1;
                    if (gf.getGeneStrand(index) == 1) {
                        start = gf.getGeneStart(index);
                    } else {
                        start = gf.getGeneEnd(index);
                    }
                    int end = start + 1;
                    StringBuilder sb = new StringBuilder();
                    StringBuilder sb1 = new StringBuilder();
                    sb.append(chr + "\t" + start + "\t" + end + "\t" + geneName + "\t");
                    sb1.append(chr + "\t" + start + "\t" + end + "\t" + geneName + "\t");
                    double[] exp = new double[temps.length - 1];
                    double[] zscore = new double[temps.length - 1];
                    for (int i = 1; i < temps.length; i++) {
                        exp[i - 1] = Double.parseDouble(temps[i]);
                        sb.append(exp[i - 1] + "\t");
                    }
                    zscore = MathUtils.getzscore(exp);
                    for (int i = 0; i < zscore.length; i++) {
                        sb1.append(zscore[i] + "\t");
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", ""));
                    bw1.write(sb1.toString().replaceAll("\\s+$", ""));
                    bw.newLine();
                    bw1.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            bw1.flush();
            bw1.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void SplitPhenoBychr() {
        System.out.println("*********Split by chr **********************************************************************");
        String infileS1 = new File(outputDir, plate + "expression_hapscanner_donor02_zscore.txt").getAbsolutePath();
        String output = new File(outputDir, "DE_log2_bed").getAbsolutePath();
        File f = new File(output);
        f.mkdir();
        try {
            BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(infileS1);
            BufferedWriter[] bw = new BufferedWriter[43];
            BufferedWriter bw1 = pgl.infra.utils.IOUtils.getTextWriter(new File(output, "header.txt").getAbsolutePath());
            for (int i = 0; i < 43; i++) {
                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(output, plate + "pheno" + i + ".bed").getAbsolutePath());
            }
            String temp = null;
            String[] temps = null;
            temp = br.readLine();
            bw1.write("#"+temp);
            bw1.newLine();
            bw1.flush();
            bw1.close();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                int count = Integer.parseInt(temps[0]);
                bw[count].write(temp);
                bw[count].newLine();
            }
            for (int i = 0; i < 43; i++) {
                bw[i].flush();
                bw[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void sortandbgzip() {
        System.out.println("*********Sort and bgzip ********************************************************************");
        String input = new File(outputDir, "DE_log2_bed").getAbsolutePath();
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 43; i++) {
            int chr = i;
            nameSet.add(String.valueOf(chr));
        }
        nameSet.parallelStream().forEach(f -> {
            try {
                StringBuilder sb = new StringBuilder();
                sb.append("sort -k1,1n -k2,2n " + plate + "pheno" + f + ".bed");
                sb.append(" > ").append(plate + "pheno" + f + ".tempbed\n");
                sb.append("cat header.txt " + plate + "pheno" + f + ".tempbed");
                sb.append(" > ").append(plate + "pheno" + f + ".sorted.bed\n");
//                sb.append("rm " + plate + "pheno" + f + ".tempbed\n");
//                sb.append("rm " + plate + "pheno" + f + ".bed\n");
                sb.append("bgzip " + plate + "pheno" + f + ".sorted.bed\n");
                sb.append("tabix -p bed " + plate + "pheno" + f + ".sorted.bed.gz\n");
                String command = sb.toString();
                System.out.println(command);
                File dir = new File(new File(input).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public static void main(String[] args) {
        new pheno(args);
    }
}
