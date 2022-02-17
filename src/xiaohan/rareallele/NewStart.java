package xiaohan.rareallele;

import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.summa.SumTaxaDivergence;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import pgl.infra.window.SimpleWindow;
import xiaohan.Assembling.AssemblingMain;
import xiaohan.eQTL.GeneFeature;
import xiaohan.utils.Alogrithm;
import xiaohan.utils.IOUtils;
import xiaohan.utils.chrUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.InputStreamReader;
import java.util.*;
import java.util.concurrent.*;

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
//        this.outliersAllSNP(args);
//        this.annotation(args);
//        this.summary(args);
//        this.Version360(args);
//        this.splitByABD(args);
//        this.Others(args);
//        this.bedFilesToABD(args);
//        this.test(args);
//        this.pseudoDiploid(args);
//        this.ChrABDtoChr42(args);
//        this.Chr42toChrABD(args);
//        this.getIBS(args);
//        this.getExpression(args);
//        this.FilterIndel(args);
//        this.phasedVCF(args);
//        this.getExisting(args);
//        this.getSlidingMAF(args);
//        this.getGeneMAF(args);
        this.getGenomeCount(args);
    }
    
    private void getGenomeCount(String[] args){

        String[] chrABD = {"A","B","D"};
        // 1A,1B,1D
        int[] chrlength = {594102056,689851870,495453186,
                780798557,801256715,651852609,
                750843639,830829764,615552423,
                744588157,673617499,509857067,
                709773743,713149757,566080677,
                618079260,720988478,473592718,
                736706236,750620385,638686055
        };
        String temp = null;
        String[] temps = null;
        int pos = -1;
        StringBuilder sb = new StringBuilder();
        try {
            for (int i = 1; i <= 7; i++) {
                for (int j = 0; j < chrABD.length; j++) {
                    String chr = i + chrABD[i];
                    BufferedReader br = IOUtils.getInFile(new File(args[0],"chr"+i+chrABD[j]+".txt.gz").getAbsolutePath());
                    BufferedWriter[] bw = new BufferedWriter[5];
                    for (int k = 0; k < bw.length; k++) {
                        bw[i] = IOUtils.getOutFile(new File(args[1],"chr"+i+chrABD[j]+".count."+k+"txt").getAbsolutePath());
                    }
                    pgl.infra.window.SimpleWindow[] sc = new SimpleWindow[5];
                    for (int k = 0; k < sc.length; k++) {
                        sc[i] = new SimpleWindow(chrlength[(i-1)*3 + j], 100000, 100000);
                    }
                    while((temp = br.readLine())!=null){
                        temps = temp.split("\t");
                        pos = Integer.parseInt(temps[2]);
                        sc[0].addPositionCount(pos);
                        double maf = Double.parseDouble(temps[3]);
                        int index = getindex(maf);
                        if(index != -1) {
                            sc[i].addPositionCount(pos);
                        }
                    }
                    br.close();
                    for (int m = 0; m < sc.length; m++) {
                        int[] values = sc[m].getWindowValuesInt();
                        for (int k = 0; k < values.length ; k++) {
                            sb.setLength(0);
                            sb.append("k"+"\t"+values[i]+"\n");
                            bw[m].write(sb.toString());
                        }
                    }

                    for (int k = 0; k < bw.length; k++) {
                        bw[k].flush();
                        bw[k].close();
                    }
                }
            }
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    private int getindex(double maf){
        if(maf < 0.01){
            return 1;
        }else if(maf >= 0.01 && maf < 0.05){
            return 2;
        }else if(maf >= 0.05 && maf > 0.1){
            return 3;
        }else if(maf >= 0.01 && maf < 0.25){
            return 4;
        }else return -1;
    }

    public void getGeneMAF(String[] args) {
        String input = new File(args[0]).getAbsolutePath();
        String output = new File(args[1]).getAbsolutePath();
//        int windowsize = 100000;
//        int windowstep = 100000;
        String[] chrName = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
//        int[] chrLength = {594102056, 689851870, 495453186,
//                780798557, 801256715, 651852609,
//                750843639, 830829764, 615552423,
//                744588157, 673617499, 509857067,
//                709773743, 713149757, 566080677,
//                618079260, 720988478, 473592718,
//                736706236, 750620385, 638686055};
        int size1 = 1000;
        int size2 = 100;
        int[] upcount = new int[size1];
        int[] downcount = new int[size1];
        int[] genecount = new int[size2];
        double[] upfrq = new double[size1];
        double[] downfrq = new double[size1];
        double[] genefrq = new double[size2];

        GeneFeature gf = new GeneFeature("");
        for (int i = 0; i < chrName.length; i++) {
//            pgl.infra.window.SimpleWindow sv = new pgl.infra.window.SimpleWindow(chrLength[i], windowsize, windowstep);
//            pgl.infra.window.SimpleWindow sc = new pgl.infra.window.SimpleWindow(chrLength[i], windowsize, windowstep);
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getInFile(new File(input, "chr" + chrName[i] + ".txt").getAbsolutePath());
            BufferedWriter bw = IOUtils.getOutFile(new File(output, "chr" + chrName[i] + ".windowmaf.txt").getAbsolutePath());
            try {
                int position = -1;
                double value = -1.0;
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    position = Integer.parseInt(temps[0]);
                    value = Double.parseDouble(temps[1]);

//                    sv.addDoubleValue(position, value);
//                    sc.addPositionCount(position);
                }
                br.close();
//                double[] values = sv.getWindowValuesDouble();
//                int[] counts = sc.getWindowValuesInt();
//                for (int j = 0; j < values.length; j++) {
//                    int w = j * windowsize;
//                    bw.write(w + "\t" + values[i] / counts[i] + "\n");
//                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getSlidingMAF(String[] args) {
        String input = new File(args[0]).getAbsolutePath();
        String output = new File(args[1]).getAbsolutePath();
        int windowsize = 100000;
        int windowstep = 100000;
        String[] chrName = {"1A", "1B", "1D", "2A", "2B", "2D", "3A", "3B", "3D", "4A", "4B", "4D", "5A", "5B", "5D", "6A", "6B", "6D", "7A", "7B", "7D"};
        int[] chrLength = {594102056, 689851870, 495453186,
                780798557, 801256715, 651852609,
                750843639, 830829764, 615552423,
                744588157, 673617499, 509857067,
                709773743, 713149757, 566080677,
                618079260, 720988478, 473592718,
                736706236, 750620385, 638686055};
        for (int i = 0; i < chrName.length; i++) {
            pgl.infra.window.SimpleWindow sv = new pgl.infra.window.SimpleWindow(chrLength[i], windowsize, windowstep);
            pgl.infra.window.SimpleWindow sc = new pgl.infra.window.SimpleWindow(chrLength[i], windowsize, windowstep);
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getInFile(new File(input, "chr" + chrName[i] + ".txt").getAbsolutePath());
            BufferedWriter bw = IOUtils.getOutFile(new File(output, "chr" + chrName[i] + ".windowmaf.txt").getAbsolutePath());
            try {
                int position = -1;
                double value = -1.0;
                while ((temp = br.readLine()) != null) {
                    temps = temp.split("\t");
                    position = Integer.parseInt(temps[0]);
                    value = Double.parseDouble(temps[1]);
                    sv.addDoubleValue(position, value);
                    sc.addPositionCount(position);
                }
                br.close();
                double[] values = sv.getWindowValuesDouble();
                int[] counts = sc.getWindowValuesInt();
                for (int j = 0; j < values.length; j++) {
                    int w = j * windowsize;
                    bw.write(w + "\t" + values[i] / counts[i] + "\n");
                }
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }

    public void getExisting(String[] args) {
        String infor = new File(args[0]).getAbsolutePath();
        String infile = new File(args[1]).getAbsolutePath();
        String outfile = new File(args[2]).getAbsolutePath();

        BufferedReader brinfor = IOUtils.getInFile(infor);
        BufferedReader br = IOUtils.getInFile(infile);
        BufferedWriter bw = IOUtils.getOutFile(outfile);

        String temp = null;
        String[] temps = null;
        HashSet<String> nameSet = new HashSet<>();

        try {
            while ((temp = brinfor.readLine()) != null) {
                if (temp.startsWith("EID")) continue;
                temps = temp.split("\t");
                nameSet.add(temps[1]);
            }
            brinfor.close();

            System.out.println(nameSet.size());

            int[] index = new int[346];

            int countline = 0;
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 5000 == 0) {
                    System.out.println(countline);
                }
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    int l = 0;
                    for (int i = 9; i < temps.length; i++) {
                        if (nameSet.contains(temps[i])) {
                            index[l] = i;
                            l++;
                        }
                    }
                    continue;
                }
                temps = temp.split("\t");

                for (int i = 0; i < index.length; i++) {
                    if (temps[index[i]].startsWith("0/1") || temps[index[i]].startsWith("1/1")) {
                        int end = Integer.parseInt(temps[1]);
                        int start = end - 1;
                        bw.write(temps[0] + "\t" + start + "\t" + end + "\t" + "1\n");
                        break;
                    }
                }
            }
            br.close();
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void phasedVCF(String[] args) {
        String input = new File(args[0]).getAbsolutePath();
        String output = new File(args[1]).getAbsolutePath();
        String inforDir = new File(args[2]).getAbsolutePath();
        int threadsNum = Integer.parseInt(args[3]);
        try {
            ExecutorService pool = Executors.newFixedThreadPool(threadsNum);
            for (int i = 0; i < 42; i++) {
                int chr = i + 1;
                String infile = new File(input, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf.gz").getAbsolutePath();
                String outfile = new File(output, "chr" + PStringUtils.getNDigitNumber(3, chr) + ".vcf.gz").getAbsolutePath();
                String infor = new File(inforDir, "chr" + PStringUtils.getNDigitNumber(3, chr) + "_secer_hv_ancestral.txt.gz").getAbsolutePath();
                Command cd = new Command(infile, outfile, infor);
                Future<Command> f = pool.submit(cd);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.SECONDS);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    class Command implements Callable<Command> {
        String infile = null;
        String outfile = null;
        String infor = null;

        public Command(String infile, String outfile, String infor) {
            this.infile = infile;
            this.outfile = outfile;
            this.infor = infor;
        }

        @Override
        public Command call() throws Exception {
            BufferedReader br = IOUtils.getInFile(infile);
            BufferedWriter bw = IOUtils.getOutFile(outfile);
            BufferedReader brinfo = IOUtils.getInFile(infor);
            String temp = null;
            String[] temps = null;
            StringBuilder sb = new StringBuilder();

            String reference = null;
            String alternative = null;
            try {
                int countline = 0;
                LinkedList<Long> positions = new LinkedList<>();
                LinkedList<String> ancestor = new LinkedList<>();
                while ((temp = brinfo.readLine()) != null) {
                    countline++;
                    if (countline % 500 == 0) {
                        System.out.println("Pasring infor : " + countline);
                    }
                    if (temp.startsWith("chr")) continue;
                    temps = temp.split("\t");
                    positions.add(Long.parseLong(temps[1]));
                    ancestor.add(temps[2]);
                }
                brinfo.close();
                countline = 0;

                Long[] poslist = positions.toArray(new Long[positions.size()]);
                String[] anc = ancestor.toArray(new String[ancestor.size()]);
                int ref = -1;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 500 == 0) {
                        System.out.println("Pasring VCF : " + countline);
                    }
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    int index = Alogrithm.BinarySearch(poslist, Long.parseLong(temps[1]), 0, poslist.length - 1);
                    if (index == -1) continue;
                    if (!anc[index].equals(temps[3]) && !anc[index].equals(temps[4])) continue;
                    if (anc[index] == temps[3]) {
                        reference = temps[3];
                        alternative = temps[4];
                    } else {
                        reference = temps[4];
                        alternative = temps[3];
                    }
                    sb.setLength(0);
                    sb.append(temps[0] + "\t" + temps[1] + "\t" + temps[2] + "\t" + reference + "\t" + alternative + "\t" + ".\t.\t.\tGT\t");
                    for (int i = 9; i < temps.length; i++) {
                        if (anc[index].equals(temps[3])) {
                            ref = 0;
                        } else {
                            ref = 1;
                        }
                        String phase1 = phase(temps[i].substring(0, 3), ref);
                        sb.append(phase1 + "\t");
                        String out = sb.toString();
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
                }
                br.close();
                bw.flush();
                bw.close();
            } catch (Exception e) {
                e.printStackTrace();
            }
            return null;
        }
    }

    public void FilterIndel(String[] args) {

        String inputDir = new File(args[0]).getAbsolutePath();
        String outputDir = new File(args[1]).getAbsolutePath();
        HashSet<String> chrSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            System.out.println(chr);
            chrSet.add(PStringUtils.getNDigitNumber(3, chr));
        }
        chrSet.parallelStream().forEach(f -> {
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getInFile(new File(inputDir, "chr" + f + ".vcf.gz").getAbsolutePath());
            BufferedWriter bw = IOUtils.getOutFile(new File(outputDir, "chr" + f + ".vcf.gz").getAbsolutePath());
            int countline = 0;
            try {
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 5000 == 0) {
                        System.out.println(countline);
                    }
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    if (temps[4].equals("A") && temps[4].equals("T") && temps[4].equals("G") && temps[4].equals("C")) {
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

    public void getExpression(String[] args) {
        String infile = new File(args[0]).getAbsolutePath();
        String outfile = new File(args[1]).getAbsolutePath();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        int samples = 0;
        try {
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                if (temp.startsWith("Gene")) {
                    samples = temps.length - 1;
                    continue;
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void getIBS(String[] args) {
        System.out.println("This is getting IBS matrix *****************************************************************");
        String infileDir = new File(args[0]).getAbsolutePath();
        String infileS1 = new File(infileDir, "RNA/RNAall.vcf.gz").getAbsolutePath();
        System.out.println(infileS1);
        String infileS2 = new File(infileDir, "DNA/DNAall.vcf.gz").getAbsolutePath();
        String ibsOutfileS = new File(infileDir, "Summary/check.txt").getAbsolutePath();
        GenotypeGrid g1 = new GenotypeGrid(infileS1, GenoIOFormat.VCF_GZ);
        System.out.println(g1 instanceof GenotypeGrid);
        GenotypeGrid g2 = new GenotypeGrid(infileS2, GenoIOFormat.VCF_GZ);
        GenotypeGrid g = GenotypeOperation.mergeGenotypesByTaxon(g1, g2);
        SumTaxaDivergence std = new SumTaxaDivergence(g);
        std.writeDxyMatrix(ibsOutfileS, IOFileFormat.Text);
        g.getIBSDistanceMatrix();
    }

    public void test(String[] args) {
        String ABD = "6D";
        String chr42 = "14";
        String chr21 = "15";
        int pos1 = 11;
        int pos2 = 71049525;
//        System.out.println(chrUtils.getChrABDtoChr42(ABD,pos1));
//        System.out.println(chrUtils.getChrABDtoChr42(ABD, pos2));
//        System.out.println(chrUtils.getChr21toChrABD(chr21));
//        System.out.println(chrUtils.getChr42toChrABD(chr42));
//        System.out.println(chrUtils.getChrABDpostoChr42pos(ABD,pos1));
//        System.out.println(chrUtils.getChrABDpostoChr42pos(ABD, pos2));
//        System.out.println(chrUtils.getChr42postoChrABDpos(chr42,pos1));
//        System.out.println(chrUtils.getChr42postoChrABDpos(chr42,pos2));
    }

    public void Chr42toChrABD(String[] args) {
        BufferedReader br = IOUtils.getInFile(new File(args[0]).getAbsolutePath());
        BufferedWriter bw = IOUtils.getOutFile(new File(args[1]).getAbsolutePath());
        int[] index = new int[3];
        index[0] = Integer.parseInt(args[2]);
        index[1] = Integer.parseInt(args[3]);
        index[2] = -1;
        if (args.length == 5) {
            index[2] = Integer.parseInt(args[4]);
        }
        String temp = null;
        String[] temps = null;
        try {
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                sb.setLength(0);
                String chr = null;
                for (int j = 0; j < temps.length; j++) {
                    if (j == index[0]) {
                        chr = chrUtils.getChr42toChrABD(temps[index[0]].replace("chr", ""));
                        sb.append("chr" + chr + "\t");
                        continue;
                    }
                    if (j == index[1]) {
                        sb.append(chrUtils.getChr42postoChrABDpos(temps[index[0]].replace("chr", ""), Integer.parseInt(temps[index[1]])) + "\t");
                        continue;
                    }
                    if (j == index[2] && index[2] != -1) {
                        sb.append(chrUtils.getChr42postoChrABDpos(temps[index[0]].replace("chr", ""), Integer.parseInt(temps[index[2]])) + "\t");
                        continue;
                    }
                    sb.append(temps[j] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void ChrABDtoChr42(String[] args) {
        BufferedReader br = IOUtils.getInFile(new File(args[0]).getAbsolutePath());
        BufferedWriter bw = IOUtils.getOutFile(new File(args[1]).getAbsolutePath());
        int[] index = new int[3];
        index[0] = Integer.parseInt(args[2]);
        index[1] = Integer.parseInt(args[3]);
        index[2] = -1;
        if (args.length == 5) {
            index[2] = Integer.parseInt(args[4]);
        }
        String temp = null;
        String[] temps = null;
        try {
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                sb.setLength(0);
                int chr = -1;
                for (int j = 0; j < temps.length; j++) {
                    if (j == index[0]) {
                        chr = chrUtils.getChrABDtoChr42Int(temps[index[0]].replace("chr", ""), Integer.parseInt(temps[index[1]]));
                        sb.append(chr + "\t");
                        continue;
                    }
                    if (j == index[1]) {
                        sb.append(chrUtils.getChrABDpostoChr42pos(temps[index[0]].replace("chr", ""), Integer.parseInt(temps[index[1]])) + "\t");
                        continue;
                    }
                    if (j == index[2] && index[2] != -1) {
                        sb.append(chrUtils.getChrABDpostoChr42pos(temps[index[0]].replace("chr", ""), Integer.parseInt(temps[index[2]])) + "\t");
                        continue;
                    }
                    sb.append(temps[j] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void pseudoDiploid(String[] args) {
        String infile = new File(args[0]).getAbsolutePath();
        String outfile = new File(args[1]).getAbsolutePath();
        String infor = new File(args[2]).getAbsolutePath();
        int goalNumber = Integer.parseInt(args[3]);
//        String infile = "/Users/yxh/Desktop/test.vcf";
//        String outfile = "/Users/yxh/Desktop/out.vcf";
//        String infor = "/Users/yxh/Desktop/anc.txt";

        BufferedReader br = IOUtils.getTextGzipReader(infile);
        BufferedReader brinfo = IOUtils.getTextGzipReader(infor);
        BufferedWriter bw = IOUtils.getTextGzipWriter(outfile);
        String temp = null;
        String[] temps = null;
        int sampleNumber = 0;

//        int goalNumber = 100;
        Random r = new Random();
        String reference = null;
        String alternative = null;
        int ref = -1;
        try {
            int countline = 0;
            StringBuilder sb = new StringBuilder();
            LinkedList<Long> positions = new LinkedList<>();
            LinkedList<String> ancestor = new LinkedList<>();
            while ((temp = brinfo.readLine()) != null) {
                countline++;
                if (countline % 500 == 0) {
                    System.out.println("Pasring infor : " + countline);
                }
                if (temp.startsWith("chr")) continue;
                temps = temp.split("\t");
                positions.add(Long.parseLong(temps[1]));
                ancestor.add(temps[2]);
            }
            brinfo.close();
            countline = 0;

            Long[] poslist = positions.toArray(new Long[positions.size()]);
            String[] anc = ancestor.toArray(new String[ancestor.size()]);
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 500 == 0) {
                    System.out.println("Pasring VCF : " + countline);
                }
                if (temp.startsWith("##")) {
                    bw.write(temp + "\n");
                    continue;
                }
                if (temp.startsWith("#C")) {
                    sb.setLength(0);
                    temps = temp.split("\t");
                    sampleNumber = temps.length - 9;
                    for (int i = 0; i < 9; i++) {
                        sb.append(temps[i] + "\t");
                    }
                    for (int i = 1; i <= goalNumber; i++) {
                        sb.append("PD" + PStringUtils.getNDigitNumber(4, i) + "\t");
                    }
                    bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
                    continue;
                }
                temps = temp.split("\t");
                int index = Alogrithm.BinarySearch(poslist, Long.parseLong(temps[1]), 0, poslist.length - 1);
                if (index == -1) continue;
                if (!anc[index].equals(temps[3]) && !anc[index].equals(temps[4])) continue;
                if (anc[index] == temps[3]) {
                    reference = temps[3];
                    alternative = temps[4];
                } else {
                    reference = temps[4];
                    alternative = temps[3];
                }
                sb.setLength(0);
                sb.append(temps[0] + "\t" + temps[1] + "\t" + temps[2] + "\t" + reference + "\t" + alternative + "\t" + ".\t.\t.\tGT\t");
                for (int i = 0; i < goalNumber; i++) {
                    if (anc[index].equals(temps[3])) {
                        ref = 0;
                    } else {
                        ref = 1;
                    }
                    String phase1 = phase(temps[i].substring(0, 3), ref);
                    sb.append(phase1 + "\t");
                    String out = sb.toString();
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static String phase(String gt, int ref) {
        if (gt.equals("0/1") || gt.equals("./.")) return ".|.";
        if (ref == 0) {
            return gt.replace("/", "|");
        } else {
            if (gt.equals("0/0")) {
                return gt.replace("/", "|").replace("0", "1");
            } else return gt.replace("/", "|").replace("1", "0");
        }

    }

    public void test() {
        String temp = "chr2\t309688\t309689";
        String[] temps = temp.split("\t");
        String chr = chrUtils.getChr42toChrABD(temps[0]);
        System.out.println(chr);
        int start = chrUtils.getChr42postoChrABDpos(temps[0], Integer.parseInt(temps[1]));
        int end = chrUtils.getChr42postoChrABDpos(temps[0], Integer.parseInt(temps[2]));
        System.out.println(chr + "\t" + start + "\t" + end);
    }

    public void bedFilesToABD(String[] args) {
        String infile = new File(args[0]).getAbsolutePath();
        String outfile = new File(args[1]).getAbsolutePath();
        BufferedReader br = IOUtils.getTextReader(infile);
        BufferedWriter bw = IOUtils.getTextWriter(outfile);
        String temp = null;
        String[] temps = null;
        String chr = null;
        int start = -1;
        int end = -1;
        try {
            StringBuilder sb = new StringBuilder();
            int countline = 0;
            while ((temp = br.readLine()) != null) {
                countline++;
                if (countline % 50 == 0) {
                    System.out.println(countline);
                }
                sb.setLength(0);
                temps = temp.split("\t");
                chr = chrUtils.getChr42toChrABD(temps[0]);
                start = chrUtils.getChr42postoChrABDpos(temps[0], Integer.parseInt(temps[1]));
                end = chrUtils.getChr42postoChrABDpos(temps[0], Integer.parseInt(temps[2]));
                sb.append(chr + "\t" + start + "\t" + end + "\t");
                for (int i = 3; i < temps.length; i++) {
                    sb.append(temps[i] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", "") + "\n");
            }
            br.close();
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void Others(String[] args) {
        new AssemblingMain(args);
    }

    public void splitByABD(String[] args) {
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
                    String CHR = temps[index].substring(8, 9);
                    if (CHR.equals("A")) {
                        bw[0].write(temp + "\n");
                    } else if (CHR.equals("B")) {
                        bw[1].write(temp + "\n");
                    } else {
                        bw[2].write(temp + "\n");
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

    public void outliersAllSNP(String[] args) {
        String infile = args[0];
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
            brinfo = IOUtils.getTextReader(new File(infile).getAbsolutePath());
            br = IOUtils.getTextGzipReader(new File(output, "all_rare_variants_SNPs_chr42.txt.gz").getAbsolutePath());
            bw1 = IOUtils.getTextGzipWriter(new File(output, "Underoutliers_rare_variants_SNPs_chr42.txt.gz").getAbsolutePath());
            bw2 = IOUtils.getTextGzipWriter(new File(output, "Overoutliers_rare_variants_SNPs_chr42.txt.gz").getAbsolutePath());
            bw3 = IOUtils.getTextGzipWriter(new File(output, "Nonoutliers_rare_variants_SNPs_chr42.txt.gz").getAbsolutePath());
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
//                if (gf.getGeneStrand(i) == 1) {
//                    start = gf.getGeneStart(i) - window - 1;
//                    end = start + 1;
                start = gf.getGeneStart(i);
                end = gf.getGeneEnd(i);
//                } else {
//                    start = gf.getGeneEnd(i) + window - 1;
//                    end = start + 1;
//                    start = gf.getGeneEnd(i);
//                    end = gf.getGeneStart(i);
//                }
                if (start >= 0 && end >= 0) {
                    bw.write("chr" + chr + "\t" + start + "\t" + end + "\t" + gf.getGeneStrand(i) + "\t" + genelist[i] + "\n");
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