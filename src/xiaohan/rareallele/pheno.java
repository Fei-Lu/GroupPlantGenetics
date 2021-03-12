package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

public class pheno {
    String plate = "46";
    public pheno(String[] args){
        this.countExpDonor02();
        this.SplitPhenoBychr();
    }

    public void SplitPhenoBychr() {
        String plate = this.plate;
//        String infileS = "/Users/yxh/Documents/eQTL/data/pheno/"+plate+"expression_donor02.txt";
        String infileS1 = "/Users/yxh/Documents/eQTL/data/pheno/"+plate+"expression_donor02_log2.txt";
        String outputDir = "/Users/yxh/Documents/eQTL/data/pheno/"+plate;
//        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionboxcox7.txt";
//        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/splitexpression";
        try {
            BufferedReader br = pgl.infra.utils.IOUtils.getTextReader(infileS1);
            BufferedWriter[] bw = new BufferedWriter[43];
            BufferedWriter[] bw1 = new BufferedWriter[43];
            for (int i = 0; i < 43; i++) {
//                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".bed").getAbsolutePath());
//                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, "S7expression" + i + ".txt").getAbsolutePath());
                bw[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, plate +"pheno" + i + ".bed").getAbsolutePath());
                bw1[i] = pgl.infra.utils.IOUtils.getTextWriter(new File(outputDir, plate + "pheno" +i + ".txt").getAbsolutePath());
            }
            String temp = null;
            String[] temps = null;
            temp = br.readLine();
//            for (int i = 0; i < 43; i++) {
//                bw[i].write(temp);
//                bw[i].newLine();
//                bw1[i].write(temp);
//                bw1[i].newLine();
//            }
            while ((temp = br.readLine()) != null) {
                temps = temp.split("\t");
                int count = Integer.parseInt(temps[0]);
                bw[count].write(temp);
                bw[count].newLine();
                bw1[count].write(temp);
                bw1[count].newLine();
            }
            for (int i = 0; i < 43; i++) {
                bw[i].flush();
                bw[i].close();
                bw1[i].flush();
                bw1[i].close();
            }
            br.close();
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void countExpDonor02() {
        String plate = this.plate;
        String infileS = "/Users/yxh/Documents/eQTL/data/pheno/"+plate+"expression_hapscanner.txt";
        String outfileS = "/Users/yxh/Documents/eQTL/data/pheno/"+plate+"expression_donor02.txt";
        String outfileS1 = "/Users/yxh/Documents/eQTL/data/pheno/"+plate+"expression_donor02_log2.txt";
        GeneFeature gf = new GeneFeature("/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/refer/wheat_v1.1_Lulab.gff3");
        try {
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("Gene")) {
                    bw.write("#Chr\tstart\tend\tID\t");
                    bw.write(temp.replace("Gene\t",""));
                    bw.newLine();
                    bw1.write("#Chr\tstart\tend\tID\t");
                    bw1.write(temp.replace("Gene\t",""));
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
//                System.out.println(Expcount);
                if (Expcount > threshold) {
                    String geneName = temps[0];
                    int index = gf.getGeneIndex(geneName);
                    int chr = gf.getGeneChromosome(index);
                    int start = gf.getGeneStart(index);
                    int end = gf.getGeneEnd(index);
                    bw.write(chr+"\t"+start+"\t"+end+"\t"+geneName+"\t");
                    bw1.write(chr+"\t"+start+"\t"+end+"\t"+geneName+"\t");
                    StringBuilder sb = new StringBuilder();
                    StringBuilder sb2 = new StringBuilder();
                    for (int i = 1; i <= temps.length ; i++) {
                        double exp = Double.parseDouble(temps[i]);
//                        System.out.println(exp);
                        double log2exp = Math.log10(exp + 1) / Math.log10(2);

                        bw.write(exp+"\t");
                        bw1.write(log2exp+"\t");
                    }
                    double exp = Double.parseDouble(temps[temps.length-1]);
//                        System.out.println(exp);
                    double log2exp = Math.log10(exp + 1) / Math.log10(2);
                    bw.write(String.valueOf(exp));
                    bw1.write(String.valueOf(log2exp));
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

    public static void main (String[] args){
        new pheno(args);
    }
}
