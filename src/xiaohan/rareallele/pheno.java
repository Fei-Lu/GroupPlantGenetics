package xiaohan.rareallele;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

public class pheno {
    public pheno(String[] args){
        this.countExpDonor02();
    }

    public void countExpDonor02() {
        String infileS = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/expressionTable/original/DEnorm7_87chr1-42.txt";
        String outputDir = "/Users/yxh/Documents/RareAllele/004test/SiPASpipeline/data/S7/";
        try {
            String temp = null;
            String[] temps = null;
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(new File(outputDir, "donor08GeneName.txt").getAbsolutePath());
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#") || temp.startsWith("chr") || temp.startsWith("Chr")) {
                    continue;
                }
                temps = temp.split("\t");
                int Expcount = 0;
                for (int i = 4; i < temps.length; i++) {
                    if (Double.parseDouble(temps[i]) != 0) {
                        Expcount++;
                    }
                }
                if (Expcount > 68) {
                    bw.write(temps[3]);
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

    public static void main (String[] args){
        new pheno(args);
    }
}
