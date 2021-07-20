package xiaohan.utils;

import java.io.BufferedReader;
import java.util.HashSet;

public class getHomoGene {

    public getHomoGene(String arg, String s) {
        this.ishomoGene(arg, s);
    }


    public static String ishomoGene(String infile, String gene) {
        BufferedReader br = IOUtils.getTextReader(infile);
        String temp = null;
        String[] temps = null;
        HashSet<String> Aset = new HashSet<>();
        HashSet<String> Bset = new HashSet<>();
        HashSet<String> Dset = new HashSet<>();
        String geneA = null;
        String geneB = null;
        String geneD = null;
        String out = null;
        try {
            while ((temp = br.readLine()) != null) {
                geneA = temps[0];
                geneB = temps[1];
                geneD = temps[2];
                Aset.add(geneA);
                Bset.add(geneB);
                Dset.add(geneD);
            }
            br.close();
            if (Aset.contains(gene)) {
                out = "A";
            }
            if (Bset.contains(gene)) {
                out = "B";
            }
            if (Dset.contains(gene)) {
                out = "D";
            } else {
                out = "null";
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return out;
    }

    public static void main(String[] args) {
        new getHomoGene(args[0], args[1]);
    }
}
