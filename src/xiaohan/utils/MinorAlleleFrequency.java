package xiaohan.utils;

import java.io.BufferedReader;
import java.util.HashMap;

/**
 * @author Xiaohan Yang
 */

public class MinorAlleleFrequency {
//
//    String infile = "";
//    String chr = "";

    public MinorAlleleFrequency(Object arg) {
        this.getMinorGenotype((Integer) arg);
        this.getMajorBase((String) arg);
        this.getMinorBase((String) arg);
        this.getMinorFrequency((String) arg);
    }

    public static HashMap<String, String> getMinorBase(String infile) {
        BufferedReader br = IOUtils.getTextReader(infile);
        HashMap<String, String> baseMap = new HashMap<>();
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("CHROM")) continue;
                temps = temp.split("\t");
                String site = temps[1];
                String base1 = temps[4].split(":")[0];
                String base2 = temps[5].split(":")[0];
                double frq1 = Double.parseDouble(temps[4].split(":")[1]);
                double frq2 = Double.parseDouble(temps[5].split(":")[1]);
                if (frq1 < frq2) {
                    baseMap.put(site, base1);
                } else if (frq1 >= frq2) {
                    baseMap.put(site, base2);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return baseMap;
    }

    public static HashMap<String, String> getMajorBase(String infile) {
        BufferedReader br = IOUtils.getTextReader(infile);
        HashMap<String, String> baseMap = new HashMap<>();
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("CHROM")) continue;
                temps = temp.split("\t");
                String site = temps[0] + "_" + temps[1];
                String base1 = temps[4].split(":")[0];
                String base2 = temps[5].split(":")[0];
                double frq1 = Double.parseDouble(temps[4].split(":")[1]);
                double frq2 = Double.parseDouble(temps[5].split(":")[1]);
                if (frq1 >= frq2) {
                    baseMap.put(site, base1);
                } else if (frq1 < frq2) {
                    baseMap.put(site, base2);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return baseMap;
    }

    public static HashMap<String, Double> getMinorFrequency(String infile) {
        BufferedReader br = IOUtils.getTextReader(infile);
        HashMap<String, Double> frqMap = new HashMap<>();
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("CHROM")) continue;
                temps = temp.split("\t");
//                String site = temps[0] + "_" + temps[1];
                String site = temps[1];
                double frq1 = Double.parseDouble(temps[4].split(":")[1]);
                double frq2 = Double.parseDouble(temps[5].split(":")[1]);
                if (frq1 < frq2) {
                    frqMap.put(site, frq1);
                } else if (frq1 >= frq2) {
                    frqMap.put(site, frq2);
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return frqMap;
    }

    public static HashMap<String, String> getMinorFrequency1(String infile) {
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        HashMap<String, String> frqMap = new HashMap<>();
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("CHROM")) continue;
                temps = temp.split("\t");
//                String site = temps[0] + "_" + temps[1];
                String site = temps[2];
                String frq = temps[3];
                frqMap.put(site, frq);
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return frqMap;
    }

    public static String[] getMinorGenotype(Integer s) {
        String[] genotype1 = {"0/1", "0/0"};
        String[] genotype2 = {"0/1", "1/1"};
        String[] genotype = new String[2];
        if (s == 0) {
            genotype = genotype1;
        } else if (s == 1) {
            genotype = genotype2;
        }
        return genotype;
    }

    public static void main(String[] args) {
        new MinorAlleleFrequency(args[0]);
    }
}
