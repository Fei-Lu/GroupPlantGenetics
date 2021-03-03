package xiaohan.eQTL;

import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashSet;

public class vcf {
    public vcf(String[] args){
        this.getFilteredVCF();
    }

    public void getFilteredVCF() {
        System.out.println("This is the beginning of the program............................");
        String inputdir = "/data2/junxu/genotypeMaf005";
        String outputdir = "/data2/xiaohan/genotype/genotype_eQTL";
        HashSet<String> nameSet = new HashSet<>();
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            nameSet.add(String.valueOf(chr));
            System.out.println(chr);
        }
        nameSet.stream().forEach(f -> {
            try {
                System.out.println("-----------------------This is calculating file : chr " + f + "------------------------------");
                BufferedReader br = IOUtils.getTextGzipReader(new File(inputdir, f + ".346.B18.recode.vcf.gz").getAbsolutePath());
                BufferedWriter bw = IOUtils.getTextWriter(new File(outputdir, f + ".346.B18.recode.vcf").getAbsolutePath());
                String temp = null;
                String[] temps = null;
                String[] tems = null;
                int countline = 0;
                while ((temp = br.readLine()) != null) {
                    countline++;
                    if (countline % 50000 == 0) {
                        System.out.println("This is counting line : " + countline);
                    }
                    if (temp.startsWith("#")) {
                        bw.write(temp + "\n");
                        continue;
                    }
                    temps = temp.split("\t");
                    int number = 0;
                    for (int i = 0; i < temps.length - 9; i++) {
                        tems = temps[i + 9].split(":");
                        if (tems[0].equals("0/1")) number++;
                    }
                    if (number <= 43) {
                        bw.write(temp + "\n");
                    } else continue;
                }
                br.close();
                bw.flush();
                bw.close();
                StringBuilder sb = new StringBuilder();
                sb.append("bgzip " + f + ".346.B18.recode.vcf");
                String command = sb.toString();
                File dir = new File(new File(outputdir).getAbsolutePath());
                String[] cmdarry = {"/bin/bash", "-c", command};
                Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
                p.waitFor();
            } catch (Exception e) {
                e.printStackTrace();
            }
        });
    }

    public static void main (String[] args){
        new vcf(args);
    }
}
