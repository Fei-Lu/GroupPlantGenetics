package xiaohan.eQTL;

import pgl.app.hapScanner.HapScanner;
import xiaohan.rareallele.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.concurrent.*;

public class covaraties {

    String number = "121120,35910,101373,72401,29060,2985,100575,74985,87007,81719,27649,8135,30143,71034,129106,102720,16228,7672,110193,36881,46901,51821,12885,2631,52736,57920,89566,61332,13361,7748,66198,27386,108871,86517,18530,2914,96461,73097,114553,66978,17611,9836";
    String[] numbers = number.split(",");
    String phenolist = null;
    String phenodir = null;
    String plate = null;

    public covaraties(String[] args) {
        this.parseparameter(args);
        this.PCA();
        this.peer();
        this.cov();
    }

    public void parseparameter(String[] args) {
        if (!args.equals(3)) {
            System.out.println("args1 : phenolist");
            System.out.println("args2 : phenodir ");
            System.out.println("args3 : plate");
        }
        phenolist = args[0];
        phenodir = args[1];
        plate = args[2];
    }


    public void PCA() {
        this.getsamplelist();
        this.getVCF();
        this.getheader();
        this.mergeVCF();
        this.calPCA();
    }

    private void getsamplelist() {
        StringBuilder sb = new StringBuilder();
        File file = new File(new File(phenodir, "geno").getAbsolutePath());
        file.mkdir();
        sb.append("cat ").append(phenolist).append(" | awk -F '\\t' '{print $2}' > ");
        sb.append(new File(phenodir, "geno/samplelist1.txt").getAbsolutePath()).append("\n");
        sb.append("cat ").append(new File(phenodir, "geno/samplelist1.txt").getAbsolutePath());
        sb.append(" | sort | uniq > ").append(new File(phenodir, "geno/samplelist.txt").getAbsolutePath()).append("\n");
        sb.append("rm ").append(new File(phenodir, "geno/samplelist1.txt").getAbsolutePath()).append("\n");
        String command = sb.toString();
        File dir = new File(new File(phenodir).getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void getVCF() {
        ExecutorService pool = Executors.newFixedThreadPool(10);
        for (int i = 0; i < 42; i++) {
            int chr = i + 1;
            StringBuilder sb = new StringBuilder();
            sb.append("bcftools view -S ").append(new File(phenodir, "geno/samplelist.txt").getAbsolutePath());
            sb.append(" /data2/junxu/genotypeMaf005/" + chr + ".360.vcf.gz | grep -v \"#\"  | shuf -n ");
//            sb.append(numbers[Integer.parseInt(chr)-1]);
            sb.append(numbers[i]);
            sb.append("| sort -k1,1n -k2,2n > ").append(new File(phenodir, "geno/chr" + chr + ".vcf").getAbsolutePath());
            File dir = new File(new File(phenodir).getAbsolutePath());
            String command = sb.toString();
            Command com = new Command(command,dir);
            Future<Command> chrom = pool.submit(com);
        }
        try {
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    class Command implements Callable<Command> {
        String command = null;
        File dir = null;
        public Command (String command, File dir){
            this.command = command;
            this.dir = dir;
        }
        @Override
        public Command call() throws Exception {
            try {
                System.out.println(command);
                String[] cmdarry1 = {"/bin/bash", "-c", command};
                Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
                p1.waitFor();
            }catch (Exception e){
                e.printStackTrace();
            }
            return this;
        }
    }

    private void getheader() {
        StringBuilder sb = new StringBuilder();
        sb.append("bcftools view -S ").append(new File(phenodir, "geno/samplelist.txt").getAbsolutePath());
        sb.append(" /data2/junxu/genotypeMaf005/1.360.vcf.gz -Ov | grep \"#\"  > header_chr1.txt");
        String command = sb.toString();
        File dir = new File(new File(phenodir).getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }


    private void mergeVCF() {
        StringBuilder sb = new StringBuilder();
        sb.append("cat chr1.vcf chr2.vcf chr3.vcf chr4.vcf chr5.vcf chr6.vcf chr7.vcf chr8.vcf chr9.vcf chr10.vcf chr11.vcf chr12.vcf chr13.vcf chr14.vcf chr15.vcf chr16.vcf chr17.vcf chr18.vcf chr19.vcf chr20.vcf chr21.vcf chr22.vcf chr23.vcf chr24.vcf chr25.vcf chr26.vcf chr27.vcf chr28.vcf chr29.vcf chr30.vcf chr31.vcf chr32.vcf chr33.vcf chr34.vcf chr35.vcf chr36.vcf chr37.vcf chr38.vcf chr39.vcf chr40.vcf chr41.vcf chr42.vcf > all.PCA.vcf\n");
        sb.append("cat header_chr1.txt all.PCA.vcf > allPCA.sort.vcf\n");
        String command = sb.toString();
        File dir = new File(new File(phenodir, "geno").getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void calPCA() {
        StringBuilder sb = new StringBuilder();
        sb.append("plink2 --chr-set 42 --threads 20 --out ").append(plate).append(" --pca 10 --vcf allPCA.sort.vcf &\n");
        String command = sb.toString();
        File dir = new File(new File(phenodir, "geno").getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void peer() {
        this.getpeerfile();
        this.calpeer();
    }

    private void getpeerfile() {
        StringBuilder sb = new StringBuilder();
        sb.append("cp " + plate + "expression_hapscanner_donor02_log2.txt " + plate + "expression_hapscanner_donor02_log2_forpeer.txt\n");
        sb.append("sed -i \'s/#//g\' ").append(plate).append("expression_hapscanner_donor02_log2_forpeer.txt\n");
        String command = sb.toString();
        File dir = new File(new File(phenodir).getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private void calpeer() {
        StringBuilder sb = new StringBuilder();
        sb.append("Rscript ./peer.R ").append(plate).append("expression_hapscanner_donor02_log2_forpeer.txt ./ \n");
        String command = sb.toString();
        File dir = new File(new File(phenodir).getAbsolutePath());
        try {
            System.out.println(command);
            String[] cmdarry1 = {"/bin/bash", "-c", command};
            Process p1 = Runtime.getRuntime().exec(cmdarry1, null, dir);
            p1.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void cov() {
        String prefix = plate;
//        String title = "/data2/xiaohan/genotype/hapscanner/output/S4/spike/Isec/phenolist.txt";
        String PCA = new File(phenodir, "geno/" + prefix + ".eigenvec").getAbsolutePath();
        String peer = new File(phenodir, "factors5.txt").getAbsolutePath();
        String output = new File(phenodir, prefix + "_cov.txt").getAbsolutePath();
        BufferedReader brPCA = IOUtils.getTextReader(PCA);
        BufferedReader brpeer = IOUtils.getTextReader(peer);
        BufferedWriter bw = IOUtils.getTextWriter(output);
        String temp = null;
        String[] temps = null;
        try {
            HashMap<String, String> PCAmap = new HashMap();
            while ((temp = brPCA.readLine()) != null) {
                if (temp.startsWith("E")) {
                    temps = temp.split("\t");
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 3; i++) {
                        sb.append(temps[i + 1] + "\t");
                    }
                    PCAmap.put(temps[0], sb.toString());
                }
            }
            brPCA.close();

            LinkedHashSet<String> nameSet = new LinkedHashSet<>();
            HashMap<String, String> PEERmap = new HashMap();
            while ((temp = brpeer.readLine()) != null) {
                if (temp.startsWith("E")) {
                    temps = temp.split("\t");
                    nameSet.add(temps[0].replace(".", "-"));
                    StringBuilder sb = new StringBuilder();
                    for (int i = 0; i < 5; i++) {
                        sb.append(temps[i + 2] + "\t");
                    }
                    PEERmap.put(temps[0].replace(".", "-"), sb.toString());
                }
            }
            brpeer.close();

            String[] namelist = nameSet.toArray(new String[nameSet.size()]);
            StringBuilder sb = new StringBuilder();
            sb.append("ID\t");
            for (int i = 0; i < namelist.length; i++) {
                sb.append(namelist[i] + "\t");
            }
            bw.write(sb.toString().replaceAll("\\s+$", ""));
            bw.newLine();
            sb.setLength(0);

            for (int i = 0; i < 3; i++) {
                int index = i + 1;
                sb.append("PC" + index + "\t");
                for (int j = 0; j < namelist.length; j++) {
                    System.out.println(PCAmap.get(namelist[j].substring(0, 4)));
                    sb.append(PCAmap.get(namelist[j].substring(0, 4)).split("\t")[i] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
                sb.setLength(0);
            }

            for (int i = 0; i < 5; i++) {
                int index = i + 1;
                sb.append("V" + index + "\t");
                for (int j = 0; j < namelist.length; j++) {
                    sb.append(PEERmap.get(namelist[j]).split("\t")[i] + "\t");
                }
                bw.write(sb.toString().replaceAll("\\s+$", ""));
                bw.newLine();
                sb.setLength(0);
            }
            bw.flush();
            bw.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) {

    }
}
