package xiaohan.Assembling;

import smile.math.Random;
import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.HashMap;
import java.util.HashSet;

/**
 * @ author: yxh
 * @ created: 2021-11-18 : 4:22 PM
 */
public class AssemblingMain {

    int length = 1000000;
    int readscount = 400;

    public AssemblingMain(String[] args) {
//        this.assemble();
        this.kmer(args);
    }

    public void kmer(String[] args) {
        int kl = Integer.parseInt(args[2]);
        String infile = args[0];
        String outfile = args[1];
        BufferedReader br = IOUtils.getTextReader(new File(infile).getAbsolutePath());
        String header = null;
        String read = null;
        String flag = null;
        String quality = null;
        BufferedWriter bw = IOUtils.getTextWriter(new File(outfile).getAbsolutePath());
        try {
            HashSet<String> kmers = new HashSet<>();
            while ((header = br.readLine()) != null) {
                read = br.readLine();
                flag = br.readLine();
                quality = br.readLine();
                for (int i = 0; i < read.length() - kl; i++) {
                    String temp = read.substring(i, i + kl);
                    kmers.add(temp);
                }
            }
            br.close();
            int[] count = new int[kmers.size()];
            HashMap<String, Integer> kmerIndex = new HashMap<>();
            String[] kmerArray = kmers.toArray(new String[0]);
            for (int i = 0; i < kmerArray.length; i++) {
                kmerIndex.put(kmerArray[i], i);
                count[i] = 0;
            }
            br = IOUtils.getTextReader(new File(infile).getAbsolutePath());
            while ((header = br.readLine()) != null) {
                read = br.readLine();
                flag = br.readLine();
                quality = br.readLine();
                for (int i = 0; i < read.length() - kl; i++) {
                    String temp = read.substring(i, i + kl);
                    count[kmerIndex.get(temp)] += 1;
                }
            }
            br.close();
            for (int i = 0; i < count.length; i++) {
                bw.write(kmerArray[i] + "\t" + count[i] + "\n");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void assemble() {
        String fa = this.generateSequence(length);
        String[] fq = this.getFastqs(fa, readscount);

    }

    public void isRelated(String fq1, String fq2) {
        String sub1 = fq1.substring(0, 20);
        String sub2 = fq2.substring(0, 20);
        if (fq1.contains(sub2)) {
//            return fq2,
        }
    }

    public String[] getFastqs(String fa, int readscount) {
        String[] fq = new String[readscount];
        Random r = new Random();
        int random = -1;
        for (int i = 0; i < readscount; i++) {
            random = r.nextInt(length - 25000);
            fq[i] = fa.substring(random, random + 25000);
        }
        return fq;
    }

    public String generateSequence(int length) {
        StringBuilder sb = new StringBuilder();
        Random r = new Random();
        for (int i = 0; i < length; i++) {
            int random = r.nextInt(99);
            sb.append(getBase(random));
        }
        return sb.toString();
    }

    public static String getBase(int random) {
        String base = null;
        if (random >= 0 && random < 28) {
            base = "G";
        } else if (random >= 28 && random < 56) {
            base = "C";
        } else if (random >= 56 && random < 78) {
            base = "A";
        } else {
            base = "T";
        }
        return base;
    }

    public static void main(String[] args) {
        new AssemblingMain(args);
    }
}