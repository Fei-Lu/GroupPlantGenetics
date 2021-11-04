package xiaohan.eQTL;

import xiaohan.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import static com.google.common.primitives.Shorts.max;
import static com.google.common.primitives.Shorts.min;

/**
 * GenotypeValidation
 *
 * @ author: yxh
 * @ created: 2021-10-12 : 9:23 PM
 */
public class genotypes {


    public genotypes(String[] args) {
//        this.validation(args);
        this.change(args);
    }

    public void change(String[] args) {
        String infile = args[0];
        String outfile = args[1];
        BufferedReader br;
        System.out.println("Beginning");
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(new File(infile).getAbsolutePath());
        } else {
            br = IOUtils.getTextReader(new File(infile).getAbsolutePath());
        }
        BufferedWriter bw = IOUtils.getTextWriter(new File(outfile).getAbsolutePath());
        int line = 18;
        int A = 0;
        int B = 0;
        int C = 0;
        String temp = null;
        String[] temps = null;
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    bw.write(temp);
                    bw.newLine();
                    continue;
                }
                temps = temp.split("\t");
                A = B = C = 0;
                for (int i = 0; i < temps.length; i++) {
                    if (temps[i].startsWith("0/0")) {
                        A++;
                    } else if (temps[i].startsWith("0/1")) {
                        B++;
                    } else if (temps[i].startsWith("1/1")) {
                        C++;
                    } else continue;
                }
                if (B <= 18) {
                    bw.write(temp);
                    bw.newLine();
                }else continue;
            }
            bw.flush();
            bw.close();
            StringBuilder sb = null;
            sb.append("bgzip "+new File(outfile).getAbsolutePath());
            String command = sb.toString();
//                        System.out.println(command);
            File dir = new File(new File("/data1/home/xiaohan").getAbsolutePath());
            String[] cmdarry = {"/bin/bash", "-c", command};
            Process p = Runtime.getRuntime().exec(cmdarry, null, dir);
            p.waitFor();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void validation(String[] args) {
        String infile = args[0];
        String infor = args[1];
        String outfile = args[2];
        BufferedReader br;
        BufferedReader brinfor;
        System.out.println("Beginning");
        if (infile.endsWith("gz")) {
            br = IOUtils.getTextGzipReader(new File(infile).getAbsolutePath());
        } else {
            br = IOUtils.getTextReader(new File(infile).getAbsolutePath());
        }
        if (infor.endsWith("gz")) {
            brinfor = IOUtils.getTextGzipReader(new File(infor).getAbsolutePath());
        } else {
            brinfor = IOUtils.getTextReader(new File(infor).getAbsolutePath());
        }
        BufferedWriter bw = IOUtils.getTextWriter(new File(outfile).getAbsolutePath());
        int line = 0;
        int A = 0;
        int B = 0;
        int C = 0;
        int sumA = 0;
        int sumB = 0;
        int sumC = 0;
        int sumB1 = 0;
        int sumB2 = 0;
        HashSet<String> sampleName = new HashSet<>();
        String temp = null;
        String[] temps = null;
        String[] tems = null;
        List<Integer> sampleIndex = new ArrayList<>();
        DecimalFormat dec = new DecimalFormat("0.00");
        try {
            while ((temp = brinfor.readLine()) != null) {
                temps = temp.split("\t");
                for (int i = 0; i < temps.length; i++) {
                    System.out.println(temps[i]);
                    if (temps[i].startsWith("E")) {
                        sampleName.add(temps[i]);
                    }
                }
            }
            line = sampleName.size();
            System.out.println(sampleName.size());
            brinfor.close();
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("##")) continue;
                if (temp.startsWith("#C")) {
                    temps = temp.split("\t");
                    for (int i = 0; i < temps.length; i++) {
                        if (sampleName.contains(temps[i])) {
                            sampleIndex.add(i);
                        }
                    }
                    continue;
                }
                temps = temp.split("\t");
                bw.write(temps[0] + "\t" + (Integer.parseInt(temps[1]) - 1) + "\t" + Integer.parseInt(temps[1]) + "\t");
                A = B = C = 0;
                sumA = sumB = sumB1 = sumB2 = sumC = 0;
                for (int i = 0; i < sampleIndex.size(); i++) {
                    tems = temps[sampleIndex.get(i)].split(":")[1].split(",");
                    if (temps[sampleIndex.get(i)].startsWith("0/0")) {
                        A++;
                        sumA += Integer.parseInt(tems[0]) + Integer.parseInt(tems[1]);
                    } else if (temps[sampleIndex.get(i)].startsWith("0/1")) {
                        B++;
                        sumB += Integer.parseInt(tems[0]) + Integer.parseInt(tems[1]);
                        sumB1 = sumB1 + min(Integer.parseInt(tems[0]), Integer.parseInt(tems[1]));
                        sumB2 = sumB2 + max(Integer.parseInt(tems[0]), Integer.parseInt(tems[1]));
                    } else if (temps[sampleIndex.get(i)].startsWith("1/1")) {
                        C++;
                        sumC += Integer.parseInt(tems[0]) + Integer.parseInt(tems[1]);
                    } else continue;
                }
                bw.write(A + "\t" + B + "\t" + C + "\t" + dec.format((B / (double) sampleIndex.size())) + "\t" + dec.format(sumA / (double) A) + "\t" + dec.format(sumB / (double) B) + "\t" + dec.format(sumC / (double) C) + "\t" + dec.format(sumB1 / (double) B) + "\t" + dec.format(sumB2 / (double) B) + "\n");
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public int min(int a, int b) {
        if (a <= b) {
            return a;
        } else {
            return b;
        }
    }

    public int max(int a, int b) {
        if (a >= b) {
            return a;
        } else {
            return b;
        }
    }

    public static void main(String[] args) {
        new genotypes(args);
    }
}
