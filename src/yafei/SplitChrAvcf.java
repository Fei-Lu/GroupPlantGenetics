package yafei;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * @author Yafei Guo
 * @create 2020-03-02 11:52 PM
 */
public class SplitChrAvcf {
    public static void main(String[] args) throws Exception {
        File writename1 = new File("/data1/home/yafei/all_vcf/chr1A.hmp.txt");
        BufferedWriter out1 = new BufferedWriter(new FileWriter(writename1));
        File writename2 = new File("/data1/home/yafei/all_vcf/chr1B.hmp.txt");
        BufferedWriter out2 = new BufferedWriter(new FileWriter(writename2));
        File writename3 = new File("/data1/home/yafei/all_vcf/chr1D.hmp.txt");
        BufferedWriter out3 = new BufferedWriter(new FileWriter(writename3));

        File writename4 = new File("/data1/home/yafei/all_vcf/chr2A.hmp.txt");
        BufferedWriter out4 = new BufferedWriter(new FileWriter(writename4));
        File writename5 = new File("/data1/home/yafei/all_vcf/chr2B.hmp.txt");
        BufferedWriter out5 = new BufferedWriter(new FileWriter(writename5));
        File writename6 = new File("/data1/home/yafei/all_vcf/chr2D.hmp.txt");
        BufferedWriter out6 = new BufferedWriter(new FileWriter(writename6));

        File writename7 = new File("/data1/home/yafei/all_vcf/chr3A.hmp.txt");
        BufferedWriter out7 = new BufferedWriter(new FileWriter(writename7));
        File writename8 = new File("/data1/home/yafei/all_vcf/chr3B.hmp.txt");
        BufferedWriter out8 = new BufferedWriter(new FileWriter(writename8));
        File writename9 = new File("/data1/home/yafei/all_vcf/chr3D.hmp.txt");
        BufferedWriter out9 = new BufferedWriter(new FileWriter(writename9));

        File writename10 = new File("/data1/home/yafei/all_vcf/chr4A.hmp.txt");
        BufferedWriter out10 = new BufferedWriter(new FileWriter(writename10));
        File writename11 = new File("/data1/home/yafei/all_vcf/chr4B.hmp.txt");
        BufferedWriter out11 = new BufferedWriter(new FileWriter(writename11));
        File writename12 = new File("/data1/home/yafei/all_vcf/chr4D.hmp.txt");
        BufferedWriter out12 = new BufferedWriter(new FileWriter(writename12));

        File writename13 = new File("/data1/home/yafei/all_vcf/chr5A.hmp.txt");
        BufferedWriter out13 = new BufferedWriter(new FileWriter(writename13));
        File writename14 = new File("/data1/home/yafei/all_vcf/chr5B.hmp.txt");
        BufferedWriter out14 = new BufferedWriter(new FileWriter(writename14));
        File writename15 = new File("/data1/home/yafei/all_vcf/chr5D.hmp.txt");
        BufferedWriter out15 = new BufferedWriter(new FileWriter(writename15));

        File writename16 = new File("/data1/home/yafei/all_vcf/chr6A.hmp.txt");
        BufferedWriter out16 = new BufferedWriter(new FileWriter(writename16));
        File writename17 = new File("/data1/home/yafei/all_vcf/chr6B.hmp.txt");
        BufferedWriter out17 = new BufferedWriter(new FileWriter(writename17));
        File writename18 = new File("/data1/home/yafei/all_vcf/chr6D.hmp.txt");
        BufferedWriter out18 = new BufferedWriter(new FileWriter(writename18));

        File writename19 = new File("/data1/home/yafei/all_vcf/chr7A.hmp.txt");
        BufferedWriter out19 = new BufferedWriter(new FileWriter(writename19));
        File writename20 = new File("/data1/home/yafei/all_vcf/chr7B.hmp.txt");
        BufferedWriter out20 = new BufferedWriter(new FileWriter(writename20));
        File writename21 = new File("/data1/home/yafei/all_vcf/chr7D.hmp.txt");
        BufferedWriter out21 = new BufferedWriter(new FileWriter(writename21));

        BufferedReader reader = null;
        reader = new BufferedReader(new FileReader(new File("/data1/home/yafei/all_vcf/all_0.05_0.05.hmp.txt")));
        String line;
        for (int i = 0; i < 1; i++) {
            line = reader.readLine();
            out1.write(line+"\n");
            out2.write(line+"\n");
            out3.write(line+"\n");
            out4.write(line+"\n");
            out5.write(line+"\n");
            out6.write(line+"\n");
            out7.write(line+"\n");
            out8.write(line+"\n");
            out9.write(line+"\n");
            out10.write(line+"\n");
            out11.write(line+"\n");
            out12.write(line+"\n");
            out13.write(line+"\n");
            out14.write(line+"\n");
            out15.write(line+"\n");
            out16.write(line+"\n");
            out17.write(line+"\n");
            out18.write(line+"\n");
            out19.write(line+"\n");
            out20.write(line+"\n");
            out21.write(line+"\n");
        }
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            String chr = list.get(2);
            switch(chr)
            {
                case "1" :
                    out1.write(line+"\n");
                    break;
                case "2" :
                    out2.write(line+"\n");
                    break;
                case "3" :
                    out3.write(line+"\n");
                    break;
                case "4" :
                    out4.write(line+"\n");
                    break;
                case "5" :
                    out5.write(line+"\n");
                    break;
                case "6" :
                    out6.write(line+"\n");
                    break;
                case "7" :
                    out7.write(line+"\n");
                    break;
                case "8" :
                    out8.write(line+"\n");
                    break;
                case "9" :
                    out9.write(line+"\n");
                    break;
                case "10" :
                    out10.write(line+"\n");
                    break;
                case "11" :
                    out11.write(line+"\n");
                    break;
                case "12" :
                    out12.write(line+"\n");
                    break;
                case "13" :
                    out13.write(line+"\n");
                    break;
                case "14" :
                    out14.write(line+"\n");
                    break;
                case "15" :
                    out15.write(line+"\n");
                    break;
                case "16" :
                    out16.write(line+"\n");
                    break;
                case "17" :
                    out17.write(line+"\n");
                    break;
                case "18" :
                    out18.write(line+"\n");
                    break;
                case "19" :
                    out19.write(line+"\n");
                    break;
                case "20" :
                    out20.write(line+"\n");
                    break;
                case "21" :
                    out21.write(line+"\n");
                    break;
            }
        }
        reader.close();
        out1.close();
        out2.close();
        out3.close();
        out4.close();
        out5.close();
        out6.close();
        out7.close();
        out8.close();
        out9.close();
        out10.close();
        out11.close();
        out12.close();
        out13.close();
        out14.close();
        out15.close();
        out16.close();
        out17.close();
        out18.close();
        out19.close();
        out20.close();
        out21.close();
    }
}
