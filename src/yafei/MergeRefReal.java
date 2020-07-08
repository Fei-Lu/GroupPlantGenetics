package yafei;

import java.io.*;
import java.util.Arrays;
import java.util.List;

/**
 * @author Yafei Guo
 * @create 2020-03-11 9:42 PM
 */
public class MergeRefReal {
    public static void main(String[] args) throws Exception {
        File writename = new File("/data1/home/yafei/all_test/chr1A_all.vcf");
        // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        out.write("##fileformat=VCFv4.1\n" +
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and al\n" +
                "##FORMAT=<ID=GL,Number=.,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1,\n" +
                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed st\n" +
                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or\n" +
                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
                "##ALT=<ID=DEL,Description=\"Deletion\">\n" +
                "##ALT=<ID=INS,Description=\"Insertion\">\n");

        String filename1 = "/data1/home/yafei/chr1A_Ref.vcf";
        BufferedReader reader1 = new BufferedReader(new FileReader(new File(filename1)));
        String filename2 = "/data1/home/yafei/chr1A_real.vcf";
        BufferedReader reader2 = new BufferedReader(new FileReader(new File(filename2)));
        String line1;
        String line2;
        for (int i = 0; i < 20; i++) {
            reader1.readLine();
        }
        for (int i = 0; i < 13; i++) {
            reader2.readLine();
        }
        while((line1 = reader1.readLine()) != null && (line2 = reader2.readLine()) != null){
            List<String> list = Arrays.asList(line2.split("\t"));
            StringBuilder sb = new StringBuilder();
            for (int i = 1; i < list.size(); i++) {
                sb.append(list.get(i)+"\t");
            }
            out.write(line1+"\t"+sb.toString().trim()+"\n");
        }
        reader1.close();
        reader2.close();
        out.close();

    }
}
