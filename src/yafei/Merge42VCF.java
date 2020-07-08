package yafei;

//import utils.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-02-19 7:17 PM
 */
public class Merge42VCF {
    public static void main(String[] args) throws Exception {
        File writename = new File("/data1/home/yafei/all_vcf/result_H/all.hmp.txt");
        // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
//        out.write("##fileformat=VCFv4.1\n" +
//                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n" +
//                "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the reference and al\n" +
//                "##FORMAT=<ID=GL,Number=.,Type=Integer,Description=\"Genotype likelihoods for 0/0, 0/1, 1/1,\n" +
//                "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n" +
//                "##INFO=<ID=NZ,Number=1,Type=Integer,Description=\"Number of taxa with called genotypes\">\n" +
//                "##INFO=<ID=AD,Number=.,Type=Integer,Description=\"Total allelelic depths in order listed st\n" +
//                "##INFO=<ID=AC,Number=.,Type=Integer,Description=\"Numbers of ALT alleles in order listed\">\n" +
//                "##INFO=<ID=GN,Number=.,Type=Integer,Description=\"Number of taxa with genotypes AA,AB,BB or\n" +
//                "##INFO=<ID=HT,Number=1,Type=Integer,Description=\"Number of heterozygotes\">\n" +
//                "##INFO=<ID=MAF,Number=1,Type=Float,Description=\"Minor allele frequency\">\n" +
//                "##ALT=<ID=DEL,Description=\"Deletion\">\n" +
//                "##ALT=<ID=INS,Description=\"Insertion\">\n");
        BufferedReader reader;
        File folder = new File("/data1/home/yafei/all_vcf/result_H");
        File[] files = folder.listFiles(new FileFilter() {
            @Override
            public boolean accept(File file) {
                if(file.getName().endsWith(".gz") && file.isFile()){
                    return true;
                }
                return false;   //否则过滤掉
            }
        });
        List<File> sortfiles = Arrays.asList(files);
        Collections.sort(sortfiles, new Comparator<File>(){
            public int compare(File o1, File o2) {
                return o1.getName().compareTo(o2.getName());
            }
        });
        boolean flag = true;
        for(File file : sortfiles) {
            if(flag) {
                String filename = file.getAbsolutePath();
                reader = IOUtils.getTextGzipReader(filename);
//                reader = new BufferedReader(new FileReader(new File(file.getAbsolutePath())));
//                System.out.println(file.getAbsolutePath());
                String line;
                while ((line = reader.readLine()) != null) {
//                    System.out.println(line+"\n");
                    out.write( line+ "\n");
                }
                reader.close();
                flag = false;
            }else{
                String filename = file.getAbsolutePath();
                reader = IOUtils.getTextGzipReader(filename);
//                reader = new BufferedReader(new FileReader(new File(file.getAbsolutePath())));
//                System.out.println(file.getAbsolutePath());
                String line;
                reader.readLine();
                while ((line = reader.readLine()) != null) {
//                    System.out.println(line+"\n");
                    out.write( line + "\n");
                }
                reader.close();
            }
        }
        out.close();
    }
}
