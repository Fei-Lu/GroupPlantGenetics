/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;

/**
 *
 * @author guoyafei
 */
//import utils.IOUtils;
import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-02-17 9:09 PM
 * 提取VMapII里面对应染色体的SNP行
 */
public class GetRefVmap {
    public static Set<String> CreateNumber(String Path){
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(new File(Path)));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        Set<String> numbers = new HashSet<String>();
        String line = null;
        while (true) {
            try {
                if (!((line = reader.readLine()) != null)) {
                    break;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            List<String> list = Arrays.asList(line.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            numbers.add(sb.toString());
        }
        try {
            reader.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return numbers;
    }
    public static void main(String[] args) throws Exception {
        File writename = new File("/data1/home/yafei/test/chr23D_Ref.vcf"); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        BufferedReader reader;
        String filename = "/data1/publicData/wheat/genotype/VMapII/VMap2.1/chr023_vmap2.1.vcf.gz";
        reader = IOUtils.getTextGzipReader(filename);

        Set<String> numbers = new HashSet<String>();
        numbers = CreateNumber("/data2/yafei/Hapscanner/chr023_pos.txt");
        String line;
        for (int i = 0; i < 21; i++) {
            line = reader.readLine();
            out.write(line+"\n");
        }
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            if(numbers.contains(sb.toString())){

                out.write(line+"\n");
            }
        }
        reader.close();
//
//        BufferedReader reader1;
//        String filename1 = "/data1/publicData/wheat/genotype/VMapII/VMap2.1/chr042_vmap2.1.vcf.gz";
//        reader1 = IOUtils.getTextGzipReader(filename1);
//        Set<String> numbers1 = new HashSet<String>();
//        numbers1 = CreateNumber("/data2/yafei/Hapscanner/chr042_pos.txt");
//        String line1;
//        for (int i = 0; i < 21; i++) {
//            line1 = reader1.readLine();
//        }
//        while ((line1 = reader1.readLine()) != null) {
//            List<String> list = Arrays.asList(line1.split("\t"));
//            StringBuilder sb = new StringBuilder();
//            sb.append(list.get(0)+"_"+list.get(1));
//            if(numbers1.contains(sb.toString())){
//                int a = Integer.parseInt(list.get(1));
//                list.set(1, String.valueOf(a+453812268));
//                list.set(0, "21");
//                String change = String.join("\t", list);
//                out.write(change+"\n");
//            }
//        }
//        reader1.close();
        out.close();
//
        File writename1 = new File("/data1/home/yafei/test/chr21B_Ref.vcf"); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out1 = new BufferedWriter(new FileWriter(writename1));
        BufferedReader reader1;
        String filename1 = "/data1/publicData/wheat/genotype/VMapII/VMap2.1/chr021_vmap2.1.vcf.gz";
        reader1 = IOUtils.getTextGzipReader(filename1);

        Set<String> numbers1 = new HashSet<String>();
        numbers1 = CreateNumber("/data2/yafei/Hapscanner/chr021_pos.txt");
        String line1;
        for (int i = 0; i < 21; i++) {
            line1 = reader1.readLine();
            out1.write(line1+"\n");
        }
        while ((line1 = reader1.readLine()) != null) {
            List<String> list = Arrays.asList(line1.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            if(numbers1.contains(sb.toString())){
                out1.write(line1+"\n");
            }
        }
        reader1.close();
        out1.close();

        File writename2 = new File("/data1/home/yafei/test/chr31A_Ref.vcf"); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out2 = new BufferedWriter(new FileWriter(writename2));
        BufferedReader reader2;
        String filename2 = "/data1/publicData/wheat/genotype/VMapII/VMap2.1/chr031_vmap2.1.vcf.gz";
        reader2 = IOUtils.getTextGzipReader(filename2);

        Set<String> numbers2 = new HashSet<String>();
        numbers2 = CreateNumber("/data2/yafei/Hapscanner/chr031_pos.txt");
        String line2;
        for (int i = 0; i < 21; i++) {
            line2 = reader2.readLine();
            out2.write(line2+"\n");
        }
        while ((line2 = reader2.readLine()) != null) {
            List<String> list = Arrays.asList(line2.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            if(numbers2.contains(sb.toString())){
                out2.write(line2+"\n");
            }
        }
        reader2.close();
        out2.close();
    }
}