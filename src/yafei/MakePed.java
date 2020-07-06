package yafei;
import java.io.*;
import java.util.Arrays;
import java.util.List;
/**
 * @author Yafei Guo
 * @create 2020-03-31 9:47 AM
 */
public class MakePed {
    public static void main(String[] args) throws IOException {
        File writename = new File("/data1/home/yafei/all_vcf/result_H/struc_als/GEMMA/new_plink.ped");
        // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        String filename = "/data1/home/yafei/all_vcf/result_H/struc_als/GEMMA/output.ped";
        BufferedReader reader_ped = new BufferedReader(new FileReader(new File(filename)));
        String filename1 = "/data1/home/yafei/all_vcf/result_H/struc_als/GEMMA/traitAll.txt";
        BufferedReader reader_trait = new BufferedReader(new FileReader(new File(filename1)));
        String ped;
        String trait;
        reader_trait.readLine();
        while((ped = reader_ped.readLine()) != null && (trait = reader_trait.readLine()) != null){
            List<String> list_ped = Arrays.asList(ped.split(" "));
            List<String> list_trait = Arrays.asList(trait.split("\t"));
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < 5; i++) {
                sb.append(list_ped.get(i)+" ");
            }
            sb.append(list_trait.get(3)+" ");
            for (int i = 6; i < list_ped.size(); i++) {
                sb.append(list_ped.get(i)+" ");
            }
            out.write(sb.toString().trim()+"\n");
        }
        reader_ped.close();
        reader_trait.close();
        out.close();
    }
}
