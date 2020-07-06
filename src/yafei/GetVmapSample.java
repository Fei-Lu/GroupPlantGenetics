package yafei;

//import utils.IOUtils;
import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-02-17 9:53 AM
 */
public class GetVmapSample {
    public static Set<String> CreateNumber(){
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(new File("/data2/yafei/chr036_pos.txt")));
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
        File writename = new File("/data2/yafei/output"); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));

        BufferedReader reader;
        String filename = "/data2/yafei/chr036_vmap2.1.vcf.gz";
        reader = IOUtils.getTextGzipReader(filename);
        Set<String> numbers = new HashSet<String>();
        numbers = CreateNumber();
        String line;
        for (int i = 0; i < 20; i++) {
            line = reader.readLine();
        }
        for (int i = 0; i < 13; i++) {
            out.write("\n");
        }
        line = reader.readLine();
        List<String> list3 = Arrays.asList(line.split("\t"));
        StringBuilder sb3 = new StringBuilder();
        for (int i = 9; i < list3.size(); i++) {
            sb3.append(list3.get(i)+"\t");
        }
        out.write(sb3.toString().trim()+"\n");
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            StringBuilder sb = new StringBuilder();
            sb.append(list.get(0)+"_"+list.get(1));
            if(numbers.contains(sb.toString())){
                List<String> list2 = Arrays.asList(line.split("\t"));
                StringBuilder sb2 = new StringBuilder();
                for (int i = 9; i < list2.size(); i++) {
                    sb2.append(list2.get(i)+"\t");
                }
                out.write(sb2.toString().trim()+"\n");
            }
        }
        reader.close();
        out.close();
    }
}