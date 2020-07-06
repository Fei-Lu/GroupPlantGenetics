package yafei;

//import utils.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;
/**
 * @author Yafei Guo
 * @create 2020-02-18 6:00 PM
 */
public class GetVCFSample {
    public static void main(String[] args) throws Exception {
        File writename = new File("/data1/home/yafei/chr01_real.vcf"); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        BufferedReader reader;
        String filename = "/data2/yafei/Hapscanner/out/VCF/1832_chr001.vcf.gz";
        reader = IOUtils.getTextGzipReader(filename);
        String line;
        while ((line = reader.readLine()) != null) {
            out.write(line+"\n");
        }
        reader.close();

        BufferedReader reader1;
        String filename1 = "/data2/yafei/Hapscanner/out/VCF/1832_chr002.vcf.gz";
        reader1 = IOUtils.getTextGzipReader(filename1);
        String line1;
        for (int i = 0; i < 14; i++) {
            line1 = reader1.readLine();
        }
        while ((line1 = reader1.readLine()) != null) {
            List<String> list = Arrays.asList(line1.split("\t"));
            int a = Integer.parseInt(list.get(1));
            list.set(1, String.valueOf(a+471304005));
            list.set(0, "1");
            String change = String.join("\t", list);
            out.write(change+"\n");
        }
        reader1.close();
        out.close();
    }
}
