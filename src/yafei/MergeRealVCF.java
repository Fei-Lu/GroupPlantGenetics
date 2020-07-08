package yafei;

//import utils.IOUtils;

import java.io.*;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;

/**
 * @author Yafei Guo
 * @create 2020-02-28 10:06 PM
 */
public class MergeRealVCF {
    public static void main(String[] args) throws IOException {
        getSample("/data1/home/yafei/chr1A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr001.vcf.gz", "1","/data2/yafei/Hapscanner/out/VCF/1832_chr002.vcf.gz",471304005);
        getSample("/data1/home/yafei/chr1B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr003.vcf.gz", "2","/data2/yafei/Hapscanner/out/VCF/1832_chr004.vcf.gz",438720154);
        getSample("/data1/home/yafei/chr1D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr005.vcf.gz", "3","/data2/yafei/Hapscanner/out/VCF/1832_chr006.vcf.gz",452179604);
        getSample("/data1/home/yafei/chr2A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr007.vcf.gz", "4","/data2/yafei/Hapscanner/out/VCF/1832_chr008.vcf.gz",462376173);
        getSample("/data1/home/yafei/chr2B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr009.vcf.gz", "5","/data2/yafei/Hapscanner/out/VCF/1832_chr010.vcf.gz",453218924);
        getSample("/data1/home/yafei/chr2D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr011.vcf.gz", "6","/data2/yafei/Hapscanner/out/VCF/1832_chr012.vcf.gz",462216879);
        getSample("/data1/home/yafei/chr3A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr013.vcf.gz", "7","/data2/yafei/Hapscanner/out/VCF/1832_chr014.vcf.gz",454103970);
        getSample("/data1/home/yafei/chr3B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr015.vcf.gz", "8","/data2/yafei/Hapscanner/out/VCF/1832_chr016.vcf.gz",448155269);
        getSample("/data1/home/yafei/chr3D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr017.vcf.gz", "9","/data2/yafei/Hapscanner/out/VCF/1832_chr018.vcf.gz",476235359);
        getSample("/data1/home/yafei/chr4A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr019.vcf.gz", "10","/data2/yafei/Hapscanner/out/VCF/1832_chr020.vcf.gz",452555092);
        getSample("/data1/home/yafei/chr4B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr021.vcf.gz", "11","/data2/yafei/Hapscanner/out/VCF/1832_chr022.vcf.gz",451014251);
        getSample("/data1/home/yafei/chr4D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr023.vcf.gz", "12","/data2/yafei/Hapscanner/out/VCF/1832_chr024.vcf.gz",451004620);
        getSample("/data1/home/yafei/chr5A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr025.vcf.gz", "13","/data2/yafei/Hapscanner/out/VCF/1832_chr026.vcf.gz",453230519);
        getSample("/data1/home/yafei/chr5B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr027.vcf.gz", "14","/data2/yafei/Hapscanner/out/VCF/1832_chr028.vcf.gz",451372872);
        getSample("/data1/home/yafei/chr5D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr029.vcf.gz", "15","/data2/yafei/Hapscanner/out/VCF/1832_chr030.vcf.gz",451901030);
        getSample("/data1/home/yafei/chr6A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr031.vcf.gz", "16","/data2/yafei/Hapscanner/out/VCF/1832_chr032.vcf.gz",452440856);
        getSample("/data1/home/yafei/chr6B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr033.vcf.gz", "17","/data2/yafei/Hapscanner/out/VCF/1832_chr034.vcf.gz",452077197);
        getSample("/data1/home/yafei/chr6D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr035.vcf.gz", "18","/data2/yafei/Hapscanner/out/VCF/1832_chr036.vcf.gz",450509124);
        getSample("/data1/home/yafei/chr7A_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr037.vcf.gz", "19","/data2/yafei/Hapscanner/out/VCF/1832_chr038.vcf.gz",450046986);
        getSample("/data1/home/yafei/chr7B_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr039.vcf.gz", "20","/data2/yafei/Hapscanner/out/VCF/1832_chr040.vcf.gz",453822637);
        getSample("/data1/home/yafei/chr7D_real.vcf","/data2/yafei/Hapscanner/out/VCF/1832_chr041.vcf.gz", "21","/data2/yafei/Hapscanner/out/VCF/1832_chr042.vcf.gz",453812268);
    }
    public static void getSample(String outFile, String inputFile1, String chr, String inputFile2, int num) throws IOException {
        File writename = new File(outFile); // 相对路径，如果没有则要建立一个新的output.txt文件
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        BufferedReader reader;
        String filename = inputFile1;
        reader = IOUtils.getTextGzipReader(filename);
        String line;
        for (int i = 0; i < 14; i++) {
            line = reader.readLine();
            out.write(line+"\n");
        }
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            list.set(0, chr);
            String change = String.join("\t", list);
            out.write(change+"\n");

        }
        reader.close();
        BufferedReader reader1;
        String filename1 = inputFile2;
        reader1 = IOUtils.getTextGzipReader(filename1);
        String line1;
        for (int i = 0; i < 14; i++) {
            line1 = reader1.readLine();
        }
        while ((line1 = reader1.readLine()) != null) {
            List<String> list = Arrays.asList(line1.split("\t"));
            int a = Integer.parseInt(list.get(1));
            list.set(1, String.valueOf(a+num));
            list.set(0, chr);
            String change = String.join("\t", list);
            out.write(change+"\n");
        }
        reader1.close();
        out.close();
    }
}
