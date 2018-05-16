/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;

/**
 *
 * @author yaozhou
 */
public class aoyueshared {
    public void test() {
        YaoIOUtils a  = new YaoIOUtils();
        
        File b = new File("/Users/Aoyue/Documents/NSF");
        File[] list = a.listRecursiveFiles(b);
        for (int i = 0; i < list.length; i++) {
//            System.out.println(list[i].getName());
            System.out.println(list[i].getAbsolutePath());
//            System.out.println(list[i].length());
//            System.out.println(list[i].getParent());
        }
        
        
    } 
    public void test2 () {
        String path = "/Volumes/Lulab3T_14/20171120_copy";
        File test = new File (path);
        File[] fs = YaoIOUtils.listRecursiveFiles(test);
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
        File[] subFs = YaoIOUtils.listFilesEndsWith(fs, "clean.fq.gz");
        String outfileS = "/Users/Aoyue/Documents/Datalist.txt";
        String[] header = {"FileName", "Path", "FileSize"};
        try {
            BufferedWriter bw = YaoIOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < header.length; i++) {
                sb.append(header[i]).append("\t");
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 0; i < subFs.length; i++) {
                sb = new StringBuilder();
                
                double qq = (double)subFs[i].length()/1024/1024/1024;
                DecimalFormat df = new DecimalFormat ("0.00");

                sb.append(subFs[i].getName()).append("\t").append(subFs[i].getAbsolutePath()).append("\t").append(df.format(qq)).append("G").append("\t").append("350bp");


//                sb.append(subFs[i].getName()).append("\t").append(subFs[i].getAbsolutePath()).append("\t").append((double)subFs[i].length()/1024/1024/1024).append("\t").append("350bp");
                
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println(subFs.length);
    }
}
