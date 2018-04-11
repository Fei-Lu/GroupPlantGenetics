/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import format.table.RowTable;
import format.table.TableInterface;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;
import utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class DataOrginazed {
    public DataOrginazed(){
//        this.listAllFiles();
//        this.listSpecificalFiles();
//        this.md5();
//        this.checkMd5();
        //this.test();
        this.covergage();
       
        
        
    }
    
    public void covergage () {
        String infileDirS = "/Volumes/Lulab3T_14/20171120CAAS/P101SC17081532_01zhangshaojing/data_release/cleandata";
        String outfileS = "/Users/Aoyue/Documents/Data/project/maize2k/coverage/111SM_14_coverage.txt";
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        int genomeSize = 2135098301;
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        ConcurrentHashMap<String, Integer> fileSizeMap = new ConcurrentHashMap();
        AtomicInteger acnt = new AtomicInteger();
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1_clean.fq.gz").getAbsolutePath();
            try {
                
                BufferedReader br = IOUtils.getTextGzipReader(infile1);
                String temp = br.readLine();
                temp = br.readLine();
                int len = temp.length();
                acnt.set(len);
                br.close();
                br = IOUtils.getTextGzipReader(infile1);
                int cnt = 0;
                while ((temp = br.readLine()) != null) {                  
                    cnt++;
                    br.readLine();br.readLine();br.readLine();                 
                }
                
                fileSizeMap.put(name, cnt);
                
                System.out.println(name+" has "+String.valueOf(cnt)+ " reads. Read length: " +String.valueOf(len));
            }
            catch (Exception e) {
                e.printStackTrace();
            }
        });
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Sample\tReadsNum\tReadsLength\tCoverage");
            bw.newLine();
            String[] names = nameSet.toArray(new String[nameSet.size()]);
            Arrays.sort(names);
            for (int i = 0; i < names.length; i++) {
                StringBuilder sb = new StringBuilder(names[i]);
                int readNum = fileSizeMap.get(names[i]);
                int readLength = acnt.intValue();
                sb.append("\t").append(readNum).append("\t").append(readLength).append("\t").append((double)readNum*readLength*2/genomeSize);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    public void test(){
        System.out.println("hello");
        System.out.println("try");
        System.out.println("haah");
        System.out.println("today is the last day in Mar");
    }
    private void listAllFiles() {
        IOUtils a = new IOUtils();
        File b = new File ("/Users/Aoyue/Documents/集群使用前信息整理");
        File[] list = a.listRecursiveFiles(b);
        for (int i = 0; i < list.length; i++) {
//            System.out.println(list[i].getName());
            System.out.println(list[i].getAbsolutePath());
//            System.out.println(list[i].length());
//            System.out.println(list[i].getParent());
        }       
    }      

    private void listSpecificalFiles() {
        String path = "/Volumes/Seagate Backup Plus Drive/BJ170733-01遗传所96个小麦DNA建库PCR-free测序过滤任务单/upload/Cleandata";
        File test = new File (path);
        File[] fs = IOUtils.listRecursiveFiles(test);
//        File[] fs = IOUtils.listRecursiveFiles(new File(path));
        File[] subFs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        String outfileS = "/Users/Aoyue/Documents/Datalist.txt";
        String[] header = {"FileName", "Path", "FileSize","Library","Company"};
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
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

                sb.append(subFs[i].getName()).append("\t").append(subFs[i].getAbsolutePath()).append("\t")
                  .append(df.format(qq)).append("G").append("\t").append("350bp").append("\t").append("安诺优达基因");

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
    
    private void md5() {
        String path = "/Volumes/Lulab3T_14/20171120_copy/P101SC17081532_01zhangshaojing/data_release/cleandata";
        File test = new File (path);
        File[] fs = IOUtils.listRecursiveFiles(test);
        File[] subFs = IOUtils.listFilesEndsWith(fs, "clean.fq.gz");
        String outfileS = "/Users/Aoyue/Documents/md5.txt";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < subFs.length; i++) {
                sb = new StringBuilder();
//                sb.append("md5").append(" ").append(subFs[i].getAbsolutePath()).append("\t");
                sb.append("md5").append(" ").append(subFs[i].getName()).append("\t");        
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

    private void checkMd5() {
        String des = "/Users/Aoyue/Documents/111md5.txt";
        String ori = "/Users/Aoyue/Documents/111originmd5.txt";
        HashMap<String, String> fmd5Map = new HashMap<>();
        TableInterface oT = new RowTable(ori, "  ");
        TableInterface dT = new RowTable(des, " ");
        for (int i = 0; i < oT.getRowNumber(); i++) {
            fmd5Map.put(oT.getCellAsString(i, 1), oT.getCellAsString(i, 0));
        }
        for (int i = 0; i < dT.getRowNumber(); i++) {
            String key = dT.getCellAsString(i, 1).replaceFirst("\\(", "").replaceFirst("\\)", "");
            String value = fmd5Map.get(key);
            if (value == null) {
                System.out.println(key+"\tdoesn't exist");
                continue;
            }
            if (value.equals(dT.getCellAsString(i, 3))) continue;
            System.out.println(key + "\t is incorrect");
        }
    }
    
    
}
