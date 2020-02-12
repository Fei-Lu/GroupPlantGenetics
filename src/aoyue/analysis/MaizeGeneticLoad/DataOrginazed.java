/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import pgl.infra.table.RowTable;
import pgl.infra.table.TableInterface;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.atomic.AtomicInteger;

/**
 *
 * @author Aoyue
 */
public class DataOrginazed {
    public DataOrginazed(){
        this.listAllFiles();
        this.listSpecificalFiles();
        this.md5();
        this.checkMd5();
        //this.test();
        this.coverage();
        this.sample();
        //this.insertTxt();
        //this.findlost();
        //this.findlostnull();
       
    }
    public void findlost() {
        File[] fs = IOUtils.listRecursiveFiles(new File("/Users/Aoyue/Desktop/output_data_stats_434SM500bp"));
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".html");
        String[] filename = new String[subFs.length];
        Set<String> html = new HashSet();        
        for (int i=0; i < subFs.length; i++){
          filename[i] = subFs[i].getName().replaceFirst(".bc.html", "");              
          html.add(filename[i]);              
        }           
        try{
            BufferedReader br = IOUtils.getTextReader("/Users/Aoyue/Desktop/434samplelist.txt");
            String temp = null;
            while((temp = br.readLine()) != null){
                String tem = temp;
               if(html.add(tem)) {
                   System.out.println(tem);
                }              
            }               
        }
        catch (IOException e){
            e.printStackTrace();
            System.exit(1);
        }            
    }
              
    public void findlostnull(){
        String infile = "/Users/Aoyue/Desktop/samtools_stats_upload/434samplelist.txt";
        String infileS = "/Users/Aoyue/Desktop/output_data_stats_434SM500bp";
        //读进表格的文件，必须要有表头。故这里自己要手动加入表头。
        RowTable <String> t = new RowTable<>(infile);
        List<String> nameList = new ArrayList <>();
        for (int i = 0; i < t.getRowNumber(); i++){
            nameList.add(t.getCell(i, 0));
        }
        File[] fs = new File (infileS).listFiles();
        File[] subFs = IOUtils.listFilesEndsWith(fs, ".bc.html");
        HashSet<String> sampleSet = new HashSet();
        for (int i = 0; i < subFs.length; i++){
            if (subFs[i].isHidden()) continue;
            sampleSet.add(subFs[i].getName().split(".bc.html")[0]);                        
        }
        for (int i = 0; i < nameList.size(); i++){
            if (!(sampleSet.contains(nameList.get(i)))){
                System.out.println(nameList.get(i));
                int a =3;
            }
        }                
    }
    
    public void insertTxt(){
        String infileS = "/Users/Aoyue/Desktop/K16BJS0001.yao.list";
        String outfileS = "/Users/Aoyue/Desktop/mvK16BJS0001.yao.listt";
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String temp = null;
            while ((temp = br.readLine()) != null){
                String tem = temp; //
                StringBuilder sb =new StringBuilder();
                sb.append("mv").append(" ").append(tem).append(" ").append("../plot1/");
                bw.write(sb.toString() + "\n");
                //System.out.println(sb.toString());
            }
            bw.flush();
            bw.close();
            //br.close(); 
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
  public void sample(){
        String infileS = "/Users/Aoyue/Downloads/All_Merged_1.25M_MAF0.05.sort.hmp";
        String outfileS = "/Users/Aoyue/Downloads/All_Merged_1.25M_MAF0.05.sort.sample100000.hmp";
        int length = 10000;
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for (int i = 0; i < length; i++){
                bw.write(br.readLine());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            br.close();          
        }
        catch(Exception e){
            e.printStackTrace();        
        }        
    }  
    
    public void coverage () {
        String infileDirS = "/Volumes/Lulab3T_14/20171120CAAS/P101SC17081532_01zhangshaojing/data_release/cleandata";
        String outfileS = "/Users/Aoyue/Documents/Data/project/maize2k/coverage/111SM_14_coverage.txt";
        File[] fs = new File(infileDirS).listFiles(); 
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        int genomeSize = 2135098301;
        //
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
        //
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
        //String path = "/Volumes/Lulab3T_14/20171120_copy/P101SC17081532_01zhangshaojing/data_release/cleandata";
        //File test = new File (path);
        File[] fs = IOUtils.listRecursiveFiles(new File("/Users/Aoyue/Documents/Data/pipeline/hapScanner/hapPosAllele/"));
        File[] subFs = IOUtils.listFilesEndsWith(fs, "posAllele.txt.gz"); //列出以posAllele.txt.gz结尾的文件
        String outfileS = "/Users/Aoyue/Documents/md5.txt";
        try {
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < subFs.length; i++) {
                sb = new StringBuilder();
//                sb.append("md5").append(" ").append(subFs[i].getAbsolutePath()).append("\t");
                sb.append("md5").append(" ").append(subFs[i].getName());        
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
        //ori 文件 54de998c7883a4592a1683bec2590d64  K16HL0119_1_clean.fq.gz
        //des文件 md5 (K16HL0133_2_clean.fq.gz) = 52c7a9d501d5fd7b3d0deaf3aae1c715
        String des = "/Users/Aoyue/Documents/111md5.txt"; //mac生成的md5文件
        String ori = "/Users/Aoyue/Documents/111originmd5.txt"; //原始md5
        HashMap<String, String> fmd5Map = new HashMap<>(); //建立一个键值对应的hashmap,此时hashmap为空。下文会把原始的ori文件放入hashmap中去
        TableInterface oT = new RowTable(ori, "  "); //分隔符是2个空格，将ori文件读进表格
        TableInterface dT = new RowTable(des, " "); //分隔符是1个空格，将des文件读进表格
        for (int i = 0; i < oT.getRowNumber(); i++) { //按行遍历
            fmd5Map.put(oT.getCellAsString(i, 1), oT.getCellAsString(i, 0)); //Return a string value of a cell. 将ori文件中的第2列放进hashmap中的key，
            //将第1列放进hashmap中的value
        }
        for (int i = 0; i < dT.getRowNumber(); i++) {
            //将des文件中的第2列样本名提取出来，并将括号去掉，用replaceFirst方法
            String key = dT.getCellAsString(i, 1).replaceFirst("\\(", "").replaceFirst("\\)", ""); //此处的key为 des文件中的样本名
            String value = fmd5Map.get(key); //get(key)返回fmd5Map中的value值。
            if (value == null) { //如果value为空，则为真。
                System.out.println(key+"\tdoesn't exist");
                continue;
            }
            if (value.equals(dT.getCellAsString(i, 3))) continue;
            System.out.println(key + "\t is incorrect");
        }
    }
    
    
    
    
}
