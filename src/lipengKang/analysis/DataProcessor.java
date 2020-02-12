/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;
import pgl.infra.table.RowTable;
import pgl.infra.table.TableInterface;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

/**
 *
 * @author kanglipeng
 */
/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */



/**
 *
 * @author Aoyue
 */
public class DataProcessor {
    public DataProcessor() throws IOException{
        this.sample();
       //this.mkMd51();
        //this.runtime();
        //this.mergeFile();
       //this.mkMd5();
       //this.checkMd5();
       
       
        
    }
    
    private void checkMd5() {
        //ori 文件 54de998c7883a4592a1683bec2590d64  K16HL0119_1_clean.fq.gz
        //des文件 MD5 (/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/180803_I13_V100004234_L1_WHEkapRAADT-585_1.clean.fq.gz) = 472e1633b40b0ae3866a820f5b8637a5
        String ori = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/md5.txt"; //原始md5
        String des = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/checkmd5.txt"; //mac生成的md5文件
        HashMap<String, String> fmd5Map = new HashMap<>(); //建立一个键值对应的hashmap,此时hashmap为空。下文会把原始的ori文件放入hashmap中去
        TableInterface oT = new RowTable(ori, "  "); //分隔符是2个空格，将ori文件读进表格
        TableInterface dT = new RowTable(des, " "); //分隔符是1个空格，将des文件读进表格
        try{
            BufferedReader br = IOUtils.getTextReader(ori);
            String temp;
            List<String> l = null;
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, "  ");
                fmd5Map.put(l.get(1), l.get(0));
            }
            br.close();
            int cnt = 0;
            br = IOUtils.getTextReader(des);
            while((temp = br.readLine()) != null){
                l = PStringUtils.fastSplit(temp, " ");
                String key = l.get(1).split("00010496/")[1].replaceFirst("\\)", "");
                String value = fmd5Map.get(key);
                if (value == null) { //如果value为空，则为真。
                System.out.println(key+"\tdoesn't exist");
                continue;
            }
            if (value.equals(l.get(3))) {
                cnt++;
                continue;
            }
            System.out.println(key + "\t is incorrect");
            }
            System.out.println(cnt + " key is correct");
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    
    public void mkMd5(){
        /*该方法只用将文件输入路径和输出文件名进行修改，即可使用*/
        String inputDirS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496";
        String outfileS = "/Users/Aoyue/Downloads/huadagene/WHB5EXONPEP00010496/checkmd5.txt";
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            /*1#新建存放数据的文件对象fall，将该目录下所有以fq.gz结尾的文件名列出来，并存放在fs文件数组中*/
            File fall = new File (inputDirS);
            File[] fs = IOUtils.listRecursiveFiles(fall);
            fs =  IOUtils.listFilesEndsWith(fs, ".gz");
            /*2# 根据fq.gz结尾的文件，生成md5计算的脚本，输入到StringBuilder sb中（这里用循环的方法）*/
            for (int i = 0; i < fs.length; i++) {
                sb.append("md5").append(" ").append(fs[i].getAbsoluteFile()).append("\n");
            }
            /*3# 将生成的sb脚本转换成String类型，对象名为cmd，并打印出来检查是否无误*/
            String cmd = sb.toString();
            System.out.println(cmd);
            /*4# 调用终端命令，运行md5命令，并将输出的结果写入输出文件outfileS中*/
            Runtime run = Runtime.getRuntime();
            Process p = run.exec(cmd);           
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuilder ssb = new StringBuilder();
            String line = null;
            while ((line = br.readLine()) != null) {
            ssb.append(line).append("\n");
            }
            String result = ssb.toString();
            System.out.println(result);  
            p.waitFor();
            bw.write(result);
            bw.flush();
            bw.close();            
            System.out.println("md5Check is finished at " + outfileS);
        }
        catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);       
        }       
    }
    
    
    public void mergeFile(){
        String inputDirS = "/Users/Aoyue/Documents/testCheckMd5/";
        String outfileS = "/Users/Aoyue/Documents/orimd5_merge.txt";
        File fall = new File (inputDirS);
        File[] orimd5 = IOUtils.listRecursiveFiles(fall);
        orimd5 = IOUtils.listFilesEndsWith(orimd5, ".md5");
        List<File> orimd5List = Arrays.asList(orimd5);
        try{
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            for(int i = 0; i< orimd5.length; i++){
                BufferedReader br = IOUtils.getTextReader(orimd5List.get(i).getPath());
                String temp = null;
                while((temp = br.readLine())!= null){
                    sb.append(temp).append("\n");                   
                }                
                String a = orimd5List.get(i).getName();
                System.out.println(a);
            }
            bw.write(sb.toString());
            bw.flush();
            bw.close();
            //br.close();            
        }
        catch (Exception e){
            e.printStackTrace();
            System.exit(1);
        }       
    }
    
    public void runtime() throws IOException{         
        try {
            StringBuilder sb = new StringBuilder("/sbin/md5");
            BufferedWriter bw = IOUtils.getTextWriter("/Users/Aoyue/Documents/checkMd5.txt");
            sb.append(" ").append("/Users/Aoyue/Documents/testCheckMd5/data1_00/SRR1575344_1.fq.gz").append(" ").append(">").append(" ").append("/Users/Aoyue/Documents/testCheckMd5.txt");
            String cmd = sb.toString();
            System.out.println(cmd);
            Runtime run = Runtime.getRuntime();   
            Process p = run.exec(cmd);           
            BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
            StringBuffer ssb = new StringBuffer();
            String line;
            while ((line = br.readLine()) != null) {
            ssb.append(line).append("\n");
            }
            String result = ssb.toString();
            bw.write(result);
            bw.flush();
            bw.close();
          //System.out.println(result); 
            p.waitFor();
        }
        catch (Exception e) {
        System.out.println(e.toString());
        System.exit(1);
        }       
        System.out.println("md5Check is finished");
    }
    
    public void mkMd51(){
        String infileDirS = "/Users/Aoyue/Documents/testCheckMd5";
        //String outfileS = "/Users/Aoyue/Documents/testCheckMd5.txt";
        File fsall = new File(infileDirS);
        File[] fs = IOUtils.listRecursiveFiles(fsall);
        fs = IOUtils.listFilesEndsWith(fs, "fq.gz");
        List<File> fsList = Arrays.asList(fs);
        for(int i = 0; i< fs.length; i++){
            String a = fsList.get(i).getName();
            System.out.println(a);
        }  
        fsList.parallelStream().forEach(f -> {
            StringBuilder sb = new StringBuilder("/sbin/md5");
            sb.append(" ").append(f.getPath());
            String cmd = sb.toString();
            System.out.println(cmd);       
            Runtime run = Runtime.getRuntime();        
            try {
                Process p = run.exec(cmd);
                p.waitFor();
                BufferedReader br = new BufferedReader(new InputStreamReader(p.getInputStream()));
                StringBuffer ssb = new StringBuffer();
                String line;
                while ((line = br.readLine()) != null) {
                ssb.append(line);
                }
                String result = ssb.toString();
                System.out.println(result);              
            } 
            catch (Exception e) {
            System.out.println(e.toString());
            System.exit(1);
        }            
        });        
        System.out.println("md5Check is finished");       
    }
   
    public void sample(){
        String infileS = "/Volumes/Lulab3T_14/maize_germ/001_genotype/282/282_libs_2015/uplifted_APGv4/hmp321_agpv4_chr1.vcf.gz";
        String outfileS = "/Volumes/Lulab3T_14/maize_germ/001_genotype/282/282_libs_2015/uplifted_APGv4/hmp321_agpv4_chr1_100sites.vcf";
        int length = 100;
        try{
            BufferedReader br = IOUtils.getTextGzipReader(infileS);
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
}
