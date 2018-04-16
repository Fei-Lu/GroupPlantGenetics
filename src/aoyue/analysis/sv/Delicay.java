/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import com.itextpdf.text.pdf.parser.Path;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class Delicay {
    public Delicay(){
        //this.mkDir(); 
        this.mkHapPosAllele();
        //创建文件并写入字符串
        this.mkfile1();
        
         
    }
    /*if (!file.exists()) {   //文件不存在则创建文件，先创建目录  
            File dir = new File(file.getParent());  
            dir.mkdirs();  
            file.createNewFile();  
        }  
        FileOutputStream outStream = new FileOutputStream(file);    //文件输出流用于将数据写入文件  
        outStream.write(sourceByte);  
        outStream.close();  //关闭文件输出流  
    } catch (Exception e) {  
        e.printStackTrace();  
    } 
*/
    public void mkfile1(){
        String input = "this is the first time that I have write a txt using java !";
        byte[] sourceByte = input.getBytes();
        if (null != sourceByte){
            try {
                File f = new File ("/Users/Aoyue/Documents/afile/a.txt");
                if (!f.exists()){
                    File dir = new File (f.getParent());
                    dir.mkdirs();
                    f.createNewFile();       
                }
                FileOutputStream outStream = new FileOutputStream(f);
                outStream.write(sourceByte);
                outStream.close();              
            }
            catch (Exception e){
                e.printStackTrace();        
            }
            
            
        }
        
 
        
        
    }
    public void mkHapPosAllele () {
        String infileDirS = "/Users/Aoyue/Desktop/";
        String outfileDirS = "/Users/Aoyue/Desktop/out/"; //该路径必须提前设置好
        //public String[] list() Returns an array of strings naming the files and directories in the directory denoted by this abstract pathname.
        //返回由此抽象路径名所表示的目录中的文件和目录的名称所组成字符串数组。
        //public File[] listFiles() Returns an array of abstract pathnames denoting the files in the directory denoted by this abstract pathname.
        //返回一个抽象路径名数组，这些路径名表示此抽象路径名所表示目录中的文件。
        File fsall = new File(infileDirS);
        File[] fs = fsall.listFiles();
        //String [] fs1 = new File (infileDirS).list(); //String类型
        fs = IOUtils.listFilesEndsWith(fs, ".vcf.gz");
        List<File> fsList = Arrays.asList(fs); //public static <T> List<T> asList(T... a) {} 把以vcf.gz结尾的文件放入fsList中
        fsList.parallelStream().forEach(f -> { //forEach流是Stream接口里的方法，对单个文件vcf.gz进行操作。
            BufferedReader br = IOUtils.getTextGzipReader(f.getAbsolutePath()); //把vcf.gz解压，将文件读进缓冲字符流中去
            String outfileS = new File(outfileDirS, f.getName().replaceFirst(".vcf.gz", ".posAllele.txt.gz")).getAbsolutePath();//输出文件名设置
            BufferedWriter bw = IOUtils.getTextGzipWriter(outfileS); //将文件写出来
            String temp = null;
            int cnt = 0;
            try {
                StringBuilder sb = new StringBuilder(); //new一个新的StringBuilder开始写文件（一个可变的字符序列）
                sb.append("Chr\tPos\tRef\tAlt"); //开始建立表头 chr Pos Alt
                bw.write(sb.toString()); //Writer类 是抽象类，用于写出字符流。
                bw.newLine();
                List<String> l = null;
                while ((temp = br.readLine()) != null) { //把当vcf.gz按行读不为空时，
                    if (temp.startsWith("#")) continue; //public boolean startsWith(String prefix) 以#开头时，继续读下一行
                    temp = temp.substring(0, 40); // 取每行的前0-40 字符串
                    l = PStringUtils.fastSplit(temp); //l最开始为一个空的list, 后将temp以 \t 分割后，放入l中。
                    sb = new StringBuilder(l.get(0)); //取集合List的第0列，即染色体号。上文已讲StringBuilder是一个可变字符序列
                    sb.append("\t").append(l.get(1)).append("\t").append(l.get(3)).append("\t").append(l.get(4)); // 将文件第0 1 3 4列的信息提取出来
                    bw.write(sb.toString()); 
                    bw.newLine();
                    if (cnt%2 == 0) System.out.println("Output " + String.valueOf(cnt) + " SNPs");
                    cnt++;
                }
                bw.flush(); //将字符输出流冲出来
                bw.close(); //将字符输出流关闭
                br.close(); //将字符输入流关闭
                System.out.println(String.valueOf(cnt) + " SNPs output from " + f.getAbsolutePath());
            }
            catch (Exception e) {
                System.out.println(f.getAbsolutePath()); //如果文件有异常，就将文件路径打印出来
                System.out.println(temp); 
                e.printStackTrace();
                System.exit(1);
            }
            
        });
    }
    
    public static File[] listFilesEndsWith (File[] all, String endStr) {
        ArrayList<File> al = new ArrayList(); //建立一个新的集合类 al， 并且规定该类是 File类型.
        //将以endStr结尾的文件添加到集合al 中
        for (int i = 0; i < all.length; i++) { 
            if (all[i].getName().endsWith(endStr)) al.add(all[i]);
        }
        return al.toArray(new File[al.size()]); //size()返回集合元素的数量; new File[3]; toArray返回T[] a（即返回集合中的元素）
    }
    
    
    public void mkDir (){
        String dir = "/Users/Aoyue/Documents/afile";
        File f = new File (dir); 
        if (!f.exists()) {  
            f.mkdir();  
        }  
        System.out.println("Perform the end "+ dir);  
    }  
    public static void main (String[] args){        
        new Delicay (); 
        System.out.println("I can do it !");
        
        
    }
          
    } 
     

