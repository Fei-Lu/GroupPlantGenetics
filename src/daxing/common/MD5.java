package daxing.common;

import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.Predicate;

/**
 *
 * @author Daxing Xu
 */
public class MD5 {

    /**
     * 返回输入文件的MD5值
     * @param inputFile 输入文件的绝对路径
     * @return 返回输入文件inputFile的MD5值
     */
    public static String getMD5FromFile(String inputFile){
        String md5Value=null;
        MessageDigest md;
        FileInputStream fis;
        try{
            md=MessageDigest.getInstance("md5");
            fis=new FileInputStream(inputFile);
            if(fis.available()<(1024*1024*1024)){
                md.update(Files.readAllBytes(Paths.get(inputFile)));
            } else {
                byte[] b=new byte[65536];
                int count;
                while ((count=fis.read(b))!=-1){
                    md.update(b,0,count);
                }
            }
            fis.close();
            byte[] digest=md.digest();
            StringBuilder sb=new StringBuilder();
            for(int i=0;i<digest.length;i++){
                sb.append(Integer.toString((digest[i] & 0xff) + 0x100, 16).substring(1));
            }
            md5Value=sb.toString();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        return md5Value;
    }

    /**
     * 以标准MD5文件的形式返回一个目录下所有文件的MD5值。在服务器上默认使用32个核；在个人电脑上默认使用全部核
     * @param inputDir 输入目录的绝对路径
     */
    public static void getMD5FromDir(String inputDir){
        long start = System.nanoTime();
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        File[] files = IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=e->e.getName().contains(".DS_Store");
        File[] fileArray= Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        if(fileArray.length<numThreads){
            numThreads=fileArray.length;
        }
        ConcurrentHashMap<String,String> md5PathValueMap=new ConcurrentHashMap<>();
        BufferedWriter bw=null;
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(fileArray.length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(index-> {
                        String md5Value=MD5.getMD5FromFile(fileArray[index].getAbsolutePath());
                        Path AbsolutePath = fileArray[index].toPath();
                        Path relPath=AbsolutePath.subpath(Paths.get(inputDir).getNameCount(), AbsolutePath.getNameCount());
                        md5PathValueMap.put(relPath.toString(), md5Value);
                    });
        }
        try{
            bw=IOUtils.getTextWriter(inputDir+"/md5.txt");
            for(Map.Entry<String,String> entry:md5PathValueMap.entrySet()){
                //System.out.println(entry.getKey()+"  "+entry.getValue());
                bw.write(entry.getValue()+"  "+entry.getKey());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
    }


    /**
     * 判断输入文件与MD5值是否相等
     * @param inputFile 输入文件的绝对路径
     * @param hash MD5值
     * @return 如果输入文件的MD5值与给定MD5值相等就返回true
     */
    public static boolean checkMD5(String inputFile, String hash){
        String fileHash=MD5.getMD5FromFile(inputFile);
        return fileHash.equals(hash);
    }

    /**
     * 对标准MD5文件进行校验(MD5文件需在等待检测MD5值的文件目录下)。在服务器上默认使用32个核；在个人电脑上默认使用全部核
     * @param inputMD5File 输入的标准MD5文件的绝对路径
     */
    public static void checkMD5(String inputMD5File){
        String inputDir=new File(inputMD5File).getParent();
        long start = System.nanoTime();
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        List<String> md5ValuePath=new ArrayList<>();
        String line;
        try(BufferedReader br=IOUtils.getTextReader(inputMD5File)){
            while((line=br.readLine())!=null){
                md5ValuePath.add(line);
            }
        }
        catch(Exception e){
            e.printStackTrace();
        }
        if(md5ValuePath.size()<numThreads){
            numThreads=md5ValuePath.size();
        }
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetSize(md5ValuePath.size(), numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->{
                String valuePath=md5ValuePath.get(index);
                List<String> md5Value= PStringUtils.fastSplit(valuePath,"  ");
                String value=md5Value.get(0);
                String path=inputDir+"/"+md5Value.get(1);
                boolean f=MD5.checkMD5(path, value);
                if(!f){
                    System.out.println("False  "+path+": "+value);
                }
                else{
                    //System.out.println("True  "+path+": "+value);
                }
            });
        }
        System.out.println("completed in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
    }

    public static boolean checkTwoFileMD5(String inputFile1, String inputFile2){
        String md5value1=MD5.getMD5FromFile(inputFile1);
        boolean res=MD5.checkMD5(inputFile2, md5value1);
        if(res) {
            System.out.println("These two files are equal");
        }
        else {
            System.out.println("!!! These two files are not equal");
        }
        return res;
    }

}

