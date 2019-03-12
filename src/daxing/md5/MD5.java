package daxing.md5;

import utils.Benchmark;
import utils.IOUtils;
import utils.PArrayUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.*;

/**
 *
 * @author xudaxing
 */
public class MD5 {

    public static String getMD5(String inputFile){
        String md5Value=null;
        MessageDigest md=null;
        FileInputStream fis=null;
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
            md5Value=sb.toString().toUpperCase();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        return md5Value;
    }

    public static void getMD5(String inputDir, String outfileMD5){
        long start = System.nanoTime();
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        File[] file = IOUtils.listRecursiveFiles(new File(inputDir));
        if(file.length<numThreads){
            numThreads=file.length;
        }
        Map<String,String> md5ValuePathMap=new Hashtable<>();
        BufferedWriter bw=null;
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetNumber(file.length, numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices)
                    .filter(index-> (!(file[index].getName().contains(".DS_Store"))))
                    .forEach(index-> {
                        String md5Value=MD5.getMD5(file[index].getAbsolutePath());
                        String fName = file[index].getAbsolutePath();
                        md5ValuePathMap.put(md5Value, fName);
                    });
        }
        try{
            bw=IOUtils.getTextWriter(outfileMD5);
            for(Map.Entry<String,String> entry:md5ValuePathMap.entrySet()){
                //System.out.println(entry.getKey()+"  "+entry.getValue());
                bw.write(entry.getKey()+"  "+entry.getValue());
                bw.newLine();
            }
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
        System.out.println("getMD5 from "+inputDir+" is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
    }

    public static boolean checkMD5(String inputFile, String hash){
        String fileHash=MD5.getMD5(inputFile);
        return fileHash.equals(hash);
    }

    public static void checkMD5(String inputHashFile){
        long start = System.nanoTime();
        int numThreads = 32;
        if(Runtime.getRuntime().availableProcessors()<32) {
            numThreads=Runtime.getRuntime().availableProcessors();
        }
        List<String> md5ValuePath=new ArrayList<>();
        String line;
        try(BufferedReader br=IOUtils.getTextReader(inputHashFile)){
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
        int[][] indices=PArrayUtils.getSubsetsIndicesBySubsetNumber(md5ValuePath.size(), numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            Arrays.stream(subLibIndices).forEach(index->{
                String valuePath=md5ValuePath.get(index);
                String[] md5Value=valuePath.split("  ");
                String value=md5Value[0];
                String path=md5Value[1];
                boolean f=MD5.checkMD5(path, value);
                if(!f){
                    System.out.println("False  "+path+": "+value);
                }
                else{
                    //System.out.println("True  "+path+": "+value);
                }
            });
        }
        System.out.println("CheckMD5 is complemeted in " + String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) + " minutes");
    }

}

