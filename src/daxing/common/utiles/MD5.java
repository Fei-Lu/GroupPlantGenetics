package daxing.common.utiles;

import daxing.common.sh.CommandUtils;
import org.apache.commons.lang3.tuple.ImmutablePair;
import org.apache.commons.lang3.tuple.Pair;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.security.MessageDigest;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

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
            for (byte b : digest) {
                sb.append(Integer.toString((b & 0xff) + 0x100, 16).substring(1));
            }
            md5Value=sb.toString();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        return md5Value;
    }


    private static Pair<Path, String> calculateMD5(String inputDir, File file){
        Path absolutePath = file.toPath();
        Path reaPath = absolutePath.subpath(Paths.get(inputDir).getNameCount(), absolutePath.getNameCount());
        String md5Value = MD5.getMD5FromFile(file.getAbsolutePath());
        return new ImmutablePair<>(reaPath, md5Value);
    }

    /**
     * 以标准MD5文件的形式返回一个目录下所有文件的MD5值。在服务器上默认使用32个核；在个人电脑上默认使用全部核
     * @param inputDir 输入目录的绝对路径
     */
    public static void getMD5FromDir(String inputDir){
        long start = System.nanoTime();
        int numThreads = Math.min(Runtime.getRuntime().availableProcessors(), 32);
        List<File> files = IOTool.getVisibleFileRecursiveDir(inputDir);
        if(files.size()<numThreads){
            numThreads=files.size();
        }
        List<Callable<Pair<Path, String>>> callableList = new ArrayList<>();
        for (File file : files){
            callableList.add(()->MD5.calculateMD5(inputDir, file));
        }
        System.out.println("using "+numThreads+" threads in parallel");
        List<Pair<Path, String>> results = CommandUtils.run_commands(callableList, numThreads);
        StringBuilder sb = new StringBuilder();
        try (BufferedWriter bw = IOTool.getWriter(new File(inputDir, "md5.txt"))) {
            for (Pair<Path, String> pair : results){
                sb.setLength(0);
                sb.append(pair.getValue()).append("  ").append(pair.getKey().toString());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
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

    private static Pair<File, Boolean> getCheckMD5Res(String inputFile, String hash){
        boolean res = MD5.checkMD5(inputFile, hash);
        return new ImmutablePair<>(new File(inputFile), res);
    }

    /**
     * 对标准MD5文件进行校验(MD5文件需在等待检测MD5值的文件目录下)。在服务器上默认使用32个核；在个人电脑上默认使用全部核
     * @param inputMD5File 输入的标准MD5文件的绝对路径
     */
    public static void checkMD5(String inputMD5File){
        String inputDir=new File(inputMD5File).getParent();
        long start = System.nanoTime();
        int numThreads = Math.min(Runtime.getRuntime().availableProcessors(), 32);
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
        List<Callable<Pair<File, Boolean>>> callableList = new ArrayList<>();
        List<String> temp;
        for (String md5Path : md5ValuePath){
            temp = PStringUtils.fastSplit(md5Path, "  ");
            String md5 = temp.get(0);
            String path = new File(inputDir, temp.get(1)).getAbsolutePath();
            callableList.add(()->MD5.getCheckMD5Res(path, md5));
        }
        List<Pair<File, Boolean>> results = CommandUtils.run_commands(callableList, numThreads);
        int failedCount = 0;
        for (Pair<File, Boolean> pair : results){
            if (pair.getValue()) continue;
            failedCount++;
            System.out.println(pair.getKey().getAbsolutePath()+" check md5sum failed");
        }
        if (failedCount > 0){
            System.out.println("Total failed file count is "+failedCount);
        }else {
            System.out.println("All file passed md5 check sum");
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

