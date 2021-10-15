package daxing.common.utiles;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * @author xudaxing
 */
public class IOTool extends IOUtils {

    public static BufferedReader getReader(File file){
        BufferedReader br;
        if (file.getName().endsWith("gz")){
            br=IOUtils.getTextGzipReader(file.getAbsolutePath());
        }else {
            br=IOUtils.getTextReader(file.getAbsolutePath());
        }
        return br;
    }

    public static BufferedReader getReader(String file){
        return getReader(new File(file));
    }

    public static BufferedReader getReader(Process p){
       return new BufferedReader(new InputStreamReader(p.getInputStream()));
    }

    public static BufferedReader getErrorReader(Process p){
        return new BufferedReader(new InputStreamReader(p.getErrorStream()));
    }

    public static BufferedWriter getWriter(File file){
       if (file.getName().endsWith("gz")) return IOUtils.getTextGzipWriter(file.getAbsolutePath());
       return IOUtils.getTextWriter(file.getAbsolutePath());
    }

    public static BufferedWriter getWriter(String file){
        return getWriter(new File(file));
    }

    public static DataInputStream getBinaryReader(File file){
        if (file.getName().endsWith("gz")) return IOUtils.getBinaryGzipReader(file.getAbsolutePath());
        return IOUtils.getBinaryReader(file.getAbsolutePath());
    }

    public static DataInputStream getBinaryReader(String file){
        return getBinaryReader(new File(file));
    }

    public static DataOutputStream getBinaryWriter(File file){
        if (file.getName().endsWith("gz")) return IOUtils.getBinaryGzipWriter(file.getAbsolutePath());
        return IOUtils.getBinaryWriter(file.getAbsolutePath());
    }

    public static DataOutputStream getBinaryWriter(String file){
        return getBinaryWriter(new File(file));
    }

    /**
     * 递归获取当前目录下的所有非隐藏文件
     * @param dir dir
     * @return 当前目录下的所有非隐藏文件
     */
    public static List<File> getVisibleFileRecursiveDir(String dir){
        File[] files= IOUtils.listRecursiveFiles(new File(dir));
        Predicate<File> hidden=File::isHidden;
        List<File> res= Arrays.stream(files).filter(hidden.negate()).sorted().collect(Collectors.toList());
        return res;
    }

    /**
     * 获取当前目录下的所有非隐藏文件(不包含目录)，不递归
     * @param dir dir
     * @return
     */
    public static List<File> getVisibleDir(String dir){
        File[] files=new File(dir).listFiles();
        Predicate<File> hidden=File::isHidden;
        Predicate<File> file=File::isFile;
        return Arrays.stream(files).filter(hidden.negate().and(file)).sorted().collect(Collectors.toList());
    }

    /**
     * 查看文件header
     * @param file
     * @param delimiter
     */
    public static void viewHeader(String file, String delimiter){
        try (BufferedReader br = IOTool.getReader(file)) {
            String line=br.readLine();
            List<String> header= PStringUtils.fastSplit(line, delimiter);
            for (int i = 0; i < header.size(); i++) {
                System.out.println(i+" "+header.get(i));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void viewHeader(String file){
        viewHeader(file,"\t");
    }

    public static List<String> readAllLines(String file){
        List<String> lines=new ArrayList<>();
        String line;
        try (BufferedReader br = IOTool.getReader(file)) {
            while ((line=br.readLine())!=null){
                lines.add(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        assert lines!=null : " check "+file;
        return lines;
    }

    public static void writeAllLines(File outFile, List<String> lines){
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            for (String line:lines){
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<File> getLogFile(String title, String logDir, int logFileNum){
        assert logFileNum > 0 : logFileNum + " must be greater than 0";
        int digitsNum=(int)(Math.log10(logFileNum) +1);
        List<File> logFiles= new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        File file;
        for (int i = 0; i < logFileNum; i++) {
            sb.setLength(0);
            sb.append(title).append("_").append(PStringUtils.getNDigitNumber(digitsNum, i)).append(".log");
            file = new File(logDir, sb.toString());
            logFiles.add(file);
        }
        return logFiles;
    }
}
