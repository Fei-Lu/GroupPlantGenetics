package daxing.common;

import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
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

    public static BufferedWriter getWriter(File file){
       if (file.getName().endsWith("gz")) return IOUtils.getTextGzipWriter(file.getAbsolutePath());
       return IOUtils.getTextWriter(file.getAbsolutePath());
    }

    public static BufferedWriter getWriter(String file){
        return getWriter(new File(file));
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
}
