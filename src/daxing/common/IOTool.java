package daxing.common;

import pgl.infra.utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * @author xudaxing
 */
public class IOTool extends IOUtils {

    public static BufferedReader getReader(String file){
        BufferedReader br;
        if (file.endsWith("gz")){
            br=IOUtils.getTextGzipReader(file);
        }else {
            br= IOUtils.getTextReader(file);
        }
        return br;
    }

    public static BufferedReader getReader(File file){
        BufferedReader br;
        if (file.getName().endsWith("gz")){
            br=IOUtils.getTextGzipReader(file.getAbsolutePath());
        }else {
            br=IOUtils.getTextReader(file.getAbsolutePath());
        }
        return br;
    }

    public static BufferedWriter getTextWriter(File file){
       return IOUtils.getTextWriter(file.getAbsolutePath());
    }

    public static BufferedWriter getTextGzipWriter(File file){
        return IOUtils.getTextGzipWriter(file.getAbsolutePath());
    }

    public static List<File> getVisibleFileRecursiveDir(String dir){
        File[] files= IOUtils.listRecursiveFiles(new File(dir));
        Predicate<File> hidden=File::isHidden;
        List<File> res= Arrays.stream(files).filter(hidden.negate()).sorted().collect(Collectors.toList());
        return res;
    }

    public static List<File> getVisibleDir(String dir){
        File[] files=new File(dir).listFiles();
        Predicate<File> hidden=File::isHidden;
        return Arrays.stream(files).filter(hidden.negate()).sorted().collect(Collectors.toList());
    }
}
