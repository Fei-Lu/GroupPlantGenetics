package daxing.common;

import utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;

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
}
