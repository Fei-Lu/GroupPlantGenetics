package daxing.applets;

import utils.IOUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.function.Predicate;

public class Script {

    public static void getLastz(String wheatInputDir, String outgroupInputDir, String outMAFDir, String outSH){
        File[] files1= IOUtils.listRecursiveFiles(new File(wheatInputDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(outgroupInputDir));
        Predicate<File> p= e->e.getName().endsWith("fa");
        File[] f1= Arrays.stream(files1).filter(p).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p).toArray(File[]::new);
        try(BufferedWriter bw=IOUtils.getNIOTextWriter(outSH)){
            StringBuilder sb;
            sb=new StringBuilder();
            for (int i = 0; i < f1.length; i++) {
                for (int j = 0; j < f2.length; j++) {
                    sb.append("lastz ").append(f1[i].getAbsolutePath()).append(" ").append(f2[j].getAbsolutePath())
                            .append(" --notransition --ambiguous=iupac --step=20 --nogapped --format=maf > ").append(outMAFDir)
                            .append("/Ta_").append(f1[i].getName().substring(0, 6)).append("_vs_At_")
                            .append(f2[j].getName().substring(0,5)).append(".maf"+" &");
                    bw.write(sb.toString());
                    bw.newLine();
                    sb=new StringBuilder();
                }
            }
            bw.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }
}
