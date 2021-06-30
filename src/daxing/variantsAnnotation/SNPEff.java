package daxing.variantsAnnotation;

import daxing.common.IOTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

/**
 * this class were used to extract effect and impact sub-fields from Ann field added by SnpEff
 */
public class SNPEff {

    public static void extractEffectAndImpact(String inputDir, String outDir){
        List<File> fileList = IOTool.getVisibleFileListInDir(inputDir);
        String[] outNames= fileList.stream().map(File::getName).map(s->s.replaceAll(".ann.txt.gz",".ann.effect.txt" +
                ".gz")).toArray(String[]::new);
        IntStream.range(0, fileList.size()).parallel().forEach(e->extractEffectAndImpact(fileList.get(e), new File(outDir,
                outNames[e])));
    }

    private static void extractEffectAndImpact(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw  = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, tem, te;
            StringBuilder sb=new StringBuilder();
            bw.write("Chr\tPos\tEffect\tImpact");
            bw.newLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                tem= PStringUtils.fastSplit(temp.get(5), ";");
                te= PStringUtils.fastSplit(tem.get(9), "|");
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(temp.get(2)).append("\t");
                sb.append(te.get(1)).append("\t").append(te.get(2));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
