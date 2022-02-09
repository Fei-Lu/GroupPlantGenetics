package daxing.v2.annotation;

import daxing.common.utiles.IOTool;
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
        String[] outNames= fileList.stream().map(File::getName).map(s->s.replaceAll("vcf.gz",
                "effect.txt.gz")).toArray(String[]::new);
        IntStream.range(0, fileList.size()).parallel().forEach(e->extractEffectAndImpact(fileList.get(e), new File(outDir,
                outNames[e])));
    }

    private static void extractEffectAndImpact(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw  = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, tem, te, t;
            StringBuilder sb=new StringBuilder();
            //  Exon or Intron rank / total number of exons or introns.
            bw.write("Chr\tPos\tGene\tTranscript\tEffect\tImpact\tRank_Total");
            bw.newLine();
            while ((line=br.readLine()).startsWith("#")) continue;
            temp= PStringUtils.fastSplit(line);
            tem= PStringUtils.fastSplit(temp.get(7), ";");
            te= PStringUtils.fastSplit(tem.get(7), ",");
            for (int i = 0; i < te.size(); i++) {
                t= PStringUtils.fastSplit(te.get(i), "|");
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                sb.append(t.get(4)).append("\t").append(t.get(6)).append("\t");
                sb.append(t.get(1)).append("\t").append(t.get(2)).append("\t");
                sb.append(t.get(8));
                bw.write(sb.toString());
                bw.newLine();
            }
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                tem= PStringUtils.fastSplit(temp.get(7), ";");
                te= PStringUtils.fastSplit(tem.get(7), ",");
                for (int i = 0; i < te.size(); i++) {
                    t= PStringUtils.fastSplit(te.get(i), "|");
                    sb.setLength(0);
                    sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                    sb.append(t.get(4)).append("\t").append(t.get(6)).append("\t");
                    sb.append(t.get(1)).append("\t").append(t.get(2)).append("\t");
                    sb.append(t.get(8));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
