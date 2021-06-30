package daxing.variantsAnnotation;

import daxing.common.IOTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.IntStream;

public class VEP {

    public static void extractEffectAndImpact(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".effect.txt");
        String[] outNames= files.stream().map(File::getName).map(s->s.replaceAll(".variant",".ann.variant")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->extractConsequenceAndImpact(files.get(e), new File(outDir,
                outNames[e])));
    }

    /**
     * Consequence is effect in SnpEff
     * @param inputFile
     * @param outFile
     */
    private static void extractConsequenceAndImpact(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp, tem, te;
            StringBuilder sb= new StringBuilder();
            bw.write("Chr\tPos\tConsequence\tImpact");
            bw.newLine();
            while ((line = br.readLine()).startsWith("##")) continue;
            while ((line = br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
//                if (temp.get(6).equals("intergenic_variant")) continue;
//                if (temp.get(6).equals("downstream_gene_variant")) continue;
//                if (temp.get(6).equals("intron_variant")) continue;
//                if (temp.get(6).equals("upstream_gene_variant")) continue;
                tem = PStringUtils.fastSplit(temp.get(1), ":");
                sb.setLength(0);
                sb.append(String.join("\t", tem)).append("\t");
                sb.append(temp.get(6)).append("\t");
                tem= PStringUtils.fastSplit(temp.get(13), ";");
                te= PStringUtils.fastSplit(tem.get(0), "=");
                sb.append(te.get(1));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void filterCDS(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".effect.txt.gz");
        String[] outNames= files.stream().map(File::getName).map(s->s.replaceAll(".ann",".vep.CDS")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->filterCDS(files.get(e), new File(outDir, outNames[e])));

    }

    private static void filterCDS(File inputFile, File outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            bw.write("Chr\tPos\tConsequence\tImpact");
            bw.newLine();
            br.readLine();
            while ((line = br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (temp.get(2).equals("intergenic_variant")) continue;
                if (temp.get(2).equals("upstream_gene_variant")) continue;
                if (temp.get(2).equals("downstream_gene_variant")) continue;
                if (temp.get(2).equals("intron_variant")) continue;
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
