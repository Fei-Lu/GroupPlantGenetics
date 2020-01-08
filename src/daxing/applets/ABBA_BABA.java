package daxing.applets;

import com.google.common.collect.Table;
import daxing.common.RowTableTool;
import org.apache.commons.lang.StringUtils;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ABBA_BABA {

    public static void convertVCFToGenoFormat(String vcfDir, String ancestralVCFDir, String genoOutDir){
        File[] files1= IOUtils.listRecursiveFiles(new File(vcfDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(ancestralVCFDir));
        Predicate<File> hidden=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(hidden.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(hidden.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f1).map(File::getName).map(str->str.replaceAll("vcf\\.gz", "geno")).toArray(String[]::new);
        IntStream.range(0, f1.length).parallel().forEach(e-> convertVCFToGenoFormat(f1[e], f2[e], new File(genoOutDir, outNames[e])));
    }

    private static void convertVCFToGenoFormat(File vcfFile, File ancestralFile, File outFile){
        Table<String, String, String> outGroupChrPosAllele= RowTableTool.getTable(ancestralFile.getAbsolutePath(), 5);
        try (BufferedReader br1 = IOUtils.getTextGzipReader(vcfFile.getAbsolutePath());
             BufferedWriter bw= IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> temp;
            StringBuilder sb=new StringBuilder(100);
            while ((line=br1.readLine()).startsWith("##")){}
            temp=PStringUtils.fastSplit(line);
            temp=temp.stream().skip(9).collect(Collectors.toList());
            sb.append("#CHROM").append("\t").append("POS").append("\t").append(StringUtils.join(temp, "\t"));
            bw.write(sb.toString());
            bw.newLine();
            String refAllele, altAllele, outgroupAllele;
            while ((line=br1.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if (!outGroupChrPosAllele.contains(temp.get(0), temp.get(1))) continue;
                refAllele=temp.get(3);
                altAllele=temp.get(4);
                sb=new StringBuilder(1000);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                temp=temp.stream().skip(9).map(str->StringUtils.split(str, ":")[0])
                        .map(str->str.replaceAll("\\./\\.", "N/N")).collect(Collectors.toList());
                Collections.replaceAll(temp, "0/0", refAllele+"/"+refAllele);
                Collections.replaceAll(temp, "0/1", refAllele+"/"+altAllele);
                Collections.replaceAll(temp, "1/1", altAllele+"/"+altAllele);
                sb.append(StringUtils.join(temp, "\t"));
                outgroupAllele=outGroupChrPosAllele.get(temp.get(0), temp.get(1));
                sb.append("\t").append(outgroupAllele+"/"+outgroupAllele);
                bw.write(sb.toString());
                bw.newLine();
            }
            br1.close();
            bw.flush();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
}
