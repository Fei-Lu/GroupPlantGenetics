package daxing.applets;

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

    public static void convertVCFToGenoFormat(String vcfDir, String genoOutDir){
        File[] files= IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> hidden=File::isHidden;
        File[] f= Arrays.stream(files).filter(hidden.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f).map(File::getName).map(str->str.replaceAll("vcf\\.gz", "geno")).toArray(String[]::new);
        IntStream.range(0, f.length).parallel().forEach(e->{
            try (BufferedReader br = IOUtils.getTextGzipReader(f[e].getAbsolutePath());
                 BufferedWriter bw= IOUtils.getTextWriter(new File(genoOutDir, outNames[e]).getAbsolutePath())) {
                String line;
                List<String> temp;
                StringBuilder sb=new StringBuilder(100);
                while ((line=br.readLine()).startsWith("##")){}
                temp=PStringUtils.fastSplit(line);
                temp=temp.stream().skip(9).collect(Collectors.toList());
                sb.append("#CHROM").append("\t").append("POS").append("\t").append(StringUtils.join(temp, "\t"));
                bw.write(sb.toString());
                bw.newLine();
                String refAllele, altAllele;
                while ((line=br.readLine())!=null){
                    temp= PStringUtils.fastSplit(line);
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
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                bw.flush();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        });
    }
}
