package daxing.applets;

import com.google.common.collect.Table;
import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import org.apache.commons.lang.StringUtils;
import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ABBA_BABA {

    public static void convertVCFToGenoFormat(String vcfDir, String ancestralVCFDir, String genoOutDir,
                                              int indexOfAncestralAllele){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        long start=System.nanoTime();
        List<File> files1= IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> files2= IOUtils.getVisibleFileListInDir(ancestralVCFDir);
        String[] outNames=files1.stream().map(File::getName).map(str->str.replaceAll("vcf\\.gz", "geno")).toArray(String[]::new);
        IntStream.range(0, files1.size()).forEach(e-> convertVCFToGenoFormat(files1.get(e), files2.get(e),
                new File(genoOutDir, outNames[e]), indexOfAncestralAllele));
        System.out.println("completed in "+Benchmark.getTimeSpanHours(start)+" hours");
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    private static void convertVCFToGenoFormat(File vcfFile, File ancestralFile, File outFile,
                                               int indexOfAncestralAllele){
        long start=System.nanoTime();
        Table<String, String, String> outGroupChrPosAllele= RowTableTool.getTable(ancestralFile.getAbsolutePath(), 2);
        try (BufferedReader br1 = IOTool.getReader(vcfFile);
             BufferedWriter bw= IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> temp, tem;
            StringBuilder sb=new StringBuilder(100);
            while ((line=br1.readLine()).startsWith("##")){}
            temp=PStringUtils.fastSplit(line);
            temp=temp.stream().skip(9).collect(Collectors.toList());
            sb.append("#CHROM").append("\t").append("POS").append("\t").append(StringUtils.join(temp, "\t"));
            sb.append("\t").append("ancestral");
            bw.write(sb.toString());
            bw.newLine();
            String refAllele, altAllele, outgroupAllele;
            while ((line=br1.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if (!outGroupChrPosAllele.contains(temp.get(0), temp.get(1))) continue;
                outgroupAllele=outGroupChrPosAllele.get(temp.get(0), temp.get(1));
                if (outgroupAllele.equals("-")) continue;
                if (outgroupAllele.equals("NA")) continue;
                refAllele=temp.get(3);
                altAllele=temp.get(4);
                sb=new StringBuilder(1000);
                sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t");
                tem=temp.stream().skip(9).map(str->StringUtils.split(str, ":")[0])
                        .map(str->str.replaceAll("\\./\\.", "N/N")).collect(Collectors.toList());
                Collections.replaceAll(tem, "0/0", refAllele+"/"+refAllele);
                Collections.replaceAll(tem, "0/1", refAllele+"/"+altAllele);
                Collections.replaceAll(tem, "1/1", altAllele+"/"+altAllele);
                sb.append(StringUtils.join(tem, "\t"));
                sb.append("\t").append(outgroupAllele+"/"+outgroupAllele);
                bw.write(sb.toString());
                bw.newLine();
            }
            br1.close();
            bw.flush();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        System.out.println(outFile.getName()+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
    }
}
