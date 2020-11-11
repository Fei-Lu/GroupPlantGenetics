package daxing.applets;

import com.google.common.collect.Table;
import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import org.apache.commons.lang.StringUtils;
import org.apache.commons.lang.math.NumberUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ABBA_BABA {

    public static void convertVCFToGenoFormat(String vcfDir, String ancestralVCFDir, String genoOutDir,
                                              int indexOfAncestralAllele){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        long start=System.nanoTime();
        List<File> files1= IOUtils.getFileListInDirEndsWith(vcfDir,"gz");
        Comparator<File> comparator= Comparator.comparing(file -> Integer.parseInt(PStringUtils.fastSplit(file.getName(),".").get(0).substring(3)));
        Collections.sort(files1, comparator);
        List<File> files2= IOUtils.getVisibleFileListInDir(ancestralVCFDir);
        String[] outNames= files1.stream().map(File::getName).map(str->str.replaceAll("vcf.gz", "geno")).toArray(String[]::new);
        IntStream.range(0, files1.size()).forEach(e-> convertVCFToGenoFormat(files1.get(e), files2.get(e),
                new File(genoOutDir, outNames[e]), indexOfAncestralAllele));
        System.out.println("completed in "+Benchmark.getTimeSpanHours(start)+" hours");
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    private static void convertVCFToGenoFormat(File vcfFile, File ancestralFile, File outFile,
                                               int indexOfAncestralAllele){
        long start=System.nanoTime();
        Table<String, String, String> outGroupChrPosAllele= RowTableTool.getTable(ancestralFile.getAbsolutePath(), indexOfAncestralAllele);
        try (BufferedReader br1 = IOTool.getReader(vcfFile);
             BufferedWriter bw= IOTool.getWriter(outFile.getAbsolutePath())) {
            String line, genotype;
            List<String> temp;
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
                for (int i = 9; i < temp.size(); i++) {
                    genotype=temp.get(i).substring(0, 3);
                    switch (genotype){
                        case "0/1":
                            sb.append(refAllele).append("/").append(altAllele).append("\t");
                            break;
                        case "0/0":
                            sb.append(refAllele).append("/").append(refAllele).append("\t");
                            break;
                        case "1/1":
                            sb.append(altAllele).append("/").append(altAllele).append("\t");
                            break;
                        case "./.":
                            sb.append("N/N").append("\t");
                            break;
                    }
                }
                sb.append(outgroupAllele).append("/").append(outgroupAllele);
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

    /**
     * remove rows which has nan or its D value is negative
     * remove rows which fd < 0 or fd > 1
     * @param inputFdResDir
     * @param outDir
     */
    public static void prepareData_transform0(String inputFdResDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputFdResDir);
        String[] outNames= files.stream().map(File::getName).map(str->str.replaceAll("csv$", "txt")).toArray(String[]::new);
        RowTableTool<String> rowTable;
        List<String> d_List;
        List<String> fd_List;
        Predicate<List<String>> p=l->l.contains("nan");
        for (int i = 0; i < files.size(); i++) {
            rowTable=new RowTableTool<>(files.get(i).getAbsolutePath(), ",");
            rowTable.removeIf(p);
            d_List=rowTable.getColumn(8);
            fd_List=rowTable.getColumn(9);
            for (int j = 0; j < d_List.size(); j++) {
                if (!(NumberUtils.isNumber(d_List.get(j)) && NumberUtils.isNumber(fd_List.get(j)))){
                    rowTable.removeRow(j);
                    d_List.remove(j);
                    fd_List.remove(j);
                    j--;
                    continue;
                }
                double fd=Double.parseDouble(fd_List.get(j));
                double d=Double.parseDouble(d_List.get(j));
                if (d<0 || fd<0){
                    rowTable.removeRow(j);
                    d_List.remove(j);
                    fd_List.remove(j);
                    j--;
                    continue;
                }
                if ( fd > 1){
                    rowTable.removeRow(j);
                    d_List.remove(j);
                    fd_List.remove(j);
                    j--;
                    continue;
                }
            }
            rowTable.setColumn(9, fd_List);
            rowTable.writeTextTable(new File(outDir, outNames[i]).getAbsolutePath(), IOFileFormat.Text);
        }
    }

    private enum P2{
//        Landrace, Cultivar
        Landrace_Europe, Landrace_WestAsia, Landrace_EastAsia, Cultivar
    }

    private enum P3{
        WildEmmer,  DomesticatedEmmer, FreeThreshTetraploid, tauschii
    }

    public static void prepareData_merge(String inputFdResDir, String  outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputFdResDir);
        List<Predicate<File>> predicateList=new ArrayList<>();
        List<String> outNames=new ArrayList<>();
        for (P2 p2:P2.values()) {
            for(P3 p3: P3.values()){
                predicateList.add(f->f.getName().contains(p2.name()) && f.getName().contains(p3.name()));
                outNames.add(p2.name()+"-"+p3.name()+".txt");
            }
        }
        List<File> f;
        RowTableTool<String> rowTable1, rowTable2;
        for (int i = 0; i < predicateList.size(); i++) {
            f=files.stream().filter(predicateList.get(i)).collect(Collectors.toList());
            rowTable1=new RowTableTool<>(f.get(0).getAbsolutePath());
            for (int j = 1; j < f.size(); j++) {
                rowTable2=new RowTableTool<>(f.get(j).getAbsolutePath());
                rowTable1.add(rowTable2);
            }
            rowTable1.write(new File(outDir, outNames.get(i)).getAbsolutePath());
        }
    }
}
