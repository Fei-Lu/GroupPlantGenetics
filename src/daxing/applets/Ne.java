package daxing.applets;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.*;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Ne {

    public static void retainAncestralVmapII(String vcfDir, String ancestralDir, String outDir, int numThreads){
        List<File> vcfFiles= IOUtils.getFileListInDirEndsWith(vcfDir, "gz");
        List<File> ancestralFiles= IOUtils.getFileListInDirEndsWith(ancestralDir, "gz");
        String[] outNames= vcfFiles.stream().map(File::getName).map(s->s.replaceAll("vcf\\.gz","ancestral.vcf")).toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(vcfFiles.size(), numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList= Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->
                    retainAncestralVmapII(vcfFiles.get(index), ancestralFiles.get(index), new File(outDir, outNames[index])));
        }
    }

    public static void retainAncestralVmapII(File vcfFile, File ancestralFile, File outFile){
        Table<Integer, Integer, String> ancestral= getAncestral(ancestralFile.getAbsolutePath());
        try (BufferedReader brVCF = IOTool.getReader(vcfFile);
             BufferedWriter bw=IOTool.getTextWriter(outFile)) {
            String line;
            List<String> temp;
            while ((line=brVCF.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            bw.write(line);
            bw.newLine();
            int chr, pos;
            int count=0;
            int total=0;
            while ((line=brVCF.readLine())!=null){
                total++;
                temp=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(temp.get(0));
                pos=Integer.parseInt(temp.get(1));
                if (!ancestral.contains(chr, pos)) continue;
                count++;
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
            System.out.println(vcfFile.getName()+" total: "+total+" count: "+count);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static Table<Integer, Integer, String> getAncestral(String ancestralFile){
        Table<Integer, Integer, String> table=HashBasedTable.create();
        try (BufferedReader br = IOTool.getReader(ancestralFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                table.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)), temp.get(2));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return table;
    }

    public static void transformRefToAncestral(String vcfDir, String ancestralDir, String outDir, int numThreads){
        List<File> vcfFiles=IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> ancestralFiles=IOUtils.getVisibleFileListInDir(ancestralDir);
        String[] outNames=vcfFiles.stream().map(File::getName).map(s->s.replaceAll("ancestral\\.vcf","refToAncestral" +
                ".vcf")).toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(vcfFiles.size(), numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList= Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(index->
                    transformRefToAncestral(vcfFiles.get(index), ancestralFiles.get(index), new File(outDir, outNames[index])));
        }
    }

    public static void transformRefToAncestral(File vcfFile, File ancestralFile, File outFile){
        Table<Integer, Integer, String> ancestralTable=getAncestral(ancestralFile.getAbsolutePath());
        try (BufferedReader br = IOTool.getReader(vcfFile);
             BufferedWriter bw=IOTool.getTextWriter(outFile)) {
            String line, refAllele, ancestralAllele;
            List<String> temp;
            int chr, pos;
            int total=0;
            int count=1;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            bw.write(line);
            bw.newLine();
            String te;
            while ((line=br.readLine())!=null){
                total++;
                temp=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(temp.get(0));
                pos=Integer.parseInt(temp.get(1));
                refAllele=temp.get(3);
                ancestralAllele=ancestralTable.get(chr, pos);
                if (refAllele.equals(ancestralAllele)){
                    count++;
                    bw.write(line);
                    bw.newLine();
                }else {
                    temp.set(3, ancestralAllele);
                    temp.set(4, refAllele);
                    for (int i = 9; i < temp.size(); i++) {
                        if (temp.get(i).startsWith("0/0")){
                            te=temp.get(i);
                            temp.set(i, te.replaceFirst("0/0", "1/1"));
                        }else if (temp.get(i).startsWith("1/1")){
                            te=temp.get(i);
                            temp.set(i, te.replaceFirst("1/1", "0/0"));
                        }
                    }
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                }
            }
            bw.flush();
            System.out.println(vcfFile.getName()+" total sites: "+total+" refAllele equal ancestral allele sites: "+ count);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
