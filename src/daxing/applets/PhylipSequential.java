package daxing.applets;

import daxing.common.DateTime;
import daxing.common.MD5;
import daxing.common.NumberTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class PhylipSequential {

    public static void toPhylipSequentialFormat(String vcfInputDir, String chrPosRefTaxonDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        File outFileDir=new File(outDir, "chr");
        outFileDir.mkdir();
        File[] files1=IOUtils.listRecursiveFiles(new File(vcfInputDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(chrPosRefTaxonDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String[] outNames= Arrays.stream(f1).map(File::getName).map(str->str.substring(0,6)+".txt").toArray(String[]::new);
        if (f1.length!=f2.length){
            System.out.println("please check your "+vcfInputDir+" and "+chrPosRefTaxonDir);
            System.exit(1);
        }
        int[] chrF1=Arrays.stream(f1).map(File::getName).map(str->str.substring(3,6)).mapToInt(Integer::parseInt).toArray();
        int[] chrF2=Arrays.stream(f2).map(File::getName).map(str->str.substring(3,6)).mapToInt(Integer::parseInt).toArray();
        for (int i = 0; i < chrF1.length; i++) {
            if (chrF1[i]==chrF2[i]) continue;
            System.out.println("please check your "+vcfInputDir+" and "+chrPosRefTaxonDir);
        }
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, 21);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(e-> PhylipSequential.toPhylipSequentialFormat(f1[e], f2[e], new File(outFileDir, outNames[e])));
        }
        PhylipSequential.merge(outFileDir.getAbsolutePath());
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    /**
     *
     * @param vcfInputFile chr
     */
    private static void toPhylipSequentialFormat(File vcfInputFile, File chrPosRefTaxonFile, File outFile){
        StringBuilder sb=new StringBuilder();
        try (BufferedReader br1 = IOUtils.getTextGzipReader(vcfInputFile.getAbsolutePath());
             BufferedReader br2 = IOUtils.getTextGzipReader(chrPosRefTaxonFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            TIntHashSet chrs=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList();
            List<String> alleleList=new ArrayList<>();
            br2.readLine();
            String allele;
            while ((line=br2.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                allele=lineList.get(5);
                if (!(allele.equals("A") || allele.equals("T") || allele.equals("C") || allele.equals("G"))) continue;
                chrs.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
                alleleList.add(allele);
            }
            if (chrs.size()>1){
                System.out.println("please check your "+chrPosRefTaxonFile.getName()+" it has duplicated chromosome");
                System.exit(1);
            }
            int chr=-1;
            int pos=-1;
            int index=-1;
            double count=0;
            int total=0;
            while ((line=br1.readLine()).startsWith("##")){}
            while ((line= br1.readLine())!=null){
                total++;
                if (sb.length()%50==0 && sb.length()!=0){
                    sb.append("\n");
                    bw.write(sb.toString());
                    sb=new StringBuilder();
                }
                lineList= PStringUtils.fastSplit(line);
                chr=Integer.parseInt(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                if (chr!=chrs.iterator().next()){
                    System.out.println("please check your "+vcfInputFile.getName()+" and "+chrPosRefTaxonFile.getName());
                    System.out.println(chr+"\t"+chrs.iterator().next());
                    System.exit(1);
                }
                index=posList.binarySearch(pos);
                if (index<0){
                    sb.append("N");
                }else {
                    count++;
                    sb.append(alleleList.get(index).toUpperCase());
                }
            }
            System.out.println(vcfInputFile.getName()+" "+NumberTool.format(count/total, 2)+"("+count+"/"+total+") " +
                    "snp having barley allele");
            bw.write(sb.toString());
            bw.newLine();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param outDir
     * @param chrs
     * @param outFile
     */
    private static void merge(String outDir, int[] chrs, String outFile){
        File[] files=IOUtils.listRecursiveFiles(new File(outDir));
        List<Integer> chrList= Arrays.stream(chrs).boxed().collect(Collectors.toList());
        Predicate<File> p=File::isHidden;
        Predicate<File> p2=p.negate().and(f->chrList.contains(Integer.parseInt(f.getName().substring(3,6))));
        File[] f1=Arrays.stream(files).filter(p2).toArray(File[]::new);
        BufferedReader[] brs=new BufferedReader[f1.length];
        for (int i = 0; i < brs.length; i++) {
            brs[i]=IOUtils.getTextReader(f1[i].getAbsolutePath());
        }
        StringBuilder sb=new StringBuilder();
        StringBuilder temp;
        try(BufferedWriter bw=IOUtils.getTextWriter(outFile)){
            String line;
            String[] lines;
            for (int i = 0; i < brs.length; i++) {
                while ((line=brs[i].readLine())!=null){
                    if (sb.length()>=1000){
                        lines=PStringUtils.getMultilineString(50, sb.toString());
                        for (int j = 0; j < lines.length-1; j++) {
                            temp=new StringBuilder();
                            bw.write(temp.append("               ").append(lines[j]).toString());
                            bw.newLine();
                        }
                        sb=new StringBuilder();
                        sb.append(lines[lines.length-1]);
                    }
                    sb.append(line);
                }
            }
            lines=PStringUtils.getMultilineString(50, sb.toString());
            for (int i = 0; i < lines.length; i++) {
                temp=new StringBuilder();
                bw.write(temp.append("               ").append(lines[i]).toString());
                bw.newLine();
            }
            for (int i = 0; i < brs.length; i++) {
                brs[i].close();
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private static void merge(String outDir){
        String out=new File(outDir).getParent();
        merge(outDir, WheatLineage.ablineage(), new File(out, "chrAB.subgenome.txt").getAbsolutePath());
        merge(outDir, WheatLineage.valueOf("D").getChrID(), new File(out, "chrD.subgenome.txt").getAbsolutePath());
    }

//    public static void main(String[] args) {
//        PhylipSequential.toPhylipSequentialFormat(args[0], args[1], args[2]);
//    }
}
