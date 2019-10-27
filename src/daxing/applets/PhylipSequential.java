package daxing.applets;

import daxing.common.DateTime;
import daxing.common.NumberTool;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.IntStream;

public class PhylipSequential {

    public static void toPhylipSequentialFormat(String vcfInputDir, String chrPosRefTaxonDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
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
        IntStream.range(0, f1.length).parallel().forEach(e->
                PhylipSequential.toPhylipSequentialFormat(f1[e], f2[e], new File(outDir, outNames[e])));
        PhylipSequential.merge(outDir);
        System.out.println(DateTime.getDateTimeOfNow());
    }

    /**
     *
     * @param vcfInputFile chr
     */
    private static void toPhylipSequentialFormat(File vcfInputFile, File chrPosRefTaxonFile, File outFile){
        StringBuilder sb=new StringBuilder();
        try (BufferedReader br1 = IOUtils.getTextGzipReader(vcfInputFile.getAbsolutePath());
             BufferedReader br2 = IOUtils.getTextReader(chrPosRefTaxonFile.getAbsolutePath());
             BufferedWriter bw=IOUtils.getTextWriter(outFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            TIntHashSet chrs=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList();
            while ((line=br1.readLine()).startsWith("##")){}
            while ((line= br1.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chrs.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
            }
            if (chrs.size()>1){
                System.out.println("please check your "+vcfInputFile.getName()+" it has duplicated chromosome");
                System.exit(1);
            }
            br2.readLine();
            int chr=-1;
            int pos=-1;
            int index=-1;
            int count=0;
            int total=0;
            while ((line=br2.readLine())!=null){
                total++;
                if (sb.length()%50==0 && sb.length()!=0){
                    sb.append("\n");
                    bw.write(sb.toString());
                    sb=new StringBuilder();
                }
                lineList=PStringUtils.fastSplit(line);
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
                    sb.append(lineList.get(3).toUpperCase());
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

    public static void merge(String outDir){
        File[] files=IOUtils.listRecursiveFiles(new File(outDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        BufferedReader[] brs=new BufferedReader[f1.length];
        for (int i = 0; i < brs.length; i++) {
            brs[i]=IOUtils.getTextReader(f1[i].getAbsolutePath());
        }
        StringBuilder sb=new StringBuilder();
        int i=Integer.MIN_VALUE;
        try (BufferedWriter bw = IOUtils.getTextWriter(new File(outDir, "barley.txt").getAbsolutePath())) {
            String line;
            int sbLen;
            line=brs[0].readLine();
            sb.append("Barley.........").append(line);
            if (sb.length()==65){
                bw.write(sb.toString());
                bw.newLine();
                sb=new StringBuilder();
            }
            for (i = 0; i < brs.length; i++) {
                while ((line=brs[i].readLine())!=null){
                    if (sb.length()==0){
                        sb.append("...............").append(line);
                        if (sb.length()==65){
                            bw.write(sb.toString());
                            bw.newLine();
                            sb=new StringBuilder();
                        }
                    }else if (sb.length()<65){
                        if (sb.length()+line.length()<65){
                            sb.append(line);
                            continue;
                        }
                        sb.append(line, 0, (65-sb.length()));
                        bw.write(sb.toString());
                        bw.newLine();
                        sbLen=sb.length();
                        sb=new StringBuilder();
                        sb.append("...............");
                        sb.append(line, 65-sbLen, line.length());
                    }else if (sb.length()==65){
                        bw.write(sb.toString());
                        bw.newLine();
                        sb=new StringBuilder();
                    }
                }
            }
            bw.write(sb.toString());
            bw.newLine();
        }catch (Exception e){
            System.out.println(sb.toString());
            System.out.println(i);
            e.printStackTrace();
        }
    }

//    public static void main(String[] args) {
//        PhylipSequential.toPhylipSequentialFormat(args[0], args[1], args[2]);
//    }
}
