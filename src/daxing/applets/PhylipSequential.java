package daxing.applets;

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

public class PhylipSequential {

    public static void toPhylipSequentialFormat(String vcfInputDir, String chrPosRefTaxonDir, String outFile){
        File[] files1=IOUtils.listRecursiveFiles(new File(vcfInputDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(chrPosRefTaxonDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        StringBuilder sb=new StringBuilder();
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
        String res=null;
        for (int i = 0; i < f1.length; i++) {
            res=PhylipSequential.toPhylipSequentialFormat(f1[i], f2[i]);
            sb.append(res);
        }
        String[] line50=PStringUtils.getMultilineString(50, res);
        sb=new StringBuilder();
        sb.append("barley.........").append(line50[0]);
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write(sb.toString());
            bw.newLine();
            for (int i = 1; i < line50.length; i++) {
                sb=new StringBuilder();
                sb.append("...............").append(line50[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
        }catch (Exception e){
            e.printStackTrace();
        }

    }

    /**
     *
     * @param vcfInputFile chr
     */
    private static String toPhylipSequentialFormat(File vcfInputFile, File chrPosRefTaxonFile){
        StringBuilder sb=new StringBuilder();
        try (BufferedReader br1 = IOUtils.getTextGzipReader(vcfInputFile.getAbsolutePath());
             BufferedReader br2 = IOUtils.getTextReader(chrPosRefTaxonFile.getAbsolutePath())) {
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
            while ((line=br2.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                if (chr!=chrs.iterator().next()){
                    System.out.println("please check your "+vcfInputFile.getName()+" and "+chrPosRefTaxonFile.getName());
                    System.exit(1);
                }
                index=posList.binarySearch(pos);
                if (index<0){
                    sb.append("N");
                }else {
                    sb.append(lineList.get(3).toUpperCase());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sb.toString();
    }

//    public static void main(String[] args) {
//        PhylipSequential.toPhylipSequentialFormat(args[0], args[1], args[2]);
//    }
}
