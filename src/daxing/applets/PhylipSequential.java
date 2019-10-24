package daxing.applets;

import format.position.ChrPos;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import utils.IOUtils;
import utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;

public class PhylipSequential {

    public static void toPhylipSequentialFormat(String vcfInputDir, String chrPosRefTaxonInputDir,
                                                  String phylipSequencialOutFile){
        File[] files1=IOUtils.listRecursiveFiles(new File(vcfInputDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(chrPosRefTaxonInputDir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String res=null;
        StringBuilder sb=new StringBuilder();
        for (int i = 0; i < f1.length; i++) {
            sb.append(PhylipSequential.toPhylipSequentialFormat(f1[i], f2[i]));
        }
        res=PStringUtils.getMultiplelineString(50, sb.toString());
        try (BufferedWriter bw = IOUtils.getTextWriter(phylipSequencialOutFile)) {
            bw.write("barley.........");
            bw.write(res);
            bw.newLine();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    /**
     *
     * @param chrPosRefTaxonInputFile  file     CHR     POS     ref     hordeum_vulgare
     *                                         1       195472  g       G
     *                                         1       195473  c       C
     *                                         1       195474  t       T
     * @param vcfInputFile chr
     */
    public static String toPhylipSequentialFormat(File vcfInputFile, File chrPosRefTaxonInputFile){
        StringBuilder sb=null;
        try (BufferedReader br1 = IOUtils.getTextReader(vcfInputFile.getAbsolutePath());
             BufferedReader br2 = IOUtils.getTextReader(chrPosRefTaxonInputFile.getAbsolutePath())) {
            String line;
            List<String> lineList;
            TIntHashSet chrs=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList();
            while ((line= br1.readLine()).startsWith("##")){}
            while ((line= br1.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                chrs.add(Integer.parseInt(lineList.get(0)));
                posList.add(Integer.parseInt(lineList.get(1)));
            }
            if (chrs.size()>1){
                System.out.println("please check your "+vcfInputFile.getName()+" it has duplicated chromosome");
                System.exit(1);
            }
            List<ChrPos> chrPosList=new ArrayList<>(posList.size());
            br2.readLine();
            int chr, pos;
            int index=Integer.MIN_VALUE;
            while ((line=br2.readLine())!=null){
                sb=new StringBuilder();
                lineList=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(lineList.get(0));
                pos=Integer.parseInt(lineList.get(1));
                if (chr!=chrs.iterator().next()){
                    System.out.println("please check your "+vcfInputFile.getName()+" and "+chrPosRefTaxonInputFile);
                    System.exit(1);
                }
                index=posList.binarySearch(pos);
                if (index<0) {
                    sb.append("N");
                }
                sb.append(lineList.get(3));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return sb.toString();
    }
}
