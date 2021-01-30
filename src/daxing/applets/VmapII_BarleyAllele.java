package daxing.applets;

import daxing.common.NumberTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TByteArrayList;
import gnu.trove.list.array.TIntArrayList;
import gnu.trove.set.hash.TIntHashSet;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.IntStream;


public class VmapII_BarleyAllele {

    public static void extractIntersectionBetweenVmapIIandAncestral(String ancestralDir, String vmap2Dir, String outDir){
        File[] files1= IOUtils.listRecursiveFiles(new File(ancestralDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(vmap2Dir));
        Predicate<File> p=File::isHidden;
        File[] f1= Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String[] outName=Arrays.stream(f1).map(File::getName).map(str->str.substring(0, 6)).toArray(String[]::new);
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setGroupingUsed(false);
        numberFormat.setMaximumFractionDigits(5);
        IntStream.range(0, f1.length).forEach(e->{
            String line;
            List<String> temp;
            TIntHashSet chrSet=new TIntHashSet();
            TIntArrayList posList=new TIntArrayList(4000000);
            TByteArrayList barleyAlleleList=new TByteArrayList(4000000);
            try (BufferedReader br1 = IOUtils.getTextGzipReader(f1[e].getAbsolutePath());
                 BufferedReader br2 = IOUtils.getTextReader(f2[e].getAbsolutePath());
                 BufferedWriter bw  = IOUtils.getTextGzipWriter(new File(outDir, outName[e]+"barleyVmapII.vcf.gz")
                         .getAbsolutePath())) {
                br1.readLine();
                br2.readLine();
                while ((line=br1.readLine())!=null){
                    temp= PStringUtils.fastSplit(line);
                    if (!AlleleEncoder.alleleBaseToByteMap.containsKey(temp.get(5).charAt(0))) continue;
                    byte allele= AlleleEncoder.alleleBaseToByteMap.get(temp.get(5).charAt(0));
                    chrSet.add(Integer.parseInt(temp.get(0)));
                    posList.add(Integer.parseInt(temp.get(1)));
                    barleyAlleleList.add(allele);
                }
                if (chrSet.size() >1 ) {
                    System.out.println("error");
                    System.exit(1);
                }
                System.out.print("barley had "+ NumberTool.parse(barleyAlleleList.size())+" chr pos in ancestral "+ f1[e].getName().substring(0,6));
                while ((line=br2.readLine()).startsWith("##")) {
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                int countVmap2=0;
                int countBarley=0;
                while ((line=br2.readLine())!=null){
                    countVmap2++;
                    temp=PStringUtils.fastSplit(line);
                    int chr=Integer.parseInt(temp.get(0));
                    int pos=Integer.parseInt(temp.get(1));
                    if (!chrSet.contains(chr)){
                        System.out.println("error");
                        System.exit(1);
                    }
                    int index=posList.binarySearch(pos);
                    if (index<0) continue;
                    countBarley++;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
                System.out.println(", and only "+ numberFormat.format(((double)countBarley/countVmap2))
                        +"("+NumberTool.parse(countBarley)+", "+NumberTool.parse(countVmap2)+ ") " +"chr pos had barley "
                        + "allele in vmapII "+f1[e].getName().substring(0, 6));
//                System.out.println("vamp2 "+f1[e].getName().substring(0, 6)+" "+NumberTool.format(((double)countBarley/countVmap2), 5)
//                        +"("+NumberTool.parse(countBarley)+", "+NumberTool.parse(countVmap2)+ ") " +"chr pos had barley "
//                        + "allele");

            } catch (IOException ex) {
                ex.printStackTrace();
            }
        });
    }

    public static void getBarleyFasta(String vcfDir, String ancestralDir, String outDir){
        File[] files1=IOUtils.listRecursiveFiles(new File(vcfDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(ancestralDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        String[] outName=Arrays.stream(f1).map(File::getName).map(str->str.substring(0, 6)+
                "_vmapII_50K_withBarleyAllele.fasta").toArray(String[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, 7);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(e-> {
                        try (BufferedReader br1 = IOUtils.getTextGzipReader(f1[e].getAbsolutePath());
                             BufferedReader br2=IOUtils.getTextGzipReader(f2[e].getAbsolutePath());
                             BufferedWriter bw=IOUtils.getTextWriter(new File(outDir, outName[e]).getAbsolutePath())) {
                            String line;
                            List<String> temp;
                            TIntHashSet chrSet=new TIntHashSet();
                            TIntArrayList posList=new TIntArrayList();
                            while ((line=br1.readLine()).startsWith("##")) continue;
                            while ((line=br1.readLine())!=null){
                                temp=PStringUtils.fastSplit(line);
                                chrSet.add(Integer.parseInt(temp.get(0)));
                                posList.add(Integer.parseInt(temp.get(1)));
                            }
                            if (chrSet.size()>1){
                                System.out.println("error");
                                System.exit(1);
                            }
                            br2.readLine();
                            bw.write("> "+f1[e].getName().substring(0, 6)+"_barley");
                            bw.newLine();
                            int count=0;
                            while ((line=br2.readLine())!=null){
                                temp=PStringUtils.fastSplit(line);
                                int chr=Integer.parseInt(temp.get(0));
                                if (!chrSet.contains(chr)){
                                    System.out.println("error");
                                    System.exit(1);
                                }
                                int pos=Integer.parseInt(temp.get(1));
                                int index=posList.binarySearch(pos);
                                if ( index <0) continue;
                                count++;
                                bw.write(temp.get(5));
                            }
                            bw.newLine();
                            bw.flush();
                            System.out.println(f1[e].getName().substring(0, 6)+" had "+posList.size()+" chr pos with " +
                                    "barley allele "+ count);
                        } catch (IOException ex) {
                            ex.printStackTrace();
                        }
                    });
        }

    }

    public static void mergeABsubgenome(String inputDir, String outFile) {
        File[] files=IOUtils.listRecursiveFiles(new File(inputDir));
        Predicate<File> p=File::isHidden;
        int[] a_lineage= WheatLineage.ablineage();
        TIntArrayList a_list=new TIntArrayList(a_lineage);
        Predicate<File> a=f->a_list.contains(Integer.parseInt(f.getName().substring(3,6)));
        File[] f=Arrays.stream(files).filter(p.negate().and(a)).toArray(File[]::new);
        BufferedReader br;
        BufferedWriter bw=IOUtils.getTextWriter(outFile);
        String line;
        try{
            bw.write("> chrAB_barley\n");
            for (int i = 0; i < f.length; i++) {
                br=IOUtils.getTextReader(f[i].getAbsolutePath());
                br.readLine();
                while ((line=br.readLine())!=null){
                    bw.write(line);
                }
                br.close();
            }
            bw.newLine();
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

//    public static void main(String[] args) {
//        extractIntersectionBetweenVmapIIandAncestral("", "", "");
//        getBarleyFasta("", "", "");
//        mergeABsubgenome("", "");
//    }

}
