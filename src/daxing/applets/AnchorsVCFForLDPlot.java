package daxing.applets;

import daxing.common.ArrayTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class AnchorsVCFForLDPlot {

    public static void getAnchorPoint(String poloid, String chrSizeFile, int numberOfAnchors, String outDir){
        try (BufferedReader br = IOUtils.getTextReader(chrSizeFile)) {
            String line;
            List<String> temp;
            br.readLine();
            br.readLine();
            int count=0;
            TIntArrayList chrAllSize=new TIntArrayList(21);
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                chrAllSize.add(Integer.parseInt(temp.get(1)));
                count++;
                if (count==21) break;
            }
            List<String> chrs=null;
            TIntArrayList chrSize=null;
            switch (poloid.toUpperCase()){
                case "D":
                    chrs=WheatLineage.valueOf("D").getChr();
                    chrSize=new TIntArrayList(7);
                    for (int i = 2; i < chrAllSize.size(); i=i+3) {
                        chrSize.add(chrAllSize.get(i));
                    }
                    break;
                case "AB":
                    chrs=WheatLineage.abLineage();
                    chrSize=new TIntArrayList(14);
                    for (int i = 0; i < chrAllSize.size(); i++) {
                        if ((i-2)%3==0) continue;
                        chrSize.add(chrAllSize.get(i));
                    }
                    break;
                case "ABD":
                    chrs=WheatLineage.abdLineage();
                    chrSize=new TIntArrayList(21);
                    for (int i = 0; i < chrAllSize.size(); i++) {
                        chrSize.add(chrAllSize.get(i));
                    }
                    break;
                default:
                    System.out.println("please input D AB or ABD");
                    System.exit(1);
            }
            BufferedWriter[] bws=new BufferedWriter[chrs.size()];
            for (int i = 0; i < chrs.size(); i++) {
                bws[i]= IOUtils.getTextWriter(new File(outDir, "chr"+chrs.get(i)+"_anchorPoint.txt").getAbsolutePath());
                bws[i].write("Chr\tPos\n");
            }
            double[] rate= ArrayTool.getElementPercent(chrSize.toArray());
            int[] anchorsNum=new int[chrSize.size()];
            for (int i = 0; i < rate.length; i++) {
                anchorsNum[i]= (int) (numberOfAnchors*rate[i]);
                System.out.println(anchorsNum[i]);
            }
            int[] randomNum;
            StringBuilder sb;
            for (int i = 0; i < chrSize.size(); i++) {
                randomNum=ArrayTool.getRandomNonrepetitionArray(anchorsNum[i], 0, chrSize.get(i)+1);
                Arrays.sort(randomNum);
                for (int j = 0; j < randomNum.length; j++) {
                    sb=new StringBuilder();
                    sb.append(chrs.get(i)).append("\t").append(randomNum[j]);
                    bws[i].write(sb.toString());
                    bws[i].newLine();
                }
                bws[i].flush();
                bws[i].close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void getChrPosLen(String vcfDir, String outDir){
        long start=System.nanoTime();
        File[] files=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        String[] outNames=Arrays.stream(f).map(File::getName).map(str->str.replaceAll("vcf$", "chrpos.txt")).toArray(String[]::new);
        BufferedReader br;
        BufferedWriter bw;
        try {
            String line;
            List<String> temp;
            StringBuilder sb;
            for (int i = 0; i < f.length; i++) {
                br=IOUtils.getTextReader(f[i].getAbsolutePath());
                bw=IOUtils.getTextWriter(new File(outDir, outNames[i]).getAbsolutePath());
                bw.write("Chr\tPos\tLen\n");
                while ((line=br.readLine()).startsWith("##")){}
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    sb=new StringBuilder();
                    sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("\t").append(line.length());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                System.out.println(new File(outDir, outNames[i]).getName()+" was completed");
                br.close();
                bw.flush();
                bw.close();
            }
            System.out.println("completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void getAnchors(String anchorPointDir, String chrPosLenDir, String vcfDir, String anchorsOutDir,
                                  double rate, int numThread){
        long start1=System.nanoTime();
        File[] files1=IOUtils.listRecursiveFiles(new File(anchorPointDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(chrPosLenDir));
        File[] files3=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        File[] f3=Arrays.stream(files3).filter(p.negate()).toArray(File[]::new);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(f1.length, numThread);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream()
                    .forEach(index-> {
                        long start=System.nanoTime();
                        try (BufferedReader br1 = IOUtils.getTextReader(f1[index].getAbsolutePath());
                             BufferedReader br2 = IOUtils.getTextReader(f2[index].getAbsolutePath())) {
                            BufferedReader br3;
                            br1.readLine();
                            br2.readLine();
                            String line;
                            List<String> temp;
                            Set<String> chrs=new HashSet<>();
                            TIntArrayList posList=new TIntArrayList();
                            while ((line=br2.readLine())!=null){
                                temp=PStringUtils.fastSplit(line);
                                chrs.add(temp.get(0));
                                posList.add(Integer.parseInt(temp.get(1)));
                            }
                            br2.close();
                            if (chrs.size()>1){
                                System.out.println("check your "+f2[index].getName());
                                System.exit(1);
                            }
                            String chr=null;
                            int pos;
                            int index1=Integer.MIN_VALUE;
                            int index2=Integer.MIN_VALUE;
                            int startIndex=Integer.MIN_VALUE;
                            int endIndex=Integer.MIN_VALUE;
                            BufferedWriter bw;
                            int anchorNumber=0;
                            List<String> vcfLines;
                            while ((line=br1.readLine())!=null){
                                anchorNumber++;
                                temp=PStringUtils.fastSplit(line);
                                chr=temp.get(0);
                                if (!chr.equals(chrs.iterator().next())){
                                    System.out.println("check your "+f1[index].getName()+" and "+ f2[index].getName());
                                }
                                pos=Integer.parseInt(temp.get(1));
                                Anchor anchor=new Anchor(chr, pos);
                                index1=posList.binarySearch(anchor.getStart());
                                index2=posList.binarySearch(anchor.getEnd());
                                if (index1<0){
                                    startIndex=-index1-1;
                                }
                                if (index2<0){
                                    endIndex=-index2-2;
                                }
                                bw= IOUtils.getTextWriter(new File(anchorsOutDir, "chr"+chr+"."+anchor.getStart()+"_"+
                                        anchor.getEnd()+".vcf").getAbsolutePath());
                                br3=IOUtils.getTextReader(f3[index].getAbsolutePath());
                                String lin=null;
                                while ((lin=br3.readLine()).startsWith("##")){
                                    bw.write(lin);
                                    bw.newLine();
                                }
                                bw.write(lin);
                                bw.newLine();
                                long start2=System.nanoTime();
                                Predicate<String> pre=str->Double.parseDouble(StringUtils.split(str, "\t;=")[20])<0.05;
                                Predicate<String> predicate=pre.negate().and(str->Math.random()<rate);
                                vcfLines= br3.lines().skip(startIndex).limit(endIndex-startIndex).filter(predicate).collect(Collectors.toList());
                                for (int j = 0; j < vcfLines.size(); j++) {
                                    bw.write(vcfLines.get(j));
                                    bw.newLine();
                                }
                                br3.close();
                                bw.flush();
                                bw.close();
                                System.out.println("chr"+chr+" anchor "+anchorNumber+" completed in "+Benchmark.getTimeSpanSeconds(start2)+
                                        " s, only "+rate+" were retained");
                            }
                            System.out.println("chr"+chr+" complicated in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
                        }catch (Exception ex){
                            ex.printStackTrace();
                        }
                    });
        }
        System.out.println("complicated in "+Benchmark.getTimeSpanHours(start1)+" hours");
    }




//    public static void main(String[] args) {
//        AnchorsVCFForLDPlot.getAnchorPoint(args[0], args[1], Integer.parseInt(args[2]), args[3]);
//        AnchorsVCFForLDPlot.getChrPosLen(args[0], args[1]);
//        AnchorsVCFForLDPlot.getAnchors(args[0], args[1], args[2], args[3], Double.parseDouble(args[4]), Integer.parseInt(args[5]));
//    }

}
class Anchor {

    String chr;
    int start;
    int end;

    public Anchor(String chr, int start, int end){
        this.chr=chr;
        this.start=start;
        this.end=end;
    }

    public Anchor(String chr, int anchorPoint){
        this.chr=chr;
        this.end=anchorPoint+25000000;
        if (anchorPoint<25000000){
            this.start=0;
        }else {
            this.start=anchorPoint-25000000;
        }
    }

    public String getChr() {
        return chr;
    }

    public int getStart() {
        return start;
    }

    public int getEnd() {
        return end;
    }

    public int getAnchorPoint(){
        return this.getEnd()-25000000;
    }

    @Override
    public String toString() {
        return String.valueOf(this.getAnchorPoint());
    }
}

