package daxing.applets;

import daxing.common.ArrayTool;
import daxing.common.RandomAccessFileTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.StringUtils;
import utils.Benchmark;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.function.Predicate;

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
                    chrs=WheatLineage.dLineage();
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

    public static void getAnchors(String anchorPointDir, String chrPosMafDir, String vcfDir, String anchorsOutDir,
                                  double rate){
        long start1=System.nanoTime();
        File[] files1=IOUtils.listRecursiveFiles(new File(anchorPointDir));
        File[] files2=IOUtils.listRecursiveFiles(new File(chrPosMafDir));
        File[] files3=IOUtils.listRecursiveFiles(new File(vcfDir));
        Predicate<File> p=File::isHidden;
        File[] f1=Arrays.stream(files1).filter(p.negate()).toArray(File[]::new);
        File[] f2=Arrays.stream(files2).filter(p.negate()).toArray(File[]::new);
        File[] f3=Arrays.stream(files3).filter(p.negate()).toArray(File[]::new);
        String lin=null;
        for (int i = 0; i < f1.length; i++) {
            long start=System.nanoTime();
            try (BufferedReader br1 = IOUtils.getTextReader(f1[i].getAbsolutePath());
                 BufferedReader br2 = IOUtils.getTextReader(f2[i].getAbsolutePath());
                 RandomAccessFile rac=new RandomAccessFile(f3[i].getAbsoluteFile(), "r")) {
                br1.readLine();
                br2.readLine();
                String line;
                List<String> temp;
                Set<String> chrs=new HashSet<>();
                TIntArrayList posList=new TIntArrayList();
                TIntArrayList lenList=new TIntArrayList();
                while ((line=br2.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrs.add(temp.get(0));
                    posList.add(Integer.parseInt(temp.get(1)));
                    lenList.add(Integer.parseInt(temp.get(2)));
                }
                br2.close();
                if (chrs.size()>1){
                    System.out.println("check your "+f2[i].getName());
                    System.exit(1);
                }
                String chr=null;
                int pos;
                int index1=Integer.MIN_VALUE;
                int index2=Integer.MIN_VALUE;
                long pointer1=Integer.MIN_VALUE;
                long pointer2=Integer.MIN_VALUE;
                BufferedWriter bw;
                int anchorNumber=0;
                while ((line=br1.readLine())!=null){
                    anchorNumber++;
                    temp=PStringUtils.fastSplit(line);
                    chr=temp.get(0);
                    if (!chr.equals(chrs.iterator().next())){
                        System.out.println("check your "+f1[i].getName()+" and "+ f2[i].getName());
                    }
                    pos=Integer.parseInt(temp.get(1));
                    Anchor anchor=new Anchor(chr, pos);
                    index1=posList.binarySearch(anchor.getStart());
                    index2=posList.binarySearch(anchor.getEnd());
                    if (index1<0){
                        index1=-index1-1;
                    }
                    if (index2<0){
                        index2=-index2-1;
                    }
                    bw= IOUtils.getTextWriter(new File(anchorsOutDir, "chr"+chr+"."+anchor.getStart()+"_"+
                            anchor.getEnd()+".vcf").getAbsolutePath());
                    long size=0;
                    rac.seek(0);
                    while ((lin=rac.readLine()).startsWith("##")){
                        bw.write(lin);
                        bw.newLine();
                        size+=lin.length()+1;
                    }
                    bw.write(lin);
                    bw.newLine();
                    size=size+lin.length()+1;
                    pointer1= RandomAccessFileTool.getPointer(lenList, index1);
                    pointer2=RandomAccessFileTool.getPointer(lenList, index2);
                    pointer1=pointer1+size;
                    pointer2=pointer2+size;
                    rac.seek(pointer1);
                    String[] te;
                    long start2=System.nanoTime();
                    while ((lin=rac.readLine())!=null && rac.getFilePointer()<pointer2){
                        te= StringUtils.split(lin, "\t;=");
                        if (Double.parseDouble(te[20])<0.05) continue;
                        if (Math.random()>rate) continue;
                        bw.write(lin);
                        bw.newLine();
                    }
                    bw.flush();
                    bw.close();
                    System.out.println("anchor "+anchorNumber+" completed in "+Benchmark.getTimeSpanMilliseconds(start2)+
                            " ms");
                }
                System.out.println("chr"+chr+" complicated in "+Benchmark.getTimeSpanMinutes(start)+" minutes");
            }catch (Exception ex){
                System.out.println(lin);
                ex.printStackTrace();
            }

        }
        System.out.println("complicated in "+Benchmark.getTimeSpanHours(start1)+" hours");
    }




//    public static void main(String[] args) {
//        AnchorsVCFForLDPlot.getAnchorPoint(args[0], Integer.parseInt(args[1]), args[2]);
//        AnchorsVCFForLDPlot.getChrPosMaf(args[0], args[1]);
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

