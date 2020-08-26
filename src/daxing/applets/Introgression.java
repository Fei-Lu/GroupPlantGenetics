package daxing.applets;

import daxing.common.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Introgression {

    public static double caculate_D(double[] p1DerivedArray, double[] p2DerivedArray, double[] p3DerivedArray ){
        int len=p1DerivedArray.length;
        if (len!=p2DerivedArray.length || len!= p3DerivedArray.length){
            System.out.println("check parameter, program quit");
            System.exit(1);
        }
        double[] abba=new double[len];
        double[] baba=new double[len];
        for (int i = 0; i < len; i++) {
            abba[i]=(1-p1DerivedArray[i])*(p2DerivedArray[i])*(p3DerivedArray[i]);
            baba[i]=(p1DerivedArray[i])*(1-p2DerivedArray[i])*(p3DerivedArray[i]);
        }
        double sumABBA= Arrays.stream(abba).sum();
        double sumBABA=Arrays.stream(baba).sum();
        return (sumABBA-sumBABA)/(sumABBA+sumBABA);
    }

    public static double caculate_fd(double[] p1DerivedArray, double[] p2DerivedArray, double[] p3aDerivedArray,
                                     double[] p3bDerivedArray){
        int len=p1DerivedArray.length;
        if (len!=p2DerivedArray.length || len != p3aDerivedArray.length || len!=p3bDerivedArray.length){
            System.out.println("check parameter, program quit");
            System.exit(1);
        }
        double[] abba_numerator=new double[len];
        double[] baba_numerator=new double[len];
        double[] abba_denominator=new double[len];
        double[] baba_denominator=new double[len];
        for (int i = 0; i < len; i++) {
            abba_numerator[i]=(1-p1DerivedArray[i])*p2DerivedArray[i]*p3aDerivedArray[i];
            baba_numerator[i]=p1DerivedArray[i]*(1-p2DerivedArray[i])*p3aDerivedArray[i];
            abba_denominator[i]=(1-p1DerivedArray[i])*p3bDerivedArray[i]*p3aDerivedArray[i];
            baba_denominator[i]=p1DerivedArray[i]*(1-p3bDerivedArray[i])*p3aDerivedArray[i];
        }
        double sum_abba_numerator= Arrays.stream(abba_numerator).sum();
        double sum_baba_numerator=Arrays.stream(baba_numerator).sum();
        double sum_abba_denominator=Arrays.stream(abba_denominator).sum();
        double sum_baba_denominator=Arrays.stream(baba_denominator).sum();
        return (sum_abba_numerator-sum_baba_numerator)/(sum_abba_denominator-sum_baba_denominator);
    }

    /**
     * split p3 to p3a and p3b
     * @param treeByGroupFile
     * @param treeByGroup_p3ab
     */
    public static void regroup(String treeByGroupFile, String treeByGroup_p3ab){
        String line;
        List<String> temp;
        String[] p3Array={"WildEmmer","DomesticatedEmmer","FreeThreshTetraploid","Ae.tauschii"};
        List<String> p3List= Arrays.stream(p3Array).collect(Collectors.toList());
        try (BufferedReader bufferedReader = IOTool.getReader(treeByGroupFile);
             BufferedWriter bufferedWriter=IOTool.getWriter(treeByGroup_p3ab)) {
            String header=bufferedReader.readLine();
            bufferedWriter.write(header);
            bufferedWriter.newLine();
            StringBuilder sb;
            boolean flag=true;
            while ((line=bufferedReader.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                sb=new StringBuilder();
                if (p3List.contains(temp.get(1))){
                    if (flag){
                        sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("_a");
                        bufferedWriter.write(sb.toString());
                        bufferedWriter.newLine();
                        flag=false;
                    }else {
                        sb.append(temp.get(0)).append("\t").append(temp.get(1)).append("_b");
                        bufferedWriter.write(sb.toString());
                        bufferedWriter.newLine();
                        flag=true;
                    }
                }else {
                    bufferedWriter.write(line);
                    bufferedWriter.newLine();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void addGeneForFreq(String freqTableDir, String outDir, String hcGeneFile, String pgfFile,
                                      String chrConvertionRule){
        List<File> files= IOUtils.getVisibleFileListInDir(freqTableDir);
        String[] outNames= files.stream().map(File::getName).map(str->str.replaceAll("\\.freq","_HC_gene.freq"))
                            .toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->addGeneForFreq(files.get(e), new File(outDir, outNames[e]),
                hcGeneFile, pgfFile,chrConvertionRule));
    }

    /**
     * add gene name on derived allele freq table
     * @param freqTable
     * @param outFile
     * @param hcGeneFile
     * @param pgfFile
     * @param chrConvertionRule
     */
    public static void addGeneForFreq(File freqTable, File outFile, String hcGeneFile, String pgfFile,
                                      String chrConvertionRule){
        RowTableTool<String> hcGeneTable=new RowTableTool<>(hcGeneFile);
        List<String> hcGeneList=hcGeneTable.getColumn(0);
        PGF pgf=new PGF(pgfFile);
        Predicate<PGF.Gene> hcGenePredict2= gene -> hcGeneList.contains(gene.getGeneName());
        pgf.removeIf(hcGenePredict2.negate());
        try (BufferedReader bufferedReader = IOTool.getReader(freqTable);
             BufferedWriter bufferedWriter =IOTool.getWriter(outFile)) {
            String header=bufferedReader.readLine();
            bufferedWriter.write(header);
            bufferedWriter.newLine();
            String line;
            List<String> temp;
            String refChr;
            int refPos;
            int chr,pos;
            ChrConvertionRule c=new ChrConvertionRule(chrConvertionRule);
            int geneIndex;
            while ((line=bufferedReader.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                refChr=temp.get(0);
                refPos=Integer.parseInt(temp.get(1));
                chr=c.getVCFChrFromRefChrPos(refChr,refPos);
                pos=c.getVCFPosFromRefChrPos(refChr,refPos);
                geneIndex=pgf.getGeneIndex(chr,pos);
                if (geneIndex<0){
                    System.out.println(refChr+"("+chr+")"+"\t"+refPos+"("+pos+")"+" is not not in the range of any " +
                            "gene");
                    continue;
                }
                temp.add(pgf.getGeneName(geneIndex));
                bufferedWriter.write(String.join("\t", temp));
                bufferedWriter.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public enum P2{
        Cultivar, Landrace_Europe, Landrace_WestAsia, Landrace_EastAsia
    }

    public enum P3_subgenomeAB{
        WildEmmer, DomesticatedEmmer, FreeThreshTetraploid
    }

    public enum P3_subgenomeD{
        tauschii  //Ae.tauschii
    }

    private static String getAB_Header(){
        StringBuilder sb=new StringBuilder();
        String[] d_fd={"D","fd"};
        for(P2 p2: P2.values()){
            for(P3_subgenomeAB p3_subgenomeAB: P3_subgenomeAB.values()){
                for (int i = 0; i < d_fd.length; i++) {
                    sb.append(d_fd[i]).append("-").append(p2).append("-").append(p3_subgenomeAB).append("\t");
                }
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    private static String getD_Header(){
        StringBuilder sb=new StringBuilder();
        String[] d_fd={"D","fd"};
        for(P2 p2: P2.values()){
            for (int i = 0; i < d_fd.length; i++) {
                sb.append(d_fd[i]).append("-").append(p2).append("-").append(P3_subgenomeD.values()[0]).append("\t");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    /**
     *
     * @param geneDerivedFreqFile
     * @param outFile
     * @param windowSize
     */
    public static void slidingWindowAB(File geneDerivedFreqFile, File outFile, int windowSize){
        long startTime=System.nanoTime();
        Set<String> geneNames=RowTableTool.getColumnSet(geneDerivedFreqFile.getAbsolutePath(), 14);
        List<String> genesList=new ArrayList<>(geneNames);
        Collections.sort(genesList);
        ColumnTableTool<String> columnTable=new ColumnTableTool<>(geneDerivedFreqFile.getAbsolutePath());
        List<String> chrposGeneName=columnTable.getColumn(14);
        List<GeneWindow> geneWindowList=new ArrayList<>();
        for (int i = 0; i < genesList.size(); i++) {
            int startIndex=chrposGeneName.indexOf(genesList.get(i));
            int endIndex=chrposGeneName.lastIndexOf(genesList.get(i));
            int increment=Integer.MIN_VALUE;
            if ((endIndex-startIndex) > windowSize){
                increment=(endIndex-startIndex-windowSize)/2;
                geneWindowList.add(new GeneWindow(genesList.get(i),startIndex+increment, endIndex-increment));
            }else {
                increment = (windowSize - (endIndex - startIndex)) / 2;
                if ((startIndex - increment) < 0) {
                    startIndex = 0 + increment;
                }
                if ((endIndex + increment) >= columnTable.getRowNumber()) {
                    endIndex = columnTable.getRowNumber() - 1 - increment;
                }
                geneWindowList.add(new GeneWindow(genesList.get(i), startIndex - increment, endIndex + increment));
            }
        }
        TDoubleArrayList p1DerivedList, p2DerivedList, p3aDerivedList, p3bDerivedList;
        double[] p1, p2, p3a, p3b;
        int start=Integer.MIN_VALUE;
        int end=Integer.MIN_VALUE;
        String startPos, endPos;
        int middleIndex=Integer.MIN_VALUE;
        String geneMiddlePos=null;
        String chr=null;
        try (BufferedWriter bufferedWriter = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("chr\tgeneName\tstartPos\tendPos\tgeneMiddlePos\tsiteInWindow\t").append(getAB_Header());
            bufferedWriter.write(sb.toString());
            bufferedWriter.newLine();
            for (int i = 0; i < geneWindowList.size(); i++) {
                sb=new StringBuilder(1000);
                start=geneWindowList.get(i).getStart();
                end=geneWindowList.get(i).getEnd();
                startPos=columnTable.getColumn(1).get(start);
                endPos=columnTable.getColumn(1).get(end);
                middleIndex=geneWindowList.get(i).getMid();
                chr=geneWindowList.get(i).getChr();
                geneMiddlePos=columnTable.getColumn(1).get(middleIndex);
                p1DerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray("IndianDwarfWheat"));
                p1=p1DerivedList.subList(start, end).toArray();
                sb.append(chr).append("\t").append(geneWindowList.get(i).getGeneName()).append("\t");
                sb.append(startPos).append("\t").append(endPos).append("\t").append(geneMiddlePos).append("\t");
                sb.append(geneWindowList.get(i).getRowNum()).append("\t");
                for (int j = 0; j < P2.values().length; j++) {
                    p2DerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray(P2.values()[j].name()));
                    p2=p2DerivedList.subList(start,end).toArray();
                    for (int k = 0; k < P3_subgenomeAB.values().length; k++) {
                        p3aDerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray(P3_subgenomeAB.values()[k].name()+"_a"));
                        p3bDerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray(P3_subgenomeAB.values()[k].name()+"_b"));
                        p3a=p3aDerivedList.subList(start,end).toArray();
                        p3b=p3bDerivedList.subList(start,end).toArray();
                        double fd=caculate_fd(p1, p2, p3a, p3b);
                        double dValue=caculate_D(p1, p2, p3a);
                        sb.append(dValue).append("\t").append(fd).append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(geneDerivedFreqFile.getName()+" completed in "+ Benchmark.getTimeSpanMinutes(startTime)+" " +
                "minutes");
    }

    public static void slidingWindowD(File geneDerivedFreqFile, File outFile, int windowSize){
        long startTime=System.nanoTime();
        Set<String> geneNames=RowTableTool.getColumnSet(geneDerivedFreqFile.getAbsolutePath(), 10);
        List<String> genesList=new ArrayList<>(geneNames);
        Collections.sort(genesList);
        ColumnTableTool<String> columnTable=new ColumnTableTool<>(geneDerivedFreqFile.getAbsolutePath());
        List<String> chrposGeneName=columnTable.getColumn(10);
        List<GeneWindow> geneWindowList=new ArrayList<>();
        for (int i = 0; i < genesList.size(); i++) {
            int startIndex=chrposGeneName.indexOf(genesList.get(i));
            int endIndex=chrposGeneName.lastIndexOf(genesList.get(i));
            int increment=Integer.MIN_VALUE;
            if ((endIndex-startIndex) > windowSize){
                increment=(endIndex-startIndex-windowSize)/2;
                geneWindowList.add(new GeneWindow(genesList.get(i),startIndex+increment, endIndex-increment));
            }else {
                increment = (windowSize - (endIndex - startIndex)) / 2;
                if ((startIndex - increment) < 0) {
                    startIndex = 0 + increment;
                }
                if ((endIndex + increment) > columnTable.getRowNumber()) {
                    endIndex = columnTable.getRowNumber() - 1 - increment;
                }
                geneWindowList.add(new GeneWindow(genesList.get(i), startIndex - increment, endIndex + increment));
            }
        }
        TDoubleArrayList p1DerivedList, p2DerivedList, p3aDerivedList, p3bDerivedList;
        double[] p1, p2, p3a, p3b;
        int start=Integer.MIN_VALUE;
        int end=Integer.MIN_VALUE;
        String startPos, endPos;
        int middleIndex=Integer.MIN_VALUE;
        String geneMiddlePos=null;
        String chr=null;
        try (BufferedWriter bufferedWriter = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            sb.append("chr\tgeneName\tstartPos\tendPos\tgeneMiddlePos\tsiteInWindow\t").append(getD_Header());
            bufferedWriter.write(sb.toString());
            bufferedWriter.newLine();
            for (int i = 0; i < geneWindowList.size(); i++) {
                sb=new StringBuilder(1000);
                start=geneWindowList.get(i).getStart();
                end=geneWindowList.get(i).getEnd();
                startPos=columnTable.getColumn(1).get(start);
                endPos=columnTable.getColumn(1).get(end);
                middleIndex=geneWindowList.get(i).getMid();
                chr=geneWindowList.get(i).getChr();
                geneMiddlePos=columnTable.getColumn(1).get(middleIndex);
                p1DerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray("IndianDwarfWheat"));
                p1=p1DerivedList.subList(start, end).toArray();
                sb.append(chr).append("\t").append(geneWindowList.get(i).getGeneName()).append("\t");
                sb.append(startPos).append("\t").append(endPos).append("\t").append(geneMiddlePos).append("\t");
                sb.append(geneWindowList.get(i).getRowNum()).append("\t");
                for (int j = 0; j < P2.values().length; j++) {
                    p2DerivedList=new TDoubleArrayList(columnTable.getColumnAsDoubleArray(P2.values()[j].name()));
                    p2=p2DerivedList.subList(start,end).toArray();
                    p3aDerivedList= new TDoubleArrayList(columnTable.getColumnAsDoubleArray("Ae.tauschii_a"));
                    p3bDerivedList= new TDoubleArrayList(columnTable.getColumnAsDoubleArray("Ae.tauschii_b"));
                    p3a=p3aDerivedList.subList(start,end).toArray();
                    p3b=p3bDerivedList.subList(start,end).toArray();
                    double fd=caculate_fd(p1, p2, p3a, p3b);
                    double dValue=caculate_D(p1, p2, p3a);
                    sb.append(dValue).append("\t").append(fd).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bufferedWriter.write(sb.toString());
                bufferedWriter.newLine();
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(geneDerivedFreqFile.getName()+" completed in "+ Benchmark.getTimeSpanMinutes(startTime)+" minutes");
    }

    /**
     *
     * @param geneDerivedFreqFileDir
     * @param outDir
     * @param windowSize
     */
    public static void slidingWindowForfd(String geneDerivedFreqFileDir, String outDir, int windowSize){
        long start=System.nanoTime();
        List<File> files= IOUtils.getVisibleFileListInDir(geneDerivedFreqFileDir);
        Predicate<File> d=file -> file.getName().substring(4,5).equals("D");
        List<File> abFiles=files.stream().filter(d.negate()).collect(Collectors.toList());
        List<File> dFiles=files.stream().filter(d).collect(Collectors.toList());
        String[] abOutNames= abFiles.stream().map(File::getName)
                .map(str->str.replaceAll("freq\\.tsv","D.fd.txt")).toArray(String[]::new);
        String[] dOutNames=dFiles.stream().map(File::getName)
                .map(str->str.replaceAll("freq\\.tsv","D.fd.txt")).toArray(String[]::new);
        IntStream.range(0, abFiles.size()).parallel().forEach(e->slidingWindowAB(abFiles.get(e), new File(outDir,
                abOutNames[e]), windowSize));
        IntStream.range(0, dFiles.size()).parallel().forEach(e->slidingWindowD(dFiles.get(e), new File(outDir, dOutNames[e]),
                windowSize));
        System.out.println("completed in " +Benchmark.getTimeSpanMinutes(start)+" minutes");
    }

    /**
     *
     * @param derivedFreqDir
     * @param outDir
     */
    public static void removedDerivedAlleleFreqRowsContainsNan(String derivedFreqDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(derivedFreqDir);
        String[] outNames=
                files.stream().map(File::getName).map(str->str.replaceAll("\\.exon","_removedNan_exon")).toArray(String[]::new);
        Predicate<List<String>> nan=l->l.contains("nan");
        RowTableTool<String> rowTableTool;
        for (int i = 0; i < files.size(); i++) {
            rowTableTool=new RowTableTool<>(files.get(i).getAbsolutePath());
            rowTableTool.removeIf(nan);
            rowTableTool.write(new File(outDir, outNames[i]).getAbsolutePath());
        }
    }

    public static void splitP2P3(String removedNanHCgenefdDir, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(removedNanHCgenefdDir);
        String chr;
        ColumnTableTool<String> columnTableTool;
        List<String> header;
        ColumnTableTool columns;
        String outFile;
        for (int i = 0; i < files.size(); i++) {
            chr=files.get(i).getName().substring(0,5);
            columnTableTool=new ColumnTableTool<>(files.get(i).getAbsolutePath());
            header=columnTableTool.getHeader();
            TIntArrayList index=new TIntArrayList();
            index.add(IntStream.range(0,6).toArray());
            for (int j = 6; j < header.size(); j=j+2) {
                index.add(j);
                index.add(j+1);
                columns=columnTableTool.selectColumn(index.toArray());
                outFile=new File(outDir, chr+"_"+header.get(j).substring(2)+"_geneWindow.txt").getAbsolutePath();
                columns.writeTextTable(outFile, IOFileFormat.Text);
                index.remove(j);
                index.remove(j+1);
            }
        }
    }

//    public static void main(String[] args) {
////        regroup("","");
//        //calculate derived allele freq per chromosome and removed the rows which contains nan
////        removedDerivedAlleleFreqRowsContainsNan("", "");
//        // add gene name and retain HCgene
////        addGeneForFreq("","","","","");
//        slidingWindowForfd("","",100);
//    }
}
