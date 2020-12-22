package daxing.applets;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.Triads;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class GeneSummaryUtils {

    public static void retainTriad(String triadFile, String geneAnnotationFile){
        Triads triads =new Triads(triadFile);
        RowTableTool<String> geneAnnotationTable=new RowTableTool<>(geneAnnotationFile);
        List<String> triadGenes= triads.getAllGenes();
        Predicate<List<String>> inTriad=l->triadGenes.contains(l.get(0).substring(0,18));
        geneAnnotationTable.removeIf(inTriad.negate());
        List<String> geneNames=geneAnnotationTable.getColumn(0);
        List<String> genes=geneNames.stream().map(s->s.substring(0,18)).collect(Collectors.toList());
        List<String> triadID=new ArrayList<>(19000);
        List<String> subgenomeList=new ArrayList<>(19000);
        String id, subgenome;
        for (int i = 0; i < genes.size(); i++) {
            id= triads.getTraidID(genes.get(i));
            subgenome= triads.getSubgenome(genes.get(i)).name();
            triadID.add(id);
            subgenomeList.add(subgenome);
        }
        geneAnnotationTable.insertColumn("TriadID",0,triadID);
        geneAnnotationTable.insertColumn("sub", 28, subgenomeList);
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        geneAnnotationTable.sortBy(c);
        String outFileName=geneAnnotationFile.replaceAll("\\.txt\\.gz$", "_triad.txt.gz");
        geneAnnotationTable.write(outFileName, IOFileFormat.TextGzip);
    }

    public static void remove169(String geneSummaryTriad, String notTriad, String triadFile){
        Triads triads =new Triads(triadFile);
        RowTable<String> notTriaTable=new RowTable<>(notTriad);
        List<String> notTriadGenes=notTriaTable.getColumn(0);
        Set<String> notTriadGeneIDset= notTriadGenes.stream().map(s-> triads.getTraidID(s)).collect(Collectors.toSet());
        List<String> notTriadGeneIDList=new ArrayList<>(notTriadGeneIDset);
        RowTableTool<String> geneSummaryTable=new RowTableTool<>(geneSummaryTriad);
        Predicate<List<String>> p=l->notTriadGeneIDList.contains(l.get(0));
        geneSummaryTable.removeIf(p);
        String outFile=geneSummaryTriad.replaceAll("\\.txt\\.gz","_Removed169.txt.gz");
        geneSummaryTable.write(outFile,IOFileFormat.TextGzip);
    }

    public static void countHGDeleteriuos(String geneAnnotationABD, String outFile, int headerIndex){
        // 000,100,010,001,110,101,011,111
        String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
        int[] count=new int[8];
        RowTableTool<String> geneAnnotation=new RowTableTool<>(geneAnnotationABD);
        Predicate<List<String>> pa=l->l.get(1).substring(8,9).equals("A");
        Predicate<List<String>> pb=l->l.get(1).substring(8,9).equals("B");
        Predicate<List<String>> pd=l->l.get(1).substring(8,9).equals("D");
        RowTableTool<String> tableA=geneAnnotation.filter(pa);
        RowTableTool<String> tableB=geneAnnotation.filter(pb);
        RowTableTool<String> tableD=geneAnnotation.filter(pd);
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        tableA.sortBy(c);
        tableB.sortBy(c);
        tableD.sortBy(c);
        for (int i = 0; i < tableA.getRowNumber(); i++) {
            int cntA=Integer.parseInt(tableA.getRow(i).get(headerIndex));
            int cntB=Integer.parseInt(tableB.getRow(i).get(headerIndex));
            int cntD=Integer.parseInt(tableD.getRow(i).get(headerIndex));
            if (cntA==0 && cntB==0 && cntD==0){
                count[0]++;
            }else if (cntB==0 && cntD==0){
                count[1]++;
            }else if(cntA==0 && cntD==0) {
                count[2]++;
            }else if (cntA==0 && cntB==0){
                count[3]++;
            }else if (cntD==0){
                count[4]++;
            }else if (cntB==0){
                count[5]++;
            }else if (cntA==0){
                count[6]++;
            }else{
                count[7]++;
            }
        }
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb;
            bw.write("model\tcount");
            bw.newLine();
            for (int i = 0; i < count.length; i++) {
                sb=new StringBuilder();
                sb.append(models[i]).append("\t").append(count[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    private static int[] getCount(RowTableTool<String> geneAnnotationABD, int headerIndex){
        // 000,100,010,001,110,101,011,111
//        String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
        int[] count=new int[8];
        Predicate<List<String>> pa=l->l.get(1).substring(8,9).equals("A");
        Predicate<List<String>> pb=l->l.get(1).substring(8,9).equals("B");
        Predicate<List<String>> pd=l->l.get(1).substring(8,9).equals("D");
        RowTableTool<String> tableA=geneAnnotationABD.filter(pa);
        RowTableTool<String> tableB=geneAnnotationABD.filter(pb);
        RowTableTool<String> tableD=geneAnnotationABD.filter(pd);
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        tableA.sortBy(c);
        tableB.sortBy(c);
        tableD.sortBy(c);
        for (int i = 0; i < tableA.getRowNumber(); i++) {
            int cntA=Integer.parseInt(tableA.getRow(i).get(headerIndex));
            int cntB=Integer.parseInt(tableB.getRow(i).get(headerIndex));
            int cntD=Integer.parseInt(tableD.getRow(i).get(headerIndex));
            if (cntA==0 && cntB==0 && cntD==0){
                count[0]++;
            }else if (cntB==0 && cntD==0){
                count[1]++;
            }else if(cntA==0 && cntD==0) {
                count[2]++;
            }else if (cntA==0 && cntB==0){
                count[3]++;
            }else if (cntD==0){
                count[4]++;
            }else if (cntB==0){
                count[5]++;
            }else if (cntA==0){
                count[6]++;
            }else{
                count[7]++;
            }
        }
        return count;
    }

    public static void getCount(String geneAnnotationABD, String outFile){
        RowTableTool<String> geneAnnotation=new RowTableTool<>(geneAnnotationABD);
        Set<String> groupset=new HashSet<>(geneAnnotation.getColumn(28));
        List<String> groups=new ArrayList<>(groupset);
        RowTableTool<String>[] tables=new RowTableTool[groups.size()];
        for (int i = 0; i < groups.size(); i++) {
            int temp=i;
            tables[i]=geneAnnotation.filter(l->l.get(28).equals(groups.get(temp)));
        }
        int[][] delCount=new int[groups.size()][];
        int[][] nonsynCount=new int[groups.size()][];
        int[][] synCount=new int[groups.size()][];
        for (int i = 0; i < groups.size(); i++) {
            delCount[i]=getCount(tables[i],13);
            nonsynCount[i]=getCount(tables[i], 8);
            synCount[i]=getCount(tables[i], 6);
        }
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("model\tdelCount\tnonsynCount\tsynCount\tgroup");
            bw.newLine();
            String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
            StringBuilder sb;
            for (int i = 0; i < groups.size(); i++) {
                for (int j = 0; j < delCount[i].length; j++) {
                    sb=new StringBuilder();
                    sb.append(models[j]).append("\t");
                    sb.append(delCount[i][j]).append("\t").append(nonsynCount[i][j]).append("\t");
                    sb.append(synCount[i][j]).append("\t").append(groups.get(i));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

//    public static void main(String[] args) {
//        retainTriad("","geneSummaryFile");
//        remove169("", "notTriadFile","");
//        getCount("triad_removed169File","");
//    }
}


