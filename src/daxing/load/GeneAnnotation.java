package daxing.load;

import daxing.applets.PhylipSequential;
import daxing.common.CollectionTool;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class GeneAnnotation {

//    public static void retainVmapIIGeneCDS_AndRetainTriad(String geneAnnotationDB, String triadFile){
//        RowTableTool<String> table=new RowTableTool<>(geneAnnotationDB);
//        String outFile=geneAnnotationDB.replaceAll("\\.txt\\.gz$","_vmap2.1_cds.txt.gz");
//        table.write(outFile, IOFileFormat.TextGzip);
//        retainTriad(triadFile, outFile);
//    }

    public static void retainTriad(String triadFile, String geneAnnotationFile){
        Triad triad=new Triad(triadFile);
        RowTableTool<String> geneAnnotationTable=new RowTableTool<>(geneAnnotationFile);
        List<String> triadGenes=triad.getAllGenes();
        Predicate<List<String>> inTriad=l->triadGenes.contains(l.get(0).substring(0,18));
        geneAnnotationTable.removeIf(inTriad.negate());
        List<String> geneNames=geneAnnotationTable.getColumn(0);
        List<String> genes=geneNames.stream().map(s->s.substring(0,18)).collect(Collectors.toList());
        List<String> triadID=new ArrayList<>(19000);
        List<String> subgenomeList=new ArrayList<>(19000);
        String id, subgenome;
        for (int i = 0; i < genes.size(); i++) {
            id=triad.getTraidID(genes.get(i));
            subgenome=triad.getSubgenome(genes.get(i)).name();
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

    public static void filterNot6(String geneSummaryTriad, String notTriad, String triadFile){
        Triad triad=new Triad(triadFile);
        RowTable<String> notTriaTable=new RowTable<>(notTriad);
        List<String> notTriadGenes=notTriaTable.getColumn(0);
        Set<String> notTriadGeneIDset= notTriadGenes.stream().map(s->triad.getTraidID(s)).collect(Collectors.toSet());
        List<String> notTriadGeneIDList=new ArrayList<>(notTriadGeneIDset);
        RowTableTool<String> geneSummaryTable=new RowTableTool<>(geneSummaryTriad);
        Predicate<List<String>> p=l->notTriadGeneIDList.contains(l.get(0));
        geneSummaryTable.removeIf(p);
        String outFile=geneSummaryTriad.replaceAll("\\.txt\\.gz","_Removed169.txt.gz");
        geneSummaryTable.write(outFile,IOFileFormat.TextGzip);
    }
//    public static void retainTriad(String triadFile, String geneAnnotationFile){
//        RowTableTool<String> geneAnnotationDB=new RowTableTool<>(geneAnnotationFile);
//        Triad triad=new Triad(triadFile);
//        List<String> subgenomeGenes=triad.getAllGenes();
//        Predicate<List<String>> inTriad=l->subgenomeGenes.contains(l.get(0).substring(0, 18));
//        geneAnnotationDB.removeIf(inTriad.negate());
//        Predicate<List<String>> aSubgenome=l->l.get(0).substring(8,9).equals("A");
//        Predicate<List<String>> bSubgenome=l->l.get(0).substring(8,9).equals("B");
//        Predicate<List<String>> dSubgenome=l->l.get(0).substring(8,9).equals("D");
//        RowTableTool<String> geneAnnotationDB_A=geneAnnotationDB.filter(aSubgenome);
//        RowTableTool<String> geneAnnotationDB_B=geneAnnotationDB.filter(bSubgenome);
//        RowTableTool<String> geneAnnotationDB_D=geneAnnotationDB.filter(dSubgenome);
//        List<String> geneA=geneAnnotationDB_A.getColumn(0);
//        List<String> geneB=geneAnnotationDB_B.getColumn(0);
//        List<String> geneD=geneAnnotationDB_D.getColumn(0);
//        String[] geneNameTriadABD;
//        List<String> column0;
//        List<String> geneABD=new ArrayList<>(19000);
//        Set<String> geneAset=new HashSet<>(geneA);
//        geneA=new ArrayList<>(geneAset);
//        for (int i = 0; i < geneA.size(); i++) {
//            geneNameTriadABD=triad.getTriad(geneA.get(i).substring(0, 18));
//            column0=geneB.stream().map(s->s.substring(0,18)).collect(Collectors.toList());
//            int indexB=column0.indexOf(geneNameTriadABD[1]);
//            if (indexB < 0) continue;
//            column0=geneD.stream().map(s->s.substring(0,18)).collect(Collectors.toList());
//            int indexD=column0.indexOf(geneNameTriadABD[2]);
//            if (indexD < 0) continue;
//            geneABD.add(geneA.get(i));
//            geneABD.add(geneB.get(indexB));
//            geneABD.add(geneD.get(indexD));
//        }
//        Predicate<List<String>> p_abd=l->geneABD.contains(l.get(0));
//        geneAnnotationDB.removeIf(p_abd.negate());
//        List<String> geneNames=geneAnnotationDB.getColumn(0);
//        List<String> triadID=new ArrayList<>(19000);
//        String id;
//        for (int i = 0; i < geneNames.size(); i++) {
//            id=triad.getTraidID(geneNames.get(i).substring(0,18));
//            triadID.add(id);
//        }
//        geneAnnotationDB.insertColumn("TriadID",0, triadID);
//        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
//        geneAnnotationDB.sortBy(c);
//        String outFileName=geneAnnotationFile.replaceAll("\\.txt\\.gz$", "_triad.txt.gz");
//        geneAnnotationDB.write(outFileName, IOFileFormat.TextGzip);
//    }

    public static void countHGDeleteriuos(String geneAnnotationABD, String outFile, int headerIndex){
        // 000,100,010,001,110,101,011,111
        String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
        int[] countHGDeleterious=new int[8];
        int[] countSyn=new int[8];
        RowTableTool<String> geneAnnotation=new RowTableTool<>(geneAnnotationABD);
        Predicate<List<String>> pa=l->l.get(1).substring(8,9).equals("A");
        Predicate<List<String>> pb=l->l.get(1).substring(8,9).equals("B");
        Predicate<List<String>> pd=l->l.get(1).substring(8,9).equals("D");
        RowTableTool<String> tableA=geneAnnotation.filter(pa);
        RowTableTool<String> tableB=geneAnnotation.filter(pb);
        RowTableTool<String> tableD=geneAnnotation.filter(pd);
        for (int i = 0; i < tableA.getRowNumber(); i++) {
            int cntA=Integer.parseInt(tableA.getRow(i).get(headerIndex));
            int cntB=Integer.parseInt(tableB.getRow(i).get(headerIndex));
            int cntD=Integer.parseInt(tableD.getRow(i).get(headerIndex));
            if (cntA==0 && cntB==0 && cntD==0){
                countHGDeleterious[0]++;
            }else if (cntB==0 && cntD==0){
                countHGDeleterious[1]++;
            }else if(cntA==0 && cntD==0) {
                countHGDeleterious[2]++;
            }else if (cntA==0 && cntB==0){
                countHGDeleterious[3]++;
            }else if (cntD==0){
                countHGDeleterious[4]++;
            }else if (cntB==0){
                countHGDeleterious[5]++;
            }else if (cntA==0){
                countHGDeleterious[6]++;
            }else{
                countHGDeleterious[7]++;
            }
        }
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            StringBuilder sb;
            bw.write("model\tcount");
            bw.newLine();
            for (int i = 0; i < countHGDeleterious.length; i++) {
                sb=new StringBuilder();
                sb.append(models[i]).append("\t").append(countHGDeleterious[i]);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}


