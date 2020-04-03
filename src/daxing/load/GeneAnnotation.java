package daxing.load;

import daxing.common.IOTool;
import daxing.common.RowTableTool;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class GeneAnnotation {

    public static void retainVmapIIGene(String geneAnnotationDB, String outFile){
        RowTableTool<String> table=new RowTableTool<>(geneAnnotationDB);
        Predicate<List<String>> p=l->Integer.parseInt(l.get(2))==0;
        table.removeIf(p);
        table.write(outFile);
    }

    public static void retainTriad(String triadFile, String geneAnnotationFile, String outDir){
        RowTableTool<String> geneAnnotationDB=new RowTableTool<>(geneAnnotationFile);
        Triad triad=new Triad(triadFile);
        Predicate<List<String>> inTriad=l->triad.contain(l.get(0).substring(0, 18));
        geneAnnotationDB.removeIf(inTriad.negate());
        Predicate<List<String>> aSubgenome=l->l.get(0).substring(8,9).equals("A");
        Predicate<List<String>> bSubgenome=l->l.get(0).substring(8,9).equals("B");
        Predicate<List<String>> dSubgenome=l->l.get(0).substring(8,9).equals("D");
        RowTableTool<String> geneAnnotationDB_A=geneAnnotationDB.filter(aSubgenome);
        RowTableTool<String> geneAnnotationDB_B=geneAnnotationDB.filter(bSubgenome);
        RowTableTool<String> geneAnnotationDB_D=geneAnnotationDB.filter(dSubgenome);
        List<String> headerGeneAnnotationDB=geneAnnotationDB.getHeader();
        List<List<String>> data_A=geneAnnotationDB_A.getCells();
        List<List<String>> data_B=geneAnnotationDB_B.getCells();
        List<List<String>> data_D=geneAnnotationDB_D.getCells();
        String[] geneNameTriadABD;
        List<String> column0;
        List<List<String>> newData_A=new ArrayList<>();
        List<List<String>> newData_B=new ArrayList<>();
        List<List<String>> newData_D=new ArrayList<>();
        for (int i = 0; i < data_A.size(); i++) {
            geneNameTriadABD=triad.getTriad(data_A.get(i).get(0).substring(0, 18));
            column0=geneAnnotationDB_B.getColumn(0).stream().map(s->s.substring(0,18)).collect(Collectors.toList());
            int indexB=column0.indexOf(geneNameTriadABD[1]);
            if (indexB < 0) continue;
            column0=geneAnnotationDB_D.getColumn(0).stream().map(s->s.substring(0,18)).collect(Collectors.toList());
            int indexD=column0.indexOf(geneNameTriadABD[2]);
            if (indexD < 0) continue;
            newData_A.add(data_A.get(i));
            newData_B.add(data_B.get(indexB));
            newData_D.add(data_D.get(indexD));
        }
        RowTableTool<String> newRowTableA=new RowTableTool<>(headerGeneAnnotationDB, newData_A);
        RowTableTool<String> newRowTableB=new RowTableTool<>(headerGeneAnnotationDB, newData_B);
        RowTableTool<String> newRowTableD=new RowTableTool<>(headerGeneAnnotationDB, newData_D);
        newRowTableA.write(new File(outDir, "geneAnnotationA.txt").getAbsolutePath());
        newRowTableB.write(new File(outDir, "geneAnnotationB.txt").getAbsolutePath());
        newRowTableD.write(new File(outDir, "geneAnnotationD.txt").getAbsolutePath());
    }

    public static void countHGDeleteriuos(String geneAnnotationA, String geneAnnotationB, String geneAnnotationD,
                                   String outFile, int headerIndex){
        // 000,100,010,001,110,101,011,111
        String[] models={"M000","M100","M010","M001","M110","M101","M011","M111"};
        int[] count=new int[8];
        RowTableTool<String> tableA=new RowTableTool<>(geneAnnotationA);
        RowTableTool<String> tableB=new RowTableTool<>(geneAnnotationB);
        RowTableTool<String> tableD=new RowTableTool<>(geneAnnotationD);
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
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
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
}


