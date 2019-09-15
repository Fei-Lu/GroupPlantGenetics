package daxing.informal;

import format.table.ColumnTable;
import utils.IOUtils;
import utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

import static java.util.stream.Collectors.counting;
import static java.util.stream.Collectors.groupingBy;

public class Erbao {


    public Erbao(String dbFile, String querydir, String outDir){
        this.getRes(dbFile,querydir, outDir);
    }

    public void getRes(String dbInputFile, String queryDir, String outdir){
        File[] files=IOUtils.listRecursiveFiles(new File(queryDir));
        Predicate<File> p=File::isHidden;
        files=Arrays.stream(files).filter(p.negate()).toArray(File[]::new);
        Arrays.sort(files);
        String[] outNames=Arrays.stream(files).map(File::getName).map(str->str.replaceAll(".txt$", "")).toArray(String[]::new);
        for (int i = 0; i < files.length; i++) {
            this.extractSampleFromDB(dbInputFile, files[i].getAbsolutePath(), new File(outdir, outNames[i]).getAbsolutePath());
        }
    }

    public  void extractSampleFromDB(String db, String queryFile, String outDir){
        String line;
        String header;
        File file=new File(outDir);
        file.mkdir();
        List<String> haptypes=new ArrayList<>();
        List<List<String>> querySamples=new ArrayList<>();
        List<String> lineList;
        try {
            List<String> accessions;
            BufferedReader brQuery= IOUtils.getTextReader(queryFile);
            brQuery.readLine();
            while ((line=brQuery.readLine())!=null){
                lineList= PStringUtils.fastSplit(line);
                haptypes.add(lineList.get(0));
                accessions=PStringUtils.fastSplit(lineList.get(3), ",");
                querySamples.add(accessions);
            }
            Collections.sort(haptypes);
            brQuery.close();
            BufferedReader brDB=IOUtils.getTextReader(db);
            for (int i = 0; i < querySamples.size(); i++) {
                Collections.sort(querySamples.get(i));
            }
            header=brDB.readLine();
            BufferedWriter[] bws =new BufferedWriter[haptypes.size()];
            for (int i = 0; i < bws.length; i++) {
                bws[i]=IOUtils.getTextWriter(new File(outDir, haptypes.get(i)+".txt").getAbsolutePath());
                bws[i].write("Haplotype"+"\t"+header);
                bws[i].newLine();
            }
            int[] indexIndex;
            List<String> lines=new ArrayList<>();
            String temp;
            StringBuilder sb=new StringBuilder();
            while ((line=brDB.readLine())!=null){
                lineList=PStringUtils.fastSplit(line);
                indexIndex=this.binarySearch(querySamples, lineList.get(0));
                if (indexIndex[1] > -1){
                    sb.append(haptypes.get(indexIndex[0])).append("\t").append(line);
                    lines.add(sb.toString());
                    sb=new StringBuilder();
                }
            }
            brDB.close();
            Comparator<String> comparator=Comparator.comparing(e->PStringUtils.fastSplit(e).get(0));
            comparator.thenComparing(e->PStringUtils.fastSplit(e).get(1));
            Collections.sort(lines, comparator);
            Map<String, List<String>> haptype=lines.stream().collect(Collectors.groupingBy(e->PStringUtils.fastSplit(e).get(0)));
            List<String> key=new ArrayList<>(haptype.keySet());
            Collections.sort(key);
            for (int i = 0; i < key.size(); i++) {
                for (int j = 0; j < haptype.get(key.get(i)).size(); j++) {
                    bws[i].write(haptype.get(key.get(i)).get(j));
                    bws[i].newLine();
                }
                bws[i].flush();
                bws[i].close();
            }
            this.caculateAverage(outDir, new File(outDir, "res.txt").getAbsolutePath());
        }catch (Exception e){
            e.printStackTrace();
        }
    }


    public int[] binarySearch(List<List<String>> querySamples, String sample){
        int index=Integer.MIN_VALUE;
        int[] indexIndex=new int[2];
        indexIndex[0]= Integer.MIN_VALUE;
        indexIndex[1]= Integer.MIN_VALUE;
        for (int i = 0; i < querySamples.size(); i++) {
            index=Collections.binarySearch(querySamples.get(i), sample);
            if (index > -1){
                indexIndex[0]=i;
                indexIndex[1]=index;
            }
        }
        return indexIndex;
    }

    public void caculateAverage(String inputDir, String outFile){
        File[] files=new File(inputDir).listFiles();
        Arrays.sort(files);
        List<ColumnTable<String>> columnTables=new ArrayList<>();
        for (int i = 0; i < files.length; i++) {
            columnTables.add(new ColumnTable<>(files[i].getAbsolutePath()));
        }
        try{
            BufferedWriter bw=IOUtils.getTextWriter(outFile);
            StringBuilder sb=new StringBuilder();
            sb.append("Haplotype\tsubpopulation\tPB-HN\tPB-LN\tEPN-HN\tEPN-LN\tCUDU-HN\tCUDU-LN\tPH-HN\tPH-LN\tHNGNP\tHNFGN\tLNGNP\tLNFGN\tPL-HN\tPL-LN\tCHANGDU-HN\tCHANGDU-LN\tLV-HN\tLV-LN\n");
            bw.write(sb.toString());
            List<String> column;
            Map<Object, Long> map;
            double[] columnNum;
            double mean;
            Predicate<String> p= e->e.equals("NA") || e.equals("#VALUE!");
            for (int i = 0; i < columnTables.size(); i++) {
                sb=new StringBuilder();
                column=columnTables.get(i).getColumn(2);
                map=column.stream().collect(groupingBy(str->PStringUtils.fastSplit(str).get(0),counting()));
                sb.append(files[i].getName()).append("\t");
                for (Map.Entry<Object, Long> entry: map.entrySet()){
                    sb.append(entry.getKey()).append(":").append(entry.getValue()).append(", ");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (int j = 3; j < columnTables.get(i).getColumnNumber(); j++) {
                    column=columnTables.get(i).getColumn(j);
                    column=column.stream().filter(p.negate()).collect(Collectors.toList());
                    columnNum=column.stream().map(Double::parseDouble).mapToDouble(Double::doubleValue).toArray();
                    mean=Arrays.stream(columnNum).summaryStatistics().getAverage();
                    sb.append(mean).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }catch (Exception e){
            e.printStackTrace();
        }




    }

//    public static void main(String[] args) {
//        new Erbao(args[0], args[1], args[2]);
//    }
}
