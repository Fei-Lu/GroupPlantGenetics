package daxing.common;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import pgl.infra.table.RowTable;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class LibraryOfGRT2 {

    String[] libs = null;
    String[][] taxaNames = null;
    String[][] barcodeR1 = null;
    String[][] barcodeR2 = null;
    Table<String, String, String> barcode1_2Taxon;
    String[] libFastqsR1 = null;
    String[] libFastqsR2 = null;
    String cutter1 = null;
    String cutter2 = null;

    public LibraryOfGRT2(String barcodeFileS, String libFastqMapFileS, String cutter1, String cutter2) {
        this.parseBarcode(barcodeFileS, libFastqMapFileS);
        this.cutter1 = cutter1;
        this.cutter2 = cutter2;
    }

    public int getLibNumber () {
        return libs.length;
    }

    public String[] getLibArray () {
        String[] na = new String[libs.length];
        System.arraycopy(libs, 0, na, 0, libs.length);
        return na;
    }

    public String getCutter1 () {
        return this.cutter1;
    }

    public String getCutter2 () {
        return this.cutter2;
    }

    public String getLibName (int index) {
        return libs[index];
    }

    public String[] getTaxaNames (int index) {
        return taxaNames[index];
    }

    public String[] getLibBarcodeR1 (int index) {
        String[] na = new String[barcodeR1[index].length];
        System.arraycopy(barcodeR1[index], 0, na, 0, barcodeR1[index].length);
        return na;
    }

    public String[] getLibBarcodeR2 (int index) {
        String[] na = new String[barcodeR2[index].length];
        System.arraycopy(barcodeR2[index], 0, na, 0, barcodeR2[index].length);
        return na;
    }

    public String getFastqFileSR1 (int index) {
        return libFastqsR1[index];
    }

    public String getFastqFileSR2 (int index) {
        return libFastqsR2[index];
    }

    public void splitBarcode(String outputDir){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        long start=System.nanoTime();
        String[] libs = this.getLibArray();
        for (int i = 0; i < libs.length; i++) {
            String fastqR1 = this.getFastqFileSR1(i);
            String fastqR2 = this.getFastqFileSR2(i);
            Table<String, String, String> bacode1_2Taxon=this.barcode1_2Taxon;
            String[] taxaNames = this.getTaxaNames(i);
            String cutter1 = this.getCutter1();
            String cutter2 = this.getCutter2();
            this.splitBarcode(libs[i], fastqR1, fastqR2, bacode1_2Taxon, taxaNames, outputDir, cutter1, cutter2);
        }
        System.out.println("all libraries completed in "+Benchmark.getTimeSpanHours(start)+" hours");
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    private void splitBarcode(String lib, String fastqR1, String fastqR2, Table<String, String, String> bacode1_2Taxon,
                              String[] taxaNames, String outputDir, String cutter1, String cutter2){
        long start=System.nanoTime();
        HashMap<String, BufferedWriter[]> taxaWriterMap = new HashMap<>();
        BufferedWriter bw1, bw2;
        BufferedWriter[] bwArray;
        for (String str:taxaNames){
            bw1= IOUtils.getTextWriter(new File(outputDir, str+"-1.fq").getAbsolutePath());
            bw2= IOUtils.getTextWriter(new File(outputDir, str+"-2.fq").getAbsolutePath());
            bwArray=new BufferedWriter[2];
            bwArray[0]=bw1;
            bwArray[1]=bw2;
            taxaWriterMap.put(str, bwArray);
        }
        try {
            BufferedReader br1, br2;
            if (fastqR1.endsWith(".gz")) {
                br1 = IOUtils.getTextGzipReader(fastqR1);
            }
            else {
                br1 = IOUtils.getTextReader(fastqR1);
            }
            if (fastqR2.endsWith(".gz")) {
                br2 = IOUtils.getTextGzipReader(fastqR2);
            }
            else {
                br2 = IOUtils.getTextReader(fastqR2);
            }
            Set<String> bSetR1 = bacode1_2Taxon.rowKeySet();
            Set<String> bSetR2 = bacode1_2Taxon.columnKeySet();
            String read1ID, read2ID, seq1, seq2, des1, des2, qual1, qual2;
            String substr1, substr2, barcode1, barcode2;
            int maxBarcodeLen1= bSetR1.stream().map(String::length).mapToInt(Integer::intValue).max().getAsInt();
            int maxBarcodeLen2= bSetR2.stream().map(String::length).mapToInt(Integer::intValue).max().getAsInt();
            int subLen1=maxBarcodeLen1+cutter1.length();
            int subLen2=maxBarcodeLen2+cutter2.length();
            int totalCnt = 0;
            int processedCnt = 0;
            TIntArrayList enzymeIndex1;
            TIntArrayList enzymeIndex2;
            System.out.println("Parsing " + fastqR1 + "\t" + fastqR2);
            StringBuilder sb;
            String taxonName;
            while ((read1ID = br1.readLine()) != null) {
                read2ID = br2.readLine();
                seq1=br1.readLine();
                seq2=br2.readLine();
                des1=br1.readLine();
                des2=br2.readLine();
                qual1=br1.readLine();
                qual2=br2.readLine();
                totalCnt++;
                if (totalCnt%10000000 == 0) {
                    System.out.println("Total read count: "+totalCnt+"\tPassed read count: "+processedCnt);
                }
                substr1=seq1.substring(0, subLen1);
                substr2=seq2.substring(0, subLen2);
                enzymeIndex1= StringTool.getIndexOfSubStr(substr1, "GATCC");
                enzymeIndex2=StringTool.getIndexOfSubStr(substr2, "CGG");
                if (enzymeIndex1.isEmpty()) continue;
                if (enzymeIndex2.isEmpty()) continue;
                barcode1=substr1.substring(0, enzymeIndex1.get(0));
                barcode2=substr2.substring(0, enzymeIndex2.get(0));
                if (!bacode1_2Taxon.containsRow(barcode1)) continue;
                if (!bacode1_2Taxon.containsColumn(barcode2)) continue;
                taxonName=bacode1_2Taxon.get(barcode1, barcode2);
                if (taxonName==null) continue;
                bwArray = taxaWriterMap.get(taxonName);
                sb=new StringBuilder(400);
                sb.append(read1ID).append("\n").append(seq1.substring(barcode1.length())).append("\n").append(des1);
                sb.append("\n").append(qual1.substring(barcode1.length()));
                bwArray[0].write(sb.toString());
                bwArray[0].newLine();
                sb=new StringBuilder(400);
                sb.append(read2ID).append("\n").append(seq2.substring(barcode2.length())).append("\n").append(des2);
                sb.append("\n").append(qual2.substring(barcode2.length()));
                bwArray[1].write(sb.toString());
                bwArray[1].newLine();
                processedCnt++;
            }
            br1.close();
            br2.close();
            for (BufferedWriter[] bufferedWriters: taxaWriterMap.values()){
                bufferedWriters[0].flush();
                bufferedWriters[1].flush();
                bufferedWriters[0].close();
                bufferedWriters[1].close();
            }
            System.out.println("Total read count: "+totalCnt+"\nPassed read count: "+processedCnt);
            System.out.print("The probability of sequencing error in the reads of the "+lib+" library is "+((double)totalCnt-processedCnt)/totalCnt);
            System.out.println(", completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private void parseBarcode (String barcodeFileS, String libFastqMapFileS) {
        RowTable<String> t = new RowTable<>(barcodeFileS);
        Set<String> s = new HashSet<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(t.getCell(i, 1)).append("_").append(t.getCell(i, 2)).append("_").append(t.getCell(i, 3));
            s.add(sb.toString());
        }
        libs = s.toArray(new String[s.size()]);
        Arrays.sort(libs);
        taxaNames = new String[libs.length][];
        barcodeR1 = new String[libs.length][];
        barcodeR2 = new String[libs.length][];
        barcode1_2Taxon= HashBasedTable.create();
        libFastqsR1 = new String[libs.length];
        libFastqsR2 = new String[libs.length];
        for (int i = 0; i < libs.length; i++) {
            List<String> nameList = new ArrayList<>();
            List<String> barcodeR1List = new ArrayList<>();
            List<String> barcodeR2List = new ArrayList<>();
            for (int j = 0; j < t.getRowNumber(); j++) {
                StringBuilder sb = new StringBuilder();
                sb.append(t.getCell(j, 1)).append("_").append(t.getCell(j, 2)).append("_").append(t.getCell(j, 3));
                if (!sb.toString().equals(libs[i])) continue;
                sb = new StringBuilder();
                sb.append(t.getCell(j, 0)).append("_").append(t.getCell(j, 1)).append("_").append(t.getCell(j, 2)).append("_").append(t.getCell(j, 3)).append("_").append(t.getCell(j, 4));
                nameList.add(sb.toString());
                barcodeR1List.add(t.getCell(j, 5));
                barcodeR2List.add(t.getCell(j, 6));
                barcode1_2Taxon.put(t.getCell(j, 5), t.getCell(j, 6), sb.toString());
            }
            taxaNames[i] = nameList.toArray(new String[nameList.size()]);
            barcodeR1[i] = barcodeR1List.toArray(new String[nameList.size()]);
            barcodeR2[i] = barcodeR2List.toArray(new String[nameList.size()]);
        }
        t = new RowTable<>(libFastqMapFileS);
        List<String> lList = new ArrayList();
        for (int i = 0; i < t.getRowNumber(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(t.getCell(i, 0)).append("_").append(t.getCell(i, 1)).append("_").append(t.getCell(i, 2));
            lList.add(sb.toString());
        }
        Collections.sort(lList);
        TIntArrayList availableIndexList = new TIntArrayList();
        for (int i = 0; i < libs.length; i++) {
            int index = Collections.binarySearch(lList, libs[i]);
            if (index < 0) {
                System.out.println(libs[i] + " does not have corresponding fastqs");
            }
            else {
                if (t.getCell(index, 3).equals("NA") || t.getCell(index, 4).equals("NA")) {
                    System.out.println(libs[i] + " does not have corresponding fastqs");
                }
                else {
                    this.libFastqsR1[i] = t.getCell(index, 3);
                    this.libFastqsR2[i] = t.getCell(index, 4);
                    availableIndexList.add(i);
                }
            }
        }
        List<String> libList = new ArrayList();
        List<String[]> taxaNameList = new ArrayList();
        List<String[]> barcodeR1List = new ArrayList();
        List<String[]> barcodeR2List = new ArrayList();
        List<String> libFqR1List = new ArrayList();
        List<String> libFqR2List = new ArrayList();
        int[] aIndex = availableIndexList.toArray();
        for (int i = 0; i < aIndex.length; i++) {
            libList.add(libs[aIndex[i]]);
            taxaNameList.add(taxaNames[aIndex[i]]);
            barcodeR1List.add(barcodeR1[aIndex[i]]);
            barcodeR2List.add(barcodeR1[aIndex[i]]);
            libFqR1List.add(libFastqsR1[aIndex[i]]);
            libFqR2List.add(libFastqsR2[aIndex[i]]);
        }

        libs = libList.toArray(new String[libList.size()]);
        taxaNames = taxaNameList.toArray(new String[taxaNameList.size()][]);
        barcodeR1 = barcodeR1List.toArray(new String[barcodeR1List.size()][]);
        barcodeR2 = barcodeR2List.toArray(new String[barcodeR2List.size()][]);
        libFastqsR1 = libFqR1List.toArray(new String[libFqR1List.size()]);
        libFastqsR2 = libFqR2List.toArray(new String[libFqR2List.size()]);

        System.out.println(libs.length+" libraries will be paralell processd. They are:");
        int cnt = 0;
        for (int i = 0; i < libs.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.append(libs[i]).append(" with ").append(this.taxaNames[i].length).append(" samples");
            System.out.println(sb.toString());
            cnt+=this.taxaNames[i].length;
        }
        System.out.println("A total of " + String.valueOf(cnt) + " samples are in the current batch");
    }
}
