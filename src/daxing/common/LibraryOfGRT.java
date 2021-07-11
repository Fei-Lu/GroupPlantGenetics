package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.StringUtils;
import pgl.app.grt.LibraryInfo;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

/**
 * @author Daxing Xu
 */
public class LibraryOfGRT extends LibraryInfo {

    public LibraryOfGRT(String barcodeFileS, String libFastqMapFileS, String cutter1, String cutter2) {
        super(barcodeFileS, libFastqMapFileS, cutter1, cutter2);
    }

    public void splitBarcode(String outputDir){
        System.out.println(DateTime.getDateTimeOfNow()+" start");
        long start=System.nanoTime();
        String[] libs = this.getLibArray();
        for (int i = 0; i < libs.length; i++) {
            String fastqR1 = this.getFastqFileSR1(i);
            String fastqR2 = this.getFastqFileSR2(i);
            HashMap<String, Set<String>> barcodeR1TaxaMap = this.getbarcodeR1TaxaMap(i);
            HashMap<String, Set<String>> barcodeR2TaxaMap = this.getbarcodeR2TaxaMap(i);
            String[] taxaNames = this.getTaxaNames(i);
            String cutter1 = this.getCutter1();
            String cutter2 = this.getCutter2();
            this.splitBarcode(libs[i], fastqR1, fastqR2, barcodeR1TaxaMap, barcodeR2TaxaMap, taxaNames, outputDir, cutter1, cutter2);
        }
        System.out.println("all libraries completed in "+Benchmark.getTimeSpanHours(start)+" hours");
        System.out.println(DateTime.getDateTimeOfNow()+" end");
    }

    private void splitBarcode(String lib, String fastqR1, String fastqR2,
                              HashMap<String, Set<String>> barcodeR1TaxaMap, HashMap<String, Set<String>> barcodeR2TaxaMap, String[] taxaNames, String outputDir, String cutter1, String cutter2){
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
            Set<String> bSetR1 = barcodeR1TaxaMap.keySet();
            Set<String> bSetR2 = barcodeR2TaxaMap.keySet();
            String read1ID, read2ID, seq1, seq2, des1, des2, qual1, qual2;
            String substr1, substr2, barcode1, barcode2;
            int maxBarcodeLen1= bSetR1.stream().map(String::length).mapToInt(Integer::intValue).max().getAsInt();
            int maxBarcodeLen2= bSetR2.stream().map(String::length).mapToInt(Integer::intValue).max().getAsInt();
            int subLen1=maxBarcodeLen1+cutter1.length();
            int subLen2=maxBarcodeLen2+cutter2.length();
            Set<String> taxaSR1, taxaSR2;
            int totalCnt = 0;
            int processedCnt = 0;
            TIntArrayList enzymeIndex1;
            TIntArrayList enzymeIndex2;
            System.out.println("Parsing " + fastqR1 + "\t" + fastqR2);
            StringBuilder sb;
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
                enzymeIndex1=StringTool.getIndexOfSubStr(substr1, "GATCC");
                enzymeIndex2=StringTool.getIndexOfSubStr(substr2, "CGG");
                if (enzymeIndex1.isEmpty()) continue;
                if (enzymeIndex2.isEmpty()) continue;
                barcode1=substr1.substring(0, enzymeIndex1.get(0));
                barcode2=substr2.substring(0, enzymeIndex2.get(0));
                if (!barcodeR1TaxaMap.containsKey(barcode1)) continue;
                if (!barcodeR2TaxaMap.containsKey(barcode2)) continue;
                taxaSR1=barcodeR1TaxaMap.get(barcode1);
                taxaSR2=barcodeR2TaxaMap.get(barcode2);
                Set<String> newSet=new HashSet<>(taxaSR1);
                newSet.retainAll(taxaSR2);
                if (newSet.size()!=1) continue;
                bwArray = taxaWriterMap.get(newSet.iterator().next());
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



}
