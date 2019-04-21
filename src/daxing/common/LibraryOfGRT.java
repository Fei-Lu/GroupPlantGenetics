package daxing.common;

import analysis.pipeline.grt.LibraryInfo;
import format.dna.Read;
import utils.IOUtils;
import utils.PArrayUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.*;

public class LibraryOfGRT extends LibraryInfo {
    int numThreads=32;

    public LibraryOfGRT(String barcodeFileS, String libFastqMapFileS, String cutter1, String cutter2) {
        super(barcodeFileS, libFastqMapFileS, cutter1, cutter2);
    }

    public void splitBarcode(String outputDir){
        String[] libs = this.getLibArray();
        int[][] indices = PArrayUtils.getSubsetsIndicesBySubsetSize(libs.length, this.numThreads);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> indexList = Arrays.asList(subLibIndices);
            indexList.parallelStream().forEach(index -> {
                String fastqR1 = this.getFastqFileSR1(index);
                String fastqR2 = this.getFastqFileSR2(index);
                HashMap<String, Set<String>> barcodeR1TaxaMap = this.getbarcodeR1TaxaMap(index);
                HashMap<String, Set<String>> barcodeR2TaxaMap = this.getbarcodeR2TaxaMap(index);
                String[] taxaNames = this.getTaxaNames(index);
                String cutter1 = this.getCutter1();
                String cutter2 = this.getCutter2();
                this.splitBarcode(libs[index], fastqR1, fastqR2, barcodeR1TaxaMap, barcodeR2TaxaMap, taxaNames, outputDir, cutter1, cutter2);
            });
        }
    }

    private void splitBarcode(String lib, String fastqR1, String fastqR2,
                              HashMap<String, Set<String>> barcodeR1TaxaMap, HashMap<String, Set<String>> barcodeR2TaxaMap, String[] taxaNames, String outputDir, String cutter1, String cutter2){
        HashMap<String, BufferedWriter[]> taxaWriterMap = new HashMap<>();
        BufferedWriter bw1, bw2;
        BufferedWriter[] bwArray;
        File outputSubDir;
        for (String str:taxaNames){
            outputSubDir=new File(outputDir, str);
            outputSubDir.mkdir();
            bw1= IOUtils.getNIOTextWriter(new File(outputSubDir, str+"-1.fq").getAbsolutePath());
            bw2= IOUtils.getNIOTextWriter(new File(outputSubDir, str+"-2.fq").getAbsolutePath());
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
            String[] barcodeR1 = bSetR1.toArray(new String[bSetR1.size()]);
            Set<String> bSetR2 = barcodeR2TaxaMap.keySet();
            String[] barcodeR2 = bSetR2.toArray(new String[bSetR2.size()]);
            Arrays.sort(barcodeR1);
            Arrays.sort(barcodeR2);
            String read1ID, read2ID, seq1, seq2;
            int index1, index2;
            Set<String> taxaSR1, taxaSR2;
            int totalCnt = 0;
            int processedCnt = 0;
            int lenBarcodeR1, lenBarcodeR2;
            System.out.println("Parsing " + fastqR1 + "\t" + fastqR2);
            Read[] reads1_2;
            boolean b;
            BufferedWriter bwError1=IOUtils.getNIOTextWriter("/Users/xudaxing/Desktop/test/sequencingError-1.fq");
            BufferedWriter bwError2=IOUtils.getNIOTextWriter("/Users/xudaxing/Desktop/test/sequencingError-2.fq");
            while ((read1ID = br1.readLine()) != null) {
                read2ID = br2.readLine();
                seq1=br1.readLine();
                seq2=br2.readLine();
                totalCnt++;
                if (totalCnt%10000000 == 0) {
                    System.out.println("Total read count: "+totalCnt+"\tPassed read count: "+processedCnt);
                }
                index1 = Arrays.binarySearch(barcodeR1, seq1);
                index2 = Arrays.binarySearch(barcodeR2, seq2);
                b=(index1 == -1 || index2==-1 );
                if (b){
                    bwError1.write(String.join("\n", read1ID, seq1, br1.readLine(), br1.readLine()));
                    bwError1.newLine();
                    bwError2.write(String.join("\n", read2ID, seq2, br2.readLine(), br2.readLine()));
                    bwError2.newLine();
                    continue;
                }
                index1 = -index1 - 2;
                index2 = -index2 - 2;
                taxaSR1 = barcodeR1TaxaMap.get(barcodeR1[index1]);
                taxaSR2 = barcodeR2TaxaMap.get(barcodeR2[index2]);
                Set<String> newSet = new HashSet<>(taxaSR1);
                newSet.retainAll(taxaSR2);
                if (newSet.size() != 1) {
                    bwError1.write(String.join("\n", read1ID, seq1, br1.readLine(), br1.readLine()));
                    bwError1.newLine();
                    bwError2.write(String.join("\n", read2ID, seq2, br2.readLine(), br2.readLine()));
                    bwError2.newLine();
                    continue;
                }
                lenBarcodeR1=barcodeR1[index1].length();
                lenBarcodeR2=barcodeR2[index2].length();
                reads1_2=new Read[2];
                reads1_2[0]=new Read(read1ID, seq1.substring(lenBarcodeR1), br1.readLine(), br1.readLine().substring(lenBarcodeR1), 33);
                reads1_2[1]=new Read(read2ID, seq2.substring(lenBarcodeR2), br2.readLine(), br2.readLine().substring(lenBarcodeR2), 33);
                bwArray = taxaWriterMap.get(newSet.iterator().next());
                bwArray[0].write(String.join("\n", read1ID, reads1_2[0].getSequence(), reads1_2[0].getDescription(), reads1_2[0].getQualS(33)));
                bwArray[0].newLine();
                bwArray[1].write(String.join("\n", read2ID, reads1_2[1].getSequence(), reads1_2[1].getDescription(), reads1_2[1].getQualS(33)));
                bwArray[1].newLine();
                processedCnt++;
            }
            br1.close();
            br2.close();
            bwError1.flush();
            bwError1.close();
            bwError2.flush();
            bwError2.close();
            for (BufferedWriter[] bufferedWriters: taxaWriterMap.values()){
                bufferedWriters[0].flush();
                bufferedWriters[1].flush();
                bufferedWriters[0].close();
                bufferedWriters[1].close();
            }
            System.out.println("Total read count: "+totalCnt+"\nPassed read count: "+processedCnt);
            System.out.println("The probability of sequencing error in the reads of the "+lib+" library is "+((double)totalCnt-processedCnt)/totalCnt);
        }catch (Exception e){
            e.printStackTrace();
        }
    }



}
