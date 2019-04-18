package daxing.common;

import analysis.pipeline.grt.LibraryInfo;
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
                this.splitBarcode(fastqR1, fastqR2, barcodeR1TaxaMap, barcodeR2TaxaMap, taxaNames, outputDir, cutter1, cutter2);
            });
        }
    }

    private void splitBarcode(String fastqR1, String fastqR2,
                              HashMap<String, Set<String>> barcodeR1TaxaMap, HashMap<String, Set<String>> barcodeR2TaxaMap, String[] taxaNames, String outputDir, String cutter1, String cutter2){
        HashMap<String, BufferedWriter[]> taxaWriterMap = new HashMap<>();
        BufferedWriter bw1, bw2;
        BufferedWriter[] bwArray;
        for (String str:taxaNames){
            bw1= IOUtils.getNIOTextWriter(new File(outputDir, str+"-1.fq").getAbsolutePath());
            bw2= IOUtils.getNIOTextWriter(new File(outputDir, str+"-2.fq").getAbsolutePath());
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
            String temp1;
            String temp2;
            int index1 = -1;
            int index2 = -1;
            Set<String> taxaSR1;
            Set<String> taxaSR2;
            int totalCnt = 0;
            int processedCnt = 0;
            System.out.println("Parsing " + fastqR1 + "\t" + fastqR2);
            String[] read1;
            String[] read2;
            while ((temp1 = br1.readLine()) != null) {
                temp2 = br2.readLine();
                totalCnt++;
                if (totalCnt%10000000 == 0) {
                    System.out.println("Total read count: "+totalCnt+"\tPassed read count: "+processedCnt);
                }
                read1=new String[4];
                read2=new String[4];
                read1[0]=temp1;
                read1[1]=br1.readLine();
                read1[2]=br1.readLine();
                read1[3]=br1.readLine();
                read2[0]=temp2;
                read2[1]=br2.readLine();
                read2[2]=br2.readLine();
                read2[3]=br2.readLine();
                index1 = Arrays.binarySearch(barcodeR1, read1[1]);
                index2 = Arrays.binarySearch(barcodeR2, read2[1]);
                if (index1 == -1 || index2 == -1) {
                    System.out.println(read1[1]);
                    System.out.println(read1[2]);
                    continue;
                }
                index1 = -index1 - 2;
                index2 = -index2 - 2;
                taxaSR1 = barcodeR1TaxaMap.get(barcodeR1[index1]);
                taxaSR2 = barcodeR2TaxaMap.get(barcodeR2[index2]);
                Set<String> newSet = new HashSet<>(taxaSR1);
                newSet.retainAll(taxaSR2);
                if (newSet.size() != 1) {
                    continue;
                }
                bwArray = taxaWriterMap.get(newSet.iterator().next());
                for (int i = 0; i < read1.length; i++) {
                    if(i==1){
                        bwArray[0].write(read1[i].substring(barcodeR1[index1].length()));
                        bwArray[0].newLine();
                        bwArray[1].write(read2[i].substring(barcodeR2[index2].length()));
                        bwArray[1].newLine();
                    }
                    bwArray[0].write(read1[i]);
                    bwArray[0].newLine();
                    bwArray[1].write(read2[i]);
                    bwArray[1].newLine();
                }
                processedCnt++;
            }
        }catch (Exception e){
            e.printStackTrace();
        }
    }

}
