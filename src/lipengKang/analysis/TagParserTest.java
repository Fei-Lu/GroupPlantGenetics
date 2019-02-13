/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;
import lipengKang.analysis.LibraryInfoTest;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
public class TagParserTest {

    LibraryInfoTest li = null;

    public TagParserTest(LibraryInfoTest li, String tagBySampleDirS) {
        this.li = li;
        this.parseFastq(tagBySampleDirS);
    }

    public void parseFastq(String tagBySampleDirS) {
            String[] libs = li.getLibArray();
        for (int i = 0; i < libs.length; i++) {
            String fastqR1 = li.getFastqFileSR1(i);
            String fastqR2 = li.getFastqFileSR2(i);
            HashMap<String, Set<String>> barcodeR1TaxaMap = li.getbarcodeR1TaxaMap(i);
            HashMap<String, Set<String>> barcodeR2TaxaMap = li.getbarcodeR2TaxaMap(i);
            String[] taxaNames = li.getTaxaNames(i);
            String cutter1 = li.getCutter1();
            String cutter2 = li.getCutter2();
            this.splitFastq(fastqR1, fastqR2, barcodeR1TaxaMap, barcodeR2TaxaMap, taxaNames, tagBySampleDirS, cutter1, cutter2);
        }

    }

    private void splitFastq(String fastqR1, String fastqR2,
            HashMap<String, Set<String>> barcodeR1TaxaMap, HashMap<String, Set<String>> barcodeR2TaxaMap, String[] taxaNames, String tagBySampleDirS, String cutter1, String cutter2) {
        HashMap<String, BufferedWriter> taxaWriterMap = new HashMap<>();
        BufferedWriter[] bws = new BufferedWriter[taxaNames.length];
        for (int i = 0; i < taxaNames.length; i++) {
            String outfile = new File(tagBySampleDirS, taxaNames[i] + ".bin").getAbsolutePath();
            BufferedWriter bw = IOUtils.getTextWriter(outfile);
            taxaWriterMap.put(taxaNames[i], bw);
            bws[i] = bw;
        }
        try {
            BufferedReader br1 = null;
            BufferedReader br2 = null;
            if (fastqR1.endsWith(".gz")) {
                br1 = IOUtils.getTextGzipReader(fastqR1);
            } else {
                br1 = IOUtils.getTextReader(fastqR1);
            }
            if (fastqR2.endsWith(".gz")) {
                br2 = IOUtils.getTextGzipReader(fastqR2);
            } else {
                br2 = IOUtils.getTextReader(fastqR2);
            }
            Set<String> bSetR1 = barcodeR1TaxaMap.keySet();
            String[] barcodeR1 = bSetR1.toArray(new String[bSetR1.size()]);
            Set<String> bSetR2 = barcodeR2TaxaMap.keySet();
            String[] barcodeR2 = bSetR2.toArray(new String[bSetR2.size()]);
            Arrays.sort(barcodeR1);
            Arrays.sort(barcodeR2);
            String temp1 = null;
            String temp2 = null;
            int index1 = -1;
            int index2 = -1;
            Set<String> taxaSR1 = null;
            Set<String> taxaSR2 = null;
            BufferedWriter bw = null;
            int totalCnt = 0;
            int processedCnt = 0;
            while ((temp1 = br1.readLine()) != null) {
                totalCnt++;
                temp2 = br2.readLine();
                temp1 = br1.readLine();
                temp2 = br2.readLine();
                index1 = Arrays.binarySearch(barcodeR1, temp1);
                index2 = Arrays.binarySearch(barcodeR2, temp2);
                if (index1 == -1 || index2 == -1) {
                    br1.readLine();
                    br2.readLine();
                    br1.readLine();
                    br2.readLine();
                    continue;
                }
                index1 = -index1 - 2;
                index2 = -index2 - 2;
                taxaSR1 = barcodeR1TaxaMap.get(barcodeR1[index1]);
                taxaSR2 = barcodeR2TaxaMap.get(barcodeR2[index2]);
                Set <String> newSet=new HashSet<>(taxaSR1); 
                
                newSet.retainAll(taxaSR2);
                if (newSet.size() != 1) {
                    br1.readLine();
                    br2.readLine();
                    br1.readLine();
                    br2.readLine();
                    continue;
                }
                temp1 = temp1.substring(barcodeR1[index1].length() , temp1.length());
                temp2 = temp2.substring(barcodeR2[index2].length() , temp2.length());
               
                temp1 = this.getProcessedRead(cutter1, cutter2, temp1);
                temp2 = this.getProcessedRead(cutter1, cutter2, temp2);
                temp1 = br1.readLine();
                temp2 = br2.readLine();
                temp1 = br1.readLine();
                temp2 = br2.readLine();
                processedCnt++;
            }
            for (int i = 0; i < bws.length; i++) {
                bws[i].flush();
                bws[i].close();
            }
            br1.close();
            br2.close();
            System.out.println(totalCnt + "\t" + processedCnt);
        } catch (Exception e) {
            e.printStackTrace();
            System.out.println("Something wrong while reading ");
            System.out.println(fastqR1);
            System.out.println(fastqR2);
            System.exit(1);
            System.out.println("Program quits");
        }

    }

    private String getProcessedRead(String cutter1, String cutter2, String read) {
        int index1 = read.indexOf(cutter1);
        int index2 = read.indexOf(cutter2);
        if (index1 < 0) {
            if (index2 < 0) {
                return read;
            } else {
                return read.substring(0, index2);
            }
        } else {
            if (index2 < 0) {
                return read.substring(0, index1);
            } else {
                if (index1 < index2) {
                    return read.substring(0, index1);
                } else {
                    return read.substring(0, index2);
                }
            }
        }
    }
}
