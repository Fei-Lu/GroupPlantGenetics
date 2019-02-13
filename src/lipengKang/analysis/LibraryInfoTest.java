/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import format.table.RowTable;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 *
 * @author kanglipeng
 */
public class LibraryInfoTest {

    String[] libs = null;
    String[][] taxaNames = null;
    String[][] barcodeR1 = null;
    String[][] barcodeR2 = null;
    HashMap<String, Set<String>>[] barcodeR1TaxaMaps = null;
    HashMap<String, Set<String>>[] barcodeR2TaxaMaps = null;

    HashMap<String, String>[][] taxaOutputFileSMaps = null;
    String[] libFastqsR1 = null;
    String[] libFastqsR2 = null;
    String cutter1 = "GGATCC";
    String cutter2 = "CCGG";

    public LibraryInfoTest(String barcodeFileS, String libraryFastqMapFileS) {
        this.parseBarcode(barcodeFileS, libraryFastqMapFileS);
    }

    public int getLibraryNumber() {
        return libs.length;
    }

    public String[] getLibArray() {
        String[] na = new String[libs.length];
        System.arraycopy(libs, 0, na, 0, libs.length);
        return na;
    }

    public String getCutter1() {
        return this.cutter1;
    }

    public String getCutter2() {
        return this.cutter2;
    }

    public String[] getTaxaNames(int index) {
        return taxaNames[index];
    }

    public String[] getLibBarcodeR1(int index) {
        String[] na = new String[barcodeR1[index].length];
        System.arraycopy(barcodeR1[index], 0, na, 0, barcodeR1[index].length);
        return na;
    }

    public String[] getLibBarcodeR2(int index) {
        String[] na = new String[barcodeR2[index].length];
        System.arraycopy(barcodeR2[index], 0, na, 0, barcodeR2[index].length);
        return na;
    }

    public HashMap<String, Set<String>> getbarcodeR1TaxaMap(int index) {
        return barcodeR1TaxaMaps[index];
    }

    public HashMap<String, Set<String>> getbarcodeR2TaxaMap(int index) {
        return barcodeR2TaxaMaps[index];
    }

    public String getFastqFileSR1(int index) {
        return libFastqsR1[index];
    }

    public String getFastqFileSR2(int index) {
        return libFastqsR2[index];
    }

    private void parseBarcode(String barcodeFileS, String libraryFastqMapFileS) {
        RowTable<String> t = new RowTable<>(barcodeFileS);
        List<String> l = t.getColumn(0);
        Set<String> s = new HashSet(l);
        libs = s.toArray(new String[s.size()]);
        Arrays.sort(libs);
        System.out.println(Arrays.toString(libs));
        taxaNames = new String[libs.length][];
        barcodeR1 = new String[libs.length][];
        barcodeR2 = new String[libs.length][];
        barcodeR1TaxaMaps = new HashMap[libs.length];
        barcodeR2TaxaMaps = new HashMap[libs.length];
        libFastqsR1 = new String[libs.length];
        libFastqsR2 = new String[libs.length];

        for (int i = 0; i < libs.length; i++) {
            List<String> nameList = new ArrayList<>();
            List<String> barcodeR1List = new ArrayList<>();
            List<String> barcodeR2List = new ArrayList<>();
            for (int j = 0; j < t.getRowNumber(); j++) {
                if (!t.getCell(j, 0).equals(libs[i])) {
                    continue;
                }
                StringBuilder sb = new StringBuilder();
                sb.append(t.getCell(j, 6)).append("_").append(t.getCell(j, 0)).append("_").append(t.getCell(j, 1));
                nameList.add(sb.toString());
                barcodeR1List.add(t.getCell(j, 3));
                barcodeR2List.add(t.getCell(j, 4));
            }

            taxaNames[i] = nameList.toArray(new String[nameList.size()]);
            barcodeR1[i] = barcodeR1List.toArray(new String[nameList.size()]);
            barcodeR2[i] = barcodeR2List.toArray(new String[nameList.size()]);
            barcodeR1TaxaMaps[i] = new HashMap<>();
            barcodeR2TaxaMaps[i] = new HashMap<>();

      /*      for (int j = 0; j < barcodeR1List.size(); j++) {
                s = new HashSet<>();
                s.add(taxaNames[i][j]);
                barcodeR1TaxaMaps[i].put(barcodeR1List.get(j), s);
                s = new HashSet<>();
                s.add(taxaNames[i][j]);
                barcodeR2TaxaMaps[i].put(barcodeR2List.get(j), s);

            }*/

            for (int j = 0; j < barcodeR1List.size(); j++) {
                s = new HashSet<>();
                barcodeR1TaxaMaps[i].put(barcodeR1List.get(j), s);
                s = new HashSet<>();
                barcodeR2TaxaMaps[i].put(barcodeR2List.get(j), s);
            }

            for (int j = 0; j < taxaNames[i].length; j++) {
                s = barcodeR1TaxaMaps[i].get(barcodeR1[i][j]);
                s.add(taxaNames[i][j]);
                barcodeR1TaxaMaps[i].put(barcodeR1[i][j], s);
                s = barcodeR2TaxaMaps[i].get(barcodeR2[i][j]);
                s.add(taxaNames[i][j]);
                barcodeR2TaxaMaps[i].put(barcodeR2[i][j], s);
            }
            int m = 3;
        }

        t = new RowTable<>(libraryFastqMapFileS);
        t.sortAsText(0);
        List<String> lList = t.getColumn(0);
        for (int i = 0; i < libs.length; i++) {
            int index = Collections.binarySearch(lList, libs[i]);
            if (index < 0) {
                System.out.println("Something in the libraryFastqMapFileS. Program quits");
                System.exit(1);
            }
            this.libFastqsR1[i] = t.getCell(index, 1);
            this.libFastqsR2[i] = t.getCell(index, 2);
        }
    }

}
