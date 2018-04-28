/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import analysis.maize2k.*;
import format.table.RowTable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import utils.FArrayUtils;
import utils.IOFileFormat;

/**
 *
 * @author feilu
 */
public class HapMapTaxaProcessorcp {
    
    public HapMapTaxaProcessorcp () {
        this.mkSampleTaxaMap();
        String infileS1 = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt";
        String infileS2 = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/t.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaBam.txt";
        updateTaxaBamMap(infileS1, infileS2, outfileS);
        //some related files please check the bath /Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor, it will  give you some instructions
        this.testTaxaDuplicates();
    }
    
    private void testTaxaDuplicates () {
        String infileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaBam.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> taxaList = t.getColumn(0);
        String[] taxaArray = taxaList.toArray(new String[taxaList.size()]);
        String[] uTaxaArray = FArrayUtils.getUniqueStringArray(taxaArray);
        Arrays.sort(uTaxaArray);
        int[] cnts = new int[uTaxaArray.length];
        for (int i = 0; i < taxaArray.length; i++) {
            int index = Arrays.binarySearch(uTaxaArray, taxaArray[i]);
            cnts[index]++;
        }
        for (int i = 0; i < cnts.length; i++) {
            if (cnts[i] < 2) continue;
            System.out.println(uTaxaArray[i]);
        }
    }
    
    private void mkSampleTaxaMap () {
        String infileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/initialTaxaNameMap.txt";
        String outfileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/HapMapTaxaProcessor/taxaNameMap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = new ArrayList(); 
        for (int i = 0; i < t.getRowNumber(); i++) {
            String currentName = t.getCell(i, 1);
            if (!isStandardTaxonName(currentName)) {
                currentName = this.getStandardTaxonName(currentName);
            }
            l.add(currentName);
        }
        t.insertColumn("Taxa", 2, l);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    private String getStandardTaxonName (String taxonName) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < taxonName.length(); i++) {
            int c = (int)taxonName.charAt(i);
            if ((47 < c && c < 58) || ((64 < c && c < 91)) || (96 < c && c < 123) || c == 45 || c ==95) {
                sb.append(taxonName.charAt(i));
            }
            else {
                if (sb.length() != 0) {
                    sb.append("_");
                }
            }
        }
        return sb.toString();
    }
    
    private boolean isStandardTaxonName (String taxonName) {
        int cnt = 0;
        for (int i = 0; i < taxonName.length(); i++) {
            int c = (int)taxonName.charAt(i);
            if ((47 < c && c < 58) || ((64 < c && c < 91)) || (96 < c && c < 123) || c == 45 || c ==95) {
                cnt++;
            }
        }
        if (cnt == taxonName.length()) return true;
        return false;
    }
    
    public static void updateTaxaBamMap (String infileS1, String infileS2, String outfileS) {
        RowTable<String> t = new RowTable<>(infileS1);
        HashMap<String, String> sampleTaxaMap = new HashMap<>();
        for (int i = 0; i < t.getRowNumber(); i++) {
            sampleTaxaMap.put(t.getCell(i, 0), t.getCell(i, 2));
        }
        t = new RowTable<>(infileS2);
        for (int i = 0; i < t.getRowNumber(); i++) {
            String c = t.getCell(i, 0);
            t.setCell(i, 0, sampleTaxaMap.get(c));
        }
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
}
