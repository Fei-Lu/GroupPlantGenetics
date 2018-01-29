/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package feilu;

import format.table.RowTable;
import java.io.BufferedWriter;
import java.util.HashMap;
import zhouyao.analysis.wheatHapMap.IOUtils;

/**
 *
 * @author feilu
 */
public class Start {
    
    public Start () {
        this.test();
    }
    
    public void test () {
        String barcodeFileS = "/Users/feilu/Downloads/barcode.txt";
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
        }
        String s = null;
        int index = barcodeIndexMap.get(s.substring(0, 8));
        BufferedWriter[] bws = new BufferedWriter[rowNumber];
        
        String seq = "ATGCCGTAGC";
        
    }
    
    public static void main (String[] args) {
        new Start();
    }
    
}
