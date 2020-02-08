/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package feilu;

import com.koloboke.collect.map.hash.HashByteByteMap;
import pgl.infra.dna.BaseEncoder;
import pgl.infra.table.RowTable;
import java.io.BufferedWriter;
import java.util.HashMap;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author feilu
 */
public class Start {
    
    public Start () {
        this.test();
    }
    
    public void test () {
     /*   String barcodeFileS = "/Users/feilu/Downloads/barcode.txt";
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
        */
        String ori = "ATGCAGTC";
        byte[] oriByte = ori.getBytes();
        HashByteByteMap ascMap = BaseEncoder.getAscIIByteMap();
        byte[] tranByte = new byte[oriByte.length];
        for (int i = 0; i < oriByte.length; i++) {
            tranByte[i] = ascMap.get(oriByte[i]);
        }
        int c = BaseEncoder.getShortSeqFromByteArray(tranByte);
        System.out.println(BaseEncoder.getSequenceFromShort(c));    
        
         }
    
    public static void main (String[] args) {
        new Start();
    }
    
}
