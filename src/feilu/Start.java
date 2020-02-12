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
        this.testCallableFuture();
    }
    
    public void testCallableFuture () {
        new CallableFuture();
    }
     
    public static void main (String[] args) {
        new Start();
    }
    
}


