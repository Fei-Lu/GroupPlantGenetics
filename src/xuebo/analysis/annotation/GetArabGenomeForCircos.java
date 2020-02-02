/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import pgl.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class GetArabGenomeForCircos {
    
    GetArabGenomeForCircos(String outfileS) {
        
        this.readBed(outfileS);
        
    }
    public void readBed (String outfileS){
        try{       
            int i = 0;
            int aa = 0;
            String L1 = null;
            String L2 = null;
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            while (aa < 31000000) { 
                    int bb = aa + 500000;
                    L1 = "chr1" + "\t" + aa + "\t" + bb + "\t";
                    if(i % 8 == 0){
                        L2 = L1 + "p36.33" + "\t" + "gneg" + "\n";              
                    }
                    if(i % 8 == 1){
                        L2 = L1 + "p36.31" + "\t" + "gpos75" + "\n";              
                    }
                    if(i % 8 == 2){
                        L2 = L1 + "p36.12" + "\t" + "gpos25" + "\n";              
                    }   
                    if(i % 8 == 3){
                        L2 = L1 + "p36.11" + "\t" + "gpos25" + "\n";              
                    }   
                    if(i % 8 == 4){
                        L2 = L1 + "p36.33" + "\t" + "gneg" + "\n";              
                    }
                    if(i % 8 == 5){
                        L2 = L1 + "p36.33" + "\t" + "gpos50" + "\n";              
                    }
                    if(i % 8 == 6){
                        L2 = L1 + "p36.33" + "\t" + "red" + "\n";              
                    }
                    if(i % 8 == 7){
                        L2 = L1 + "p36.33" + "\t" + "gpos100" + "\n";              
                    }
                    aa = aa + 500000;
                    
                    i = i + 1;
                    
                    bw.write(L2); 
                }
            bw.flush(); 
            
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
