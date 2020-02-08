/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class RemoveStreamOverlap {
    
    public RemoveStreamOverlap(String infileS ,String outfileS){
        this.getRemoveStreamOverlap(infileS, outfileS);
    }
    
    
    public void getRemoveStreamOverlap(String infileS ,String outfileS) {
        
        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = IOUtils.getTextGzipReader(infileS);

            } else {

                br = IOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            
            int i = 0;
            int aa = 0;
                      
            String Sitemethy10 = null;
            
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
             
            while (( temp = br.readLine()) != null) {
                
                ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                
                String[] temm = temp.split("\t");
                String[] tem = temm[1].split(":");
                
                if(aa <  Integer.valueOf(tem[0])){
                    aa = Integer.valueOf(tem[1]);
                    bw.write(temp + "\n");
                }
                
                
                
            }
                             
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
            
    }
    
}
