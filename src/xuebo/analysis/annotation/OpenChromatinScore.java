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
public class OpenChromatinScore {
    
    public OpenChromatinScore (String infileS,String outfileS ) {

      this.getOpenChromatinScore(infileS,outfileS);
    }
    
    public void getOpenChromatinScore (String infileS,String outfileS) {
        try {
            
            BufferedReader br;
            br = IOUtils.getTextReader(infileS);
            String temp = null;
//          StringBuilder sb = new StringBuilder(temp);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            
            int i = 0;
            
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                    
                    String[] tem = temp.split("\t"); 
            
            }
            
            
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
