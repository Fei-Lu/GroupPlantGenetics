/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4C;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author xuebozhao
 */
public class BedUnSorted {
    
    BedUnSorted(String infileS,String outfileS) {
        
        this.readBed(infileS,outfileS);
        
    }
    
    public void readBed (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);
            int i = 1;
            String temp = null;
            String L1 = null;
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10000 == 0) {

                    System.out.println("Filtering " + i + "....");

                    }


                    String[] tem = temp.split("\t");
                    String[] temfirst = tem[0].split("_");
                    
                    int aa = Integer.valueOf(temfirst[1]) + Integer.valueOf(tem[1]);
                    int bb = Integer.valueOf(temfirst[1]) + Integer.valueOf(tem[2]);
                    
                    L1 = temfirst[0] + "\t" + aa + "\t" + bb + "\n";
                    
                    bw.write(L1);
                    bw.flush();
            }
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
