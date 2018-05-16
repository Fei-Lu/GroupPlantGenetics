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
public class GetCircosBedMaize {
    GetCircosBedMaize(String infileS,String outfileS) {
        
        this.readBed(infileS,outfileS);
        
    }
    
    public void readBed (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
            String temp = null;
            String L1 = null;
            String L2 = null;
            int i = 0;

            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10 == 0) {

                        System.out.println("Filtering " + i + "....");

                    }

                    String[] tem = temp.split("\t");
                    String[] tem1 = tem[0].split("_");
                    String[] tem2 = tem[1].split("_");
                    
                    String aa = tem1[0].substring(1, 2);
                    String bb = tem2[0].substring(1, 2); 
                    
                    int cc  = Integer.valueOf(tem1[1]) + 10;
                    int dd  = Integer.valueOf(tem2[1]) + 10;
                    
//                    L1 = "chr" + tem[0] + "\t" + tem[1] + "\t" + aa + "\t" + "chr2" + "\t" + "7808877" + "\t" + "7810811";
                    
                    L1 = "chr" + aa + "\t" + tem1[1] + "\t" + cc + "\t" + "chr" + bb + "\t" + tem2[1] + "\t" + dd ;
                            
                    bw.write(L1 + "\n"); 
            }
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
}
