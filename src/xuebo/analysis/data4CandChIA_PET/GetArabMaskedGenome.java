/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4CandChIA_PET;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author xuebozhao
 */
public class GetArabMaskedGenome {
    
    GetArabMaskedGenome(String infileS,String outfileS) {
        
        //this.readGenomelong(infileS,outfileS);
        this.readGenomereplace(infileS, outfileS);
        
    }
    
    private void readGenomelong (String infileS,String outfileS){
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
                    
                    if (i % 1000 == 0) {

                        System.out.println("Filtering " + i + "....");

                    }
                    if (temp.startsWith(">")){
                        bw.write(temp + "\n"); 
                        //continue;
                    }
                    else{
                        bw.write(temp); 
                    }                
            }
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    
    private void readGenomereplace (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
            String L1 = null;
            String L2 = null;
            int i = 0;

            while (( L1 = br.readLine()) != null) {
               
                    if (L1.startsWith(">")){
                        bw.write(L1 + "\n");                      
                    }
                    else{
                        //L2 = L1.substring(15036046, 15137045);
                        //L2 = L1.substring(3557930, 3658929);
                        //L2 = L1.substring(11675025, 11776024);
                        //L2 = L1.substring(13749418, 13900417);
                        //L2 = L1.substring(14158953, 14259952);
                        
//                        L2 = L1.substring(13036046, 17137045);
//                        L2 = L1.substring(1557930, 5658929);
//                        L2 = L1.substring(11749418, 15900417);
//                        L2 = L1.substring(1906022, 6007021);
                        L2 = L1.substring(9675025, 13776024);
//                        
                        
                            StringBuilder sb = new StringBuilder();
                            for (i = 9675025; i < 13776024; ++i)
                            {
                             sb.append("N");
                            }
                            String L3 = sb.toString();

                        bw.write(L1.replace(L2, L3)); 
                    }                
            }
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
    
    
    
    
    
}
