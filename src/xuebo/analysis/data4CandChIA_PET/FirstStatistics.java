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
public class FirstStatistics {
    

    FirstStatistics(String infileS,String outfileS) {
        
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
            double a1 = 0;
            int a2 = 0;
            int a3 = 0;
            double f = 0;
            double fo = 0;
            int count = 1;
            int count0 =0;
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 1000 == 0) {

                        System.out.println("Filtering " + i + "....");

                    }


                    String[] tem = temp.split("\t");
                    count0 = Integer.valueOf(tem[3]);
                    

//                if ( count > 1){

                    a1 = 7.1 * Double.parseDouble(tem[0]);
                    a2 = 100 * Integer.parseInt(tem[1]);
                    a3 = 1000 * Integer.parseInt(tem[2]);
                                   
                    f = a1 + a2 + a3;
                    
                    L1 = "chr" + tem[0] + "\t" + tem[1] + "\t" + tem[2] + "\t";
                  

                    if(f == fo){

                        count  = count + 1;
                        L2 = L1;
                        
                    }
                    else{ 

                        fo = f; 
                        
                        if(i > 1) {                       
                            
                            bw.write(L2 + count + "\n"); 
                            bw.flush(); 
                            count = Integer.valueOf(tem[3]);
                        }
                        L2 = L1;
                        count = count0;
                    }                 
    //                bw.close();
//                }
            }
            bw.write(L2 + count + "\n"); 
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
}
