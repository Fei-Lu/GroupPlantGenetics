/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import xuebo.analysis.data4CandChIA_PET.XueboIOUtils;

/**
 *
 * @author xuebozhao
 */
public class ParallelMean {
    
    public ParallelMean(String infileS ,String outfileS){
        this.getParallelMean(infileS, outfileS);
    }
    
    public void getParallelMean(String infileS ,String outfileS){
        
        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = XueboIOUtils.getTextGzipReader(infileS);

            } else {

                br = XueboIOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
                      
            String OpenchrScore = null;
            
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
             
            while (( temp = br.readLine()) != null) {
                
                    i = i + 1;
//                    
//                    if (i % 100 == 0) {
//
//                    System.out.println("MethylationAnnotation" + i + "....");
//
//                    }
                    
                String[] tem = temp.split("\t"); 
                
                int aa = Integer.valueOf(tem[1]);
                int bb = Integer.valueOf(tem[2]);
                
                for(int p = aa ;p < bb;p++){
                    
                    OpenchrScore = "1";
                    bw.write(OpenchrScore + "\n"); 
                }
//                
//                    if (Integer.valueOf(tem[1]) == i){
//                        
////                        LmethyScore = tem[0] + "\t" + i + "\t" +tem[2];
//                        LmethyScore = tem[0] + "\t" + tem[2];
//                        bw.write(LmethyScore + "\n"); 
//                        
//                    }
//                    else{
//                        
//                        for(int j = i; j< Integer.valueOf(tem[1]);j++){ 
//                            
////                        LmethyScore = tem[0] + "\t" +  i + "\t" + "0";
//                        LmethyScore = tem[0] + "\t" + "0";
//                        bw.write(LmethyScore + "\n");
//                        }
//                        
////                        LmethyScore = tem[0] + "\t" + tem[1] + "\t" +tem[2];
//                        LmethyScore = tem[0] + "\t" + tem[2];
//                        bw.write(LmethyScore + "\n"); 
//                        i = Integer.valueOf(tem[1]);
//                    }                 
            }  
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
            
    }
                    
                    
                    
                    
                    
}
