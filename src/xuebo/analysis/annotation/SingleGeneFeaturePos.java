/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class SingleGeneFeaturePos {
    
     public SingleGeneFeaturePos (String infileS,String outfileS1,
            String outfileS2,String outfileS3,String outfileS4,String outfileS5,String outfileS6) {

      this.getSingleGeneFeaturePos(infileS,outfileS1,outfileS2,outfileS3,outfileS4,outfileS5,outfileS6);
    }
    
    public void getSingleGeneFeaturePos (String infileS,String outfileS1,
            String outfileS2,String outfileS3,String outfileS4,String outfileS5,String outfileS6) {
        try {
            
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            BufferedWriter bw3 = IOUtils.getTextWriter(outfileS3);
            BufferedWriter bw4 = IOUtils.getTextWriter(outfileS4);
            BufferedWriter bw5 = IOUtils.getTextWriter(outfileS5);
            BufferedWriter bw6 = IOUtils.getTextWriter(outfileS6);
            
            String temp = null;
            int i = 0 ;
            
            String L2 = null;
            String LUpstream = null;
            String LDownsteam = null;
            String L5UTR = null;
            String LCDS = null;
            String LIntron = null;
            String L3UTR = null;
            String L3 = null;
            String L4 = null;
            String L5 = null;
            String L6 = null;
                    
            while ((temp = br.readLine()) != null) {
                
                ++i;               
                if (i % 5000 == 0) {
                    System.out.println("getSingleGeneFeaturePos " + i + "....");
                }
                
                String[] tem = temp.split("\t");
                boolean outWrite = false;
                
                if (i % 6 == 2) {
//                    L2 = tem[1] + "_" + tem[2]; 
                    L2 = tem[2];
                    int aa = Integer.valueOf(tem[3]) - 2000;
                    int bb = Integer.valueOf(tem[4]) + 2000;
                    if(tem[5].equals("1")){
                        LUpstream = L2 + "\t" + Integer.toString(aa) + ":" + tem[3] + "\n";
                        LDownsteam = L2 + "\t" + tem[4] + ":" + Integer.toString(bb) + "\n";
                    }
                    else{
                        LUpstream = L2 + "\t" + tem[4] + ":" + Integer.toString(bb) + "\n";
                        LDownsteam = L2 + "\t" + Integer.toString(aa) + ":" + tem[3] + "\n";
                    }
                    bw1.write(LUpstream);
                    bw2.write(LDownsteam);
                }
                
                if (i % 6 == 3) {
                    
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 5 ){
                                L5UTR = L2 + "\t"+ tem[1] +  "\n";
                            }                        
//                            bw3.write(L5UTR);
                        }
                    }
                    bw3.write(L5UTR);
                }
                
                if (i % 6 == 4) {
                    
                   if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 25 ){
                                LCDS = L2 + "\t"+ tem[1] +  "\n";
                            }                        
//                            bw3.write(L5UTR);
                        }
                    }
                    bw4.write(LCDS);
                }
                
                if (i % 6 == 5) {
                    
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 25){
                                LIntron = L2 + "\t"+ tem[1] +  "\n";
                            }                        
//                            bw3.write(L5UTR);
                        }
                    }
                    bw5.write(LIntron);
                }
                
                if (i % 6 == 0) {
                    
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 5 ){
                                L3UTR = L2 + "\t"+ tem[1] +  "\n";
                            }                        
//                            bw3.write(L5UTR);
                        }
                    }
                    bw6.write(L3UTR);
                }
                
                bw1.flush();
                bw2.flush();
                bw3.flush();
                bw4.flush();
                bw5.flush();
                bw6.flush();               
            }
            
           bw1.close(); 
           bw2.close(); 
           bw3.close(); 
           bw4.close(); 
           bw5.close(); 
           bw6.close(); 
           
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
}
