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
public class SingleSingleGeneFeaturePos {
    
     public SingleSingleGeneFeaturePos (String infileS,String outfileS1,
            String outfileS2,String outfileS3,String outfileS4,String outfileS5,String outfileS6,String outfileS7) {

      this.getSingleSingleGeneFeaturePos(infileS,outfileS1,outfileS2,outfileS3,outfileS4,outfileS5,outfileS6,outfileS7);
    }
    
    public void getSingleSingleGeneFeaturePos (String infileS,String outfileS1,
            String outfileS2,String outfileS3,String outfileS4,String outfileS5,String outfileS6,String outfileS7) {
        try {
            
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw1 = IOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextWriter(outfileS2);
            BufferedWriter bw3 = IOUtils.getTextWriter(outfileS3);
            BufferedWriter bw4 = IOUtils.getTextWriter(outfileS4);
            BufferedWriter bw5 = IOUtils.getTextWriter(outfileS5);
            BufferedWriter bw6 = IOUtils.getTextWriter(outfileS6);
            BufferedWriter bw7 = IOUtils.getTextWriter(outfileS7);
            
            String temp = null;
            int i = 0 ;
            
            String L2 = null;
            String Lstrand = null;
            String LUpstream = null;
            String LDownsteam = null;
            String Lall = null;
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
                //boolean outWrite = false;
          
             
                if (i % 6 == 2) {
//                    L2 = tem[1] + "_" + tem[2]; 
                    L2 = tem[2];
                    Lstrand = tem[5];
                    int aa = Integer.valueOf(tem[3]) - 2000;
                    int bb = Integer.valueOf(tem[4]) + 2000;
                    if(tem[5].equals("1")){
                        LUpstream = L2 + "\t" + Integer.toString(aa) + ":" + tem[3] + "\t"  + Lstrand + "\n";
                        LDownsteam = L2 + "\t" + tem[4] + ":" + Integer.toString(bb) + "\t"  + Lstrand + "\n";
                        if((Integer.valueOf(tem[4]) - Integer.valueOf(tem[3])) > 100 ){
                        Lall = L2 + "\t" + Integer.toString(Integer.valueOf(tem[3])) + ":" + Integer.toString(Integer.valueOf(tem[4]))+ "\t"  + Lstrand + "\n";
                        }
                    }
                         
                    else{
                        LUpstream = L2 + "\t" + tem[4] + ":" + Integer.toString(bb) + "\t"  + Lstrand + "\n";
                        LDownsteam = L2 + "\t" + Integer.toString(aa) + ":" + tem[3] + "\t"  + Lstrand + "\n";
                        if((Integer.valueOf(tem[4]) - Integer.valueOf(tem[3])) > 100 ){
                        Lall = L2 + "\t" + Integer.toString(Integer.valueOf(tem[3])) + ":" + Integer.toString(Integer.valueOf(tem[4])) + "\t"  + Lstrand + "\n";
                        }
                    }
                    bw1.write(LUpstream);
                    bw2.write(LDownsteam);
                    bw7.write(Lall);
                }
                
                if (i % 6 == 3) {
                    
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            //System.out.println(temt[k]);
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 100 ){
                                L5UTR = L2 + "\t"+ temt[k] + "\t"  + Lstrand + "\n";
                                bw3.write(L5UTR);
                            }                        
//                            bw3.write(L5UTR);
                        }
                        
                    }
                    
                }
                
                if (i % 6 == 4) {
                    
                   if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 100 ){
                                LCDS = L2 + "\t"+ temt[k] + "\t"  + Lstrand + "\n";
                                //outWrite = true;
                                bw4.write(LCDS);
                            }                        
//                            bw3.write(L5UTR);
                        }
                        
                    }
                   
                }
                
                if (i % 6 == 5) {
                    
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 100){
                                LIntron = L2 + "\t"+ temt[k] + "\t"  + Lstrand + "\n";
                                //outWrite = true;
                                bw5.write(LIntron);
                            }                        
//                            bw3.write(L5UTR);
                        }
                        
                    }
                    
                }
                
                if (i % 6 == 0) {
//                    if(i>17){
//                        System.out.printf("testing");
//                    }
                    if (!tem[1].startsWith("NA")) {
                        
                        String [] temt = tem[1].split(";");
                        for (int k = 0; k < temt.length; k++) {
                            String[] te = temt[k].split(":");
                            if((Integer.valueOf(te[1]) - Integer.valueOf(te[0])) > 100){
                                L3UTR = L2 + "\t"+ temt[k] + "\t"  + Lstrand + "\n";
                                //outWrite = true;
                                bw6.write(L3UTR);
                            }                        
//                            bw3.write(L5UTR);
                        }
                        
                    }
                    
                }
                
                bw1.flush();
                bw2.flush();
                bw3.flush();
                bw4.flush();
                bw5.flush();
                bw6.flush();   
                bw7.flush();
            }
            
           bw1.close(); 
           bw2.close(); 
           bw3.close(); 
           bw4.close(); 
           bw5.close(); 
           bw6.close(); 
           bw7.close();
           
        }
         catch (Exception e) {
            e.printStackTrace();
        }
    }
}
