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
public class MethylationAnnotation {
    
//    public MethylationAnnotation(String infileS, String outfileS1, String outfileS2,String outfileS3) {
//
//        this.getMethylationAnnotation(infileS,outfileS1, outfileS2,outfileS3);
//    }
    
    public MethylationAnnotation(String infileS ,String outfileS){
        this.getMethylationlevelScore10 (infileS, outfileS);
    }
//    
//    public MethylationAnnotation(String infileS ,String outfileS){
//        //this.getMethylationlevelScore10(infileS, outfileS);
//        this.getOpenchrScore(infileS, outfileS);
//    }
    
    
    public void getMethylationAnnotation(String infileS, String outfileS1, String outfileS2,String outfileS3) {

        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = XueboIOUtils.getTextGzipReader(infileS);

            } else {

                br = XueboIOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
            double methylationScore = 0;           
            String LCpG = null;
            String LCHH = null;
            String LCHG = null;
             
             BufferedWriter bw1 = XueboIOUtils.getTextWriter(outfileS1);
             BufferedWriter bw2 = XueboIOUtils.getTextWriter(outfileS2);
             BufferedWriter bw3 = XueboIOUtils.getTextWriter(outfileS3);
             
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                    
                    String[] tem = temp.split("\t"); 
                    
                    boolean outWrite = false;
                    
                    double methylationC = Double.valueOf(tem[5]);
                    double methylationT = Double.valueOf(tem[5])+ Double.valueOf(tem[6]);
                    methylationScore = methylationC / methylationT;
                    
                    //int mscore = Integer.valueOf(tem[4]);
//                    if(tem[2].equals("-")){
//                        
//                        tem[1] = Integer.toString(Integer.valueOf(tem[1]) -1);
//                    }
                    
                    if(tem[3].equals("CG")){
                        
                        LCpG = tem[0] + "\t" + tem[1] + "\t"  + methylationScore ;
                        bw1.write(LCpG + "\n"); 
                    }
                    
                    if(tem[3].equals("CHH")){
                        
                        LCHH = tem[0] + "\t" + tem[1] + "\t"  + methylationScore ;
                        bw2.write(LCHH + "\n");
                    }
                    
                    if(tem[3].equals("CHG")){
                        
                        LCHG = tem[0] + "\t" + tem[1] + "\t"  + methylationScore;
                        bw3.write(LCHG + "\n");
                    }
                    
//                    if (outWrite) {
//
//                    bw1.write(LCpG + "\n");
//                    bw2.write(LCHH + "\n");
//                    bw3.write(LCHG + "\n");
//                    
//                    bw1.flush();
//                    bw2.flush();
//                    bw3.flush();
//                    }
            }
            bw1.close();
            bw2.close();
            bw3.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    
    }
    
 /**
 *
 * @author xuebozhao
 */
    
    public void getMethylationlevelScore(String infileS ,String outfileS){
        
        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = XueboIOUtils.getTextGzipReader(infileS);

            } else {

                br = XueboIOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
                      
            String LmethyScore = null;
            
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
             
            while (( temp = br.readLine()) != null) {
                
                    i = i + 1;
//                    
                    if (i % 1000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                    
                String[] tem = temp.split("\t"); 
                
                    if (Integer.valueOf(tem[1]) == i){
                        
                        //LmethyScore = tem[0] + "\t" + i + "\t" +tem[2];
                        LmethyScore = tem[0] + "\t" + tem[2];
                        bw.write(LmethyScore + "\n"); 
                        
                    }
                    else{
                        
                        for(int j = i; j< Integer.valueOf(tem[1]);j++){ 
                            
                        //LmethyScore = tem[0] + "\t" +  j + "\t" + "0";
                        LmethyScore = tem[0] + "\t" + "0";
                        bw.write(LmethyScore + "\n");
                        }
                        
                        //LmethyScore = tem[0] + "\t" + i + "\t" +tem[2];
                        LmethyScore = tem[0] + "\t" + tem[2];
                        bw.write(LmethyScore + "\n"); 
                        i = Integer.valueOf(tem[1]);
                    }                 
            }  
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
            
    }
    
   
 /**
 *
 * @author xuebozhao
 */
    
    public void getMethylationlevelScore10(String infileS ,String outfileS){
        
        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = XueboIOUtils.getTextGzipReader(infileS);

            } else {

                br = XueboIOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
                      
            String Sitemethy10 = null;
            
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            
             
            while (( temp = br.readLine()) != null) {
                
                ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                
                String[] tem = temp.split("\t");
                
//                if(Integer.valueOf(tem[0]) == 10){
                    
                    Sitemethy10 = tem[1];
                    bw.write(Sitemethy10 + "\n"); 
//                }
                
                
            }
                             
            bw.flush();
            bw.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
            
    }
    
    /**
 *
 * @author xuebozhao
 */
    
    public void getOpenchrScore(String infileS ,String outfileS){
        
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
                    
                    OpenchrScore =tem[0] + "\t" + Integer.toString(p) + "\t" + "1";
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
