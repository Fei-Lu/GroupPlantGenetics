/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import xuebo.analysis.data4C.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class MethylationAnnotation {
    
    public MethylationAnnotation(String infileS, String outfileS1, String outfileS2,String outfileS3) {

        this.getMethylationAnnotation(infileS,outfileS1, outfileS2,outfileS3);
//        this.outputFiles(outfileS1, outfileS2);
    }
    
    public void getMethylationAnnotation(String infileS, String outfileS1, String outfileS2,String outfileS3) {

        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = IOUtils.getTextGzipReader(infileS);

            } else {

                br = IOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
            double methylationScore = 0;           
            String LCpG = null;
            String LCHH = null;
            String LCHG = null;
             
             BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfileS1);
             BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfileS2);
             BufferedWriter bw3 = IOUtils.getTextGzipWriter(outfileS3);
             
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                    
                    String[] tem = temp.split("\t"); 
                    
                    boolean outWrite = false;
                    
//                    int methylationC = Integer.valueOf(tem[6]);
//                    int methylationT = Integer.valueOf(tem[7]);
//                    methylationScore = methylationC / (methylationC + methylationT);
                    
                    //int mscore = Integer.valueOf(tem[4]);
                    
                    if(tem[9].equals("CpG")){
                        
                        LCpG = tem[0] + "\t" + tem[1] + "\t" + tem[2] + "\t" + tem[4] ;
                        bw1.write(LCpG + "\n"); 
                    }
                    
                    if(tem[9].equals("CHH")){
                        
                        LCHH = tem[0] + "\t" + tem[1] + "\t"  + tem[2] + "\t" + tem[4] ;
                        bw2.write(LCHH + "\n");
                    }
                    
                    if(tem[9].equals("CHG")){
                        
                        LCHG = tem[0] + "\t" + tem[1] + "\t" + tem[2] + "\t" + tem[4] ;
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
    
}
