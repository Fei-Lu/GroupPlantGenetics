/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import pgl.infra.utils.IOUtils;
//import utils.IoUtils;

/**
 *
 * @author xuebozhao
 */
public class Containing {

    Containing(String infileS1, String infileS2, String outfileS1, String outfileS2) {

        this.readFile(infileS1, infileS2, outfileS1, outfileS2);
//        this.outputFiles(outfileS1, outfileS2);
    }

    public void readFile(String infileS1, String infileS2, String outfileS1, String outfileS2) {

        try {

            BufferedReader br1;
            BufferedReader br2;

            if (infileS1.endsWith("gz")) {

                br1 = IOUtils.getTextGzipReader(infileS1);
                br2 = IOUtils.getTextGzipReader(infileS2);

            } else {

                br1 = IOUtils.getTextReader(infileS1);
                br2 = IOUtils.getTextReader(infileS2);

            }

            String temp1 = null;
            String temp2 = null;

            BufferedWriter bw1 = IOUtils.getTextGzipWriter(outfileS1);
            BufferedWriter bw2 = IOUtils.getTextGzipWriter(outfileS2);

            int i = 0;

            String L11 = null,
            L12 = null,
            L13 = null,
            L14 = null;
 
            String L21 = null,
                    L22 = null,
                    L23 = null,
                    L24 = null;

            
            int L2t = 0;
            int aa = 0;
            int bb = 0;
            int cc = 0;
            int dd = 0;
            
            while ((temp1 = br1.readLine()) != null
                    && (temp2 = br2.readLine()) != null) {
//                System.out.println(temp1.substring(0 ,1)); 
//                System.out.println(i++); 

                ++i;
                
                if (i % 500000 == 0) {

                    System.out.println("Filtering " + i + "....");

                }

                boolean outWrite = false;

                if (i % 4 == 1) {

                    L11 = temp1;
                    L21 = temp2;

                }
                if (i % 4 == 2) {
                    int len = temp1.length();
                    
                    aa = temp1.indexOf("TTAA");
                    bb = temp2.indexOf("AGATCT");
                
                    cc = temp2.indexOf("TTAA");
                    dd = temp1.indexOf("AGATCT");

                    String a1, a2, a3, a4;

                    a1 = temp1.substring(72, 78);
                    a2 = temp2.substring(72, 78);
                    a3 = temp1.substring(79, 83);
                    a4 = temp2.substring(79, 83);

                    if (a1.equals("AGATCT") && a4.equals("TTAA")) {
                        
                        if(temp1.contains("TTAA") && temp2.contains("AGATCT")){
//                            int aa = temp1.indexOf("TTAACC");
//                            int bb = temp2.indexOf("AGATCT");
                            if(aa > 78 && bb>83){
                            L12 = temp1.substring(78, aa);
                            L22 = temp2.substring(83, bb);
                            }
                            if (L12.length() == L22.length()){
                                L2t = 11;
                                
                            }                          
                        }
                        
//                        if(temp2.contains("AGATCT")){
//                            int bb = temp2.indexOf("AGATCT");
//                            L22 = temp2.substring(83, bb);
//                            L2t = 12;
//                        }
                        
                        else{ 
                            L12 = temp1.substring(78,len);
                            L22 = temp2.substring(83,len);
                            L2t = 1;
                        }
                        

                    } //                       if((temp1.contains("TTAA") && temp2.contains("AGATCT"))  || 
                    //                               (temp1.contains("AGATCT") && temp2.contains("TTAA"))){                              
                    //                       }        
                    //                       else{
                    //                           String d2 = temp1;
                    
                    else if (a2.equals("AGATCT") && a3.equals("TTAA")) {
                        
                        if(temp2.contains("TTAA") && temp1.contains("AGATCT")){
//                            int cc = temp2.indexOf("TTAACC");
//                            int dd = temp1.indexOf("AGATCT");
                            if(cc > 78 && dd > 83){
                            L12 = temp2.substring(78, cc);
                            L22 = temp1.substring(83, dd);
                            }
                            if (L12.length() == L22.length()){
                                L2t = 13;
                            }                         
                        }
                        
//                        if(temp1.contains("AGATCT")){
//                            int bb = temp1.indexOf("AGATCT");
//                            L22 = temp1.substring(83, bb);
//                            L2t = 14;
//                        }

                        else {
                         
                            L12 = temp1.substring(83, len);                       
                            L22 = temp2.substring(78, len);
                            L2t = 2;
                        }

                    }
                }

                if (i % 4 == 3) {

                    L13 = temp1;
                    L23 = temp2;

                }
                if (i % 4 == 0) {

                    int sum1 = 0;
                    int sum2 = 0;
                    int l1 = temp1.length();
                    int l2 = temp2.length();

//                       for (int j = 0; j < l1; j++) {
//                           sum1 += Integer.getInteger(L12.substring(j));
//                        }
//                       for (int k = 0; k < l2; k++) {
//                           sum2 += Integer.getInteger(L22.substring(k));
//                        }
//                       long mean1 = 0;
//                       long mean2 = 0;
//                       mean1 = sum1 / l1 ;
//                       mean2 = sum1 / l2 ;
//                       if(mean1 >=53 && mean2 >=53){
//                           L14 = 
//                       }
                    String L140 = null;
                    String L240 = null;

                    if (L2t == 1) {

                        L140 = temp1.substring(78, l1);
                        L240 = temp2.substring(83, l2);

                        for (int j = 0; j < L140.length() - 1; j++) {

                            sum1 += stringToAscii(L140.substring(j, j + 1));

                        }

                        for (int k = 0; k < L240.length() - 1; k++) {

                            sum2 += stringToAscii(L240.substring(k, k + 1));

                        }

                        long mean1 = 0;
                        long mean2 = 0;

                        mean1 = sum1 / (l1 - 78);
                        mean2 = sum2 / (l2 - 83);

//                            for (int j = 0; j < L140.length(); j++){
//                                sum1 += (Integer.getInteger(L140.substring(j,j+1)) 
//                                        + Integer.getInteger(L240.substring(j,j+1)));
//                            }
//                            int mean1 = (int) (sum1/(2*L140.length()));
                        if (mean1 > 52 && mean2 > 52) {

//                              if(mean1 > 52){
                            L14 = L140;
                            L24 = L240;

                            outWrite = true;        

                        }
                    } 
                    
                    if (L2t == 2) {

                        L140 = temp1.substring(83, l1);
                        L240 = temp2.substring(78, l2);

                        for (int j = 0; j < L140.length() - 1; j++) {

                            sum1 += stringToAscii(L140.substring(j, j + 1));

                        }
                        for (int k = 0; k < L240.length() - 1; k++) {

                            sum2 += stringToAscii(L240.substring(k, 1 + k));

                        }

                        long mean1 = 0;
                        long mean2 = 0;

                        mean1 = sum1 / (l1 - 83);
                        mean2 = sum1 / (l2 - 78);

                        if (mean1 > 52 && mean2 > 52) {

                            L14 = L140;
                            L24 = L240;

                            outWrite = true;

                        }
                    }
                    
                    if (L2t == 11){
                        
                        if(aa > 78 && bb > 83){
                            
                        L14 = temp1.substring(78, aa);
                        L24 = temp2.substring(83, bb);
                        
                        outWrite = true;
                        }
                    }
                    
                    if (L2t == 13){
                        if(cc > 78 && dd > 83){
                            
                        L14 = temp2.substring(78, cc);
                        L24 = temp1.substring(83, dd);
                        
                        outWrite = true;
                        }
                    }
                    L2t = 0;
                    
                }

                if (outWrite) {

                    bw1.write(L11 + "\n");
                    bw1.write(L12 + "\n");
                    bw1.write(L13 + "\n");
                    bw1.write(L14 + "\n");
                    bw2.write(L21 + "\n");
                    bw2.write(L22 + "\n");
                    bw2.write(L23 + "\n");
                    bw2.write(L24 + "\n");
                    bw1.flush();
                    bw2.flush();

                }
            }

            bw1.close();
            bw2.close();

        } catch (Exception e) {
            e.printStackTrace();
        }
    }

//    public void outputFiles(String outfileS1,String outfileS2){
//      BufferedWriter bw1 = IoUtils.getTextWriter(outfileS1);
//      BufferedWriter bw2 = IoUtils.getTextWriter(outfileS2); 
//      for (int m = 0; m < (); m++){
//          
//      }
//    }
    public int stringToAscii(String value) {
        int sbu = 0;
        char[] chars = value.toCharArray();
        for (int i = 0; i < chars.length; i++) {
            sbu = (int) chars[i];
        }
        return sbu;
    }
}
