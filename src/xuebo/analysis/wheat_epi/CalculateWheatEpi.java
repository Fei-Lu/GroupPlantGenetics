/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author xuebozhao
 */
public class CalculateWheatEpi {
    
    public CalculateWheatEpi(String infileS,String outfileS){
        //this.readfile(infileS,outfileS);
        this.getCalculatechr(infileS, outfileS);
        //this.readfilecoverage(infileS, outfileS);
    }
    
    public void readfile (String infileS,String outfileS){       
        //StringBuilder L3 = new StringBuilder();
        try{           
            BufferedReader br = XueboIOUtils.getTextReader(infileS);    
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            String temp = null;
            String L2 = null;  
            
            int i = 0;
            while (( temp = br.readLine()) != null) {          
                    ++i;                   
                    if (i % 10000 == 0) {
                    System.out.println("Filtering " + i + "....");
                    }
                    String[] tem = temp.split("\t"); 
                    String L1 = tem[1];
                    String []temm = L1.split("_");
                    int startpos = Integer.valueOf(temm[1]) + Integer.valueOf(tem[2]);
                    L2 = temm[0] + "\t" + temm[3] + "\t" + startpos + "\n";
                    //L3.append(L2);
                    bw.write(L2);
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }

    public void getCalculatechr (String infileS,String outfileS){       
        //StringBuilder L3 = new StringBuilder();
        try{           
            BufferedReader br = XueboIOUtils.getTextReader(infileS);    
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            String temp = null;
            int L1 = 1; 
            //String L2 = null;  
            int i =0;
            int j = 0;
            while (( temp = br.readLine()) != null) {          
                    ++i;                   
                    if (i % 10000 == 0) {
                    System.out.println("Filtering " + i + "....");
                    }
                    String[] tem = temp.split("\t");  
                    
                    if(Integer.valueOf(tem[0]) == L1){
                        j = j + 1;
                        L1 = Integer.valueOf(tem[0]);
                    }
                    else{
                       //j = j + 1 ;
                       bw.write(L1 + "\t" + j + "\n");
                       j = 1;
                       L1 = Integer.valueOf(tem[0]);
                    }
            }
            bw.write(L1 + "\t" + j);
            bw.flush();
            bw.close();
        }
        catch (Exception e){        
            e.printStackTrace();
        }
    }
    
        
    
     public void readfilecoverage (String infileS,String outfileS){       
        //StringBuilder L3 = new StringBuilder();
        try{           
            BufferedReader br = XueboIOUtils.getTextReader(infileS);    
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            String temp = null;
            String L2 = null;  
            
            int i = 0;
            while (( temp = br.readLine()) != null) {          
                    ++i;                   
                    if (i % 10000 == 0) {
                    System.out.println("Filtering " + i + "....");
                    }
                    String[] tem = temp.split("\t"); 
                    String L1 = tem[0];
                    String []temm = L1.split("_");
                    int startpos = Integer.valueOf(temm[1]) + Integer.valueOf(tem[1]);
                    L2 = temm[0] + "\t" + startpos + "\t" + tem[2] + "\n";
                    //L3.append(L2);
                    bw.write(L2);
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
    
}
