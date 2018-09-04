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
public class CountRNAseq {
    
    public CountRNAseq(String infileS,String outfileS){
        this.getCountRNAseq(infileS, outfileS);
    }
    public void getCountRNAseq (String infileS,String outfileS){       
        //StringBuilder L3 = new StringBuilder();
        try{           
            BufferedReader br = XueboIOUtils.getTextReader(infileS);    
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            String temp = null;
            int L1 = 1; 
            //String L2 = null;  
            int i =0;
            int CountA = 0;
            int CountB = 0;
            int CountD = 0;
            int j = 0;
            int p = 0;
            int q = 0;
            while (( temp = br.readLine()) != null) {          
                    ++i;                   
                    if (i % 10000 == 0) {
                    System.out.println("Filtering " + i + "....");
                    }
                    String[] tem = temp.split("\t");  
                    String chr = tem[0].substring(8, 9);
                    if(chr.equals("A")){
                        j = j + 1;
                        CountA = CountA + Integer.valueOf(tem[1]);
                    }
                    if(chr.equals("B")){
                        p = p + 1;
                        CountB = CountB + Integer.valueOf(tem[1]);
                    }
                    if(chr.equals("D")){
                        q = q + 1 ;
                        CountD = CountD + Integer.valueOf(tem[1]);
                    }
            }
            bw.write("countA"+ "\t" + j + "\t" + CountA + "\n" + "countB"+ "\t" + p + "\t" +
                     CountB + "\n" + "countD"+ "\t" +  q + "\t" + CountD);
            bw.flush();
            bw.close();
        }
        catch (Exception e){        
            e.printStackTrace();
        }
    }
    
}
