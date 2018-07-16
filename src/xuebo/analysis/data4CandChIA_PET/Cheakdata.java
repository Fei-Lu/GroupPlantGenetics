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
public class Cheakdata {
    
    Cheakdata(String infileS,String outfileS) {
        
        //this.readFasta(infileS,outfileS);
        this.spliteChrFasta(infileS,outfileS);
        
    }
    
    private void readFasta (String infileS,String outfileS){
        try{
            
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);           
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            String temp = null;
            String temp2 = null;
            int i = 0;
            boolean check = false;
            StringBuilder sb = new StringBuilder(); 
            while (( temp = br.readLine()) != null) {
                ++i;                    
                if (i % 1000000 == 0) {
                    System.out.println("Filtering " + i + "....");
                }
                if(temp.startsWith(">5")){
                    check = true;
                    continue;
                }
                if(check){
                    sb.append(temp);
                }
                if(check & temp.startsWith(">")) break;
            }
            String fa = sb.toString().substring(11755725,11762652);
            sb = new StringBuilder();
            bw.write(">chr5Sub\n");
            bw.write(fa);
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
        
        
    }
    
    
    
    private void spliteChrFasta (String infileS,String outfileS){
        try{
            BufferedReader br;            
            br = XueboIOUtils.getTextReader(infileS);           
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
//            BufferedWriter bw2 = XueboIOUtils.getTextWriter(outfileS);
//            BufferedWriter bw3 = XueboIOUtils.getTextWriter(outfileS);
//            BufferedWriter bw4 = XueboIOUtils.getTextWriter(outfileS);
//            BufferedWriter bw5 = XueboIOUtils.getTextWriter(outfileS);
            StringBuilder sb = new StringBuilder(); 
            boolean spliteChr = false;
            String temp2 = null;
            
            while (( temp2 = br.readLine()) != null) {
                if(temp2.startsWith(">4")){
                     spliteChr = true;
                     continue;
                }
                if(spliteChr){
                    sb.append(temp2 + "\n");
                }
                if(spliteChr & temp2.startsWith(">")) break;
            }
            String fa = sb.toString();
            sb = new StringBuilder();
            
            bw.write(">Splitechr4\n");
            bw.write(fa);
            
            bw.flush(); 
            bw.close();
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
}
