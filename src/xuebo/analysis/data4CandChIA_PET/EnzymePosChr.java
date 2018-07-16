/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4CandChIA_PET;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import static xuebo.analysis.data4CandChIA_PET.ReducedLibrary.getWriteStreamAppend;

/**
 *
 * @author xuebozhao
 */
public class EnzymePosChr {
    
     EnzymePosChr(String infileS,String cutter,String outfileS) {

        this.SecondFiltering(infileS,cutter,outfileS);
    }

    public void SecondFiltering(String infileS,String cutter,String outfileS) {
        try{
            
           BufferedReader br;            
           br = XueboIOUtils.getTextReader(infileS);
           
            String temp = null;
            StringBuilder temp2 = new StringBuilder();
            int i = 0;
           
            try{
                BufferedWriter bw1 = XueboIOUtils.getTextWriter(outfileS);
                bw1.write("");
                bw1.flush(); 
                bw1.close(); 
            }
            catch (Exception e){
                e.printStackTrace();
            }
            String chr ="";
            
            List <Integer> Pos = null;
            
            while (( temp = br.readLine()) != null) {
               
                if(temp.startsWith(">")){
                    
                    System.out.println("Processing chromosome " + i + "...");
                    if(i > 0){     
                        
                        chr = Integer.toString(i);  
                        
                        Pos = binaryCutter(temp2,cutter);
                        getReducedLibrary(chr,Pos,temp2,outfileS);
                    }
                    
                    i++; 
                    
                    temp2.delete(0, temp2.length());
                }
                else {                    
                    temp2.append(temp);                   
                }
                
            }
            chr = Integer.toString(i);
            Pos = binaryCutter(temp2,cutter);
            getReducedLibrary(chr,Pos,temp2,outfileS);
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
       
    public List<Integer> binaryCutter(StringBuilder inChr,String cutter){
        String query = null;
        List<Integer> pos = new ArrayList<>();
        for(int i = 0 ; i < inChr.length()-cutter.length() ; i++){
            query = inChr.substring(i, i + cutter.length());
            if(query.equals(cutter)){
                pos.add(i);
            }       
        }
        return pos;
    }
    
    public void getReducedLibrary(String Chr,List<Integer> pos,StringBuilder inChr,
            String outfileS){
        try{
        String cutterinChr = null;
        String headcutterinchr = null;
        String headcutterbed = null;
//        BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            for(int i = 0;i < pos.size();i = i+1){
                int aa = Integer.valueOf(pos.get(i)) + 5 ;
                headcutterbed = Chr + "\t" + pos.get(i) + "\t" + aa + "\t0" + "\n";
                getWriteStreamAppend(outfileS,headcutterbed);
            }
        }
        catch (Exception e){
            e.printStackTrace();
        }
        
        
    }

    public static void getWriteStreamAppend(String file, String conent) {                         
        BufferedWriter out = null;                                                   
        try {                                                                        
            out = new BufferedWriter(new OutputStreamWriter(
                    new FileOutputStream(file, true)));                              
            out.write(conent);                                                      
        } 
        catch (Exception e) {                                                     
            e.printStackTrace();                                                    
        } 
        finally {                                                                 
            try {                                                                    
                out.close();                                                        
            } catch (Exception e) {                                               
                e.printStackTrace();                                                
            }                                                                       
        }                                                                           
    }   
    
    
    
}
