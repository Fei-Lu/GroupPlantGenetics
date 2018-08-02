/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author xuebozhao
 */
    
public class GetWheatpseudoGenome {
    
    GetWheatpseudoGenome(String infileS1,String infileS2,String outfileS1) {

        this.readfasta(infileS1,infileS2,outfileS1);
    }

    public void readfasta(String infileS1,String infileS2,String outfileS1) {
        try{
            BufferedReader br;  
            if (infileS1.endsWith(".gz")) {
                br = XueboIOUtils.getTextGzipReader(infileS1);
            }else {
                br = XueboIOUtils.getTextReader(infileS1);
            }
            String temp = null;
            StringBuilder temp2 = new StringBuilder();
            int i = 1;          
            try{
                BufferedWriter bw1 = XueboIOUtils.getTextWriter(outfileS1);
                bw1.write("");
                bw1.flush();                
                bw1.close();                
            }
            catch (Exception e){
                e.printStackTrace();
            }
            String chr ="";
            List <Integer> getPos12 = null;
            while (( temp = br.readLine()) != null) {              
                if(temp.startsWith(">")){
                    System.out.println("Processing chromosome " + i + "...");
                    if(i > 0){   
                        //chr = Integer.toString(i);                          
                        getPos12 = readNumKGFNoUn(infileS2);
                        getReducedgenome(getPos12,temp2,outfileS1);
                    }
                    i++;                 
                    temp2.delete(0, temp2.length());
                }
                else {
                    temp2.append(temp);
                    
                }
                
            }
            //chr = Integer.toString(i);
            getPos12 = readNumKGFNoUn(infileS2);
            getReducedgenome(getPos12,temp2,outfileS1); 
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
       

    public List<Integer> readNumKGFNoUn(String infileS2) throws IOException{
        List<Integer> posgene12 = new ArrayList<>();
        //List<Integer> posChr = new ArrayList<>();
        String tempinfileS2 = null;
        try{
            BufferedReader br2 = XueboIOUtils.getTextReader(infileS2);
            br2.readLine();
            int ii = 0;
            //String tem22 = null;
            //List<Integer> pos12 = new ArrayList<>();
            while((tempinfileS2 = br2.readLine()) != null){
                ++ii;
                if(ii % 10000 == 0){
                    System.out.println("It's time to" + ii);
                }
                String[] teminfileS2 = tempinfileS2.split("\t");
                int posChr1 = Integer.valueOf(teminfileS2[2]);
                posgene12.add(posChr1);
                int posgene1 = Integer.valueOf(teminfileS2[3]);
                posgene12.add(posgene1);
                //posChr.add(posChr1);
                int posgene2 = Integer.valueOf(teminfileS2[4]);
                posgene12.add(posgene2);
                //posChr.add(posChr1);
            }   
           
        }
        catch(Exception e){
            e.printStackTrace();
        }
        return posgene12;
        //return posChr;
    } 
    
    public void getReducedgenome(List<Integer> posgene12 , StringBuilder inChr, String outfileS1){
        try{
        String cuttergene = null;
        String headcutter = null;
        
            for(int i = 0;i < posgene12.size();i = i+3){
                System.out.println(i);
                cuttergene = inChr.substring(posgene12.get(i+1),posgene12.get(i+2));
                
                headcutter =  ">" + posgene12.get(i) +  "\t" + posgene12.get(i+1) + "_" + posgene12.get(i+2) + "\n";
      
                getWriteStreamAppend(outfileS1,headcutter);
                getWriteStreamAppend(outfileS1,cuttergene+"\n");
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