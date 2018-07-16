/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4CandChIA_PET;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author xuebozhao
 */
public class ReducedLibrary {
    
    ReducedLibrary(String infileS,String cutter1,String cutter2, String outfileS1,String outfileS2) {

        this.SecondFiltering(infileS,cutter1,cutter2,outfileS1,outfileS2);
//        this.outputFiles(outfileS1, outfileS2);
    }

    public void SecondFiltering(String infileS, String cutter1,String cutter2,String outfileS1,String outfileS2) {
        try{
            
           BufferedReader br;            
           br = XueboIOUtils.getTextReader(infileS);
           
            String temp = null;
//            String temp1 = null;
//            String temp2 = "";
//            String temp3 = null;
            StringBuilder temp2 = new StringBuilder();
            int i = 0;
           
            try{
                BufferedWriter bw1 = XueboIOUtils.getTextWriter(outfileS1);
                BufferedWriter bw2 = XueboIOUtils.getTextWriter(outfileS2);
                bw1.write("");
                bw2.write("");
                bw1.flush(); 
                bw2.flush(); 
                bw1.close(); 
                bw2.close(); 
            }
            catch (Exception e){
                e.printStackTrace();
            }
            String chr ="";
            List <Integer> getPos = null;
            while (( temp = br.readLine()) != null) {
                
//               System.out.println(tem[0]);
               
                if(temp.startsWith(">")){
                    System.out.println("Processing chromosome " + i + "...");
                    if(i > 0){
//                       String chr = Integer.toString(Integer.valueOf(temp.substring(1,2)) - 1);       
                        chr = Integer.toString(i);  
                        
                        getPos = posSearching(temp2,cutter1,cutter2);
                        getReducedLibrary(chr,getPos,temp2,outfileS1,outfileS2);
//                       getReducedLibraryi = 0;
                    }
                    i++;
//                    String[] tem = temp.split(" ");
//                    temp1 = tem[0];                   
                    temp2.delete(0, temp2.length());
//                    bw.write(tem[0] + "\n");
                }
                else {
//                   while (temp.startsWith("A") || temp.startsWith("T")|| temp.startsWith("C")|| temp.startsWith("G")) {
//                   temp3 = temp2 + temp;
//                   temp2 = temp; 
//                   } 
//                   bw.write(temp3 + "\n");   
//                   temp2 = temp;
//                   temp3 = temp2;
                    temp2.append(temp);
                    
                }
                
            }
            chr = Integer.toString(i);
            getPos = posSearching(temp2,cutter1,cutter2);
            getReducedLibrary(chr,getPos,temp2,outfileS1,outfileS2);
//            bw.write(temp2+"\n");
//            bw.flush();  
//            bw.close(); 
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
    
//    public void binarySearching(String inChr,String cutter1,String cutter2){
//        List<Integer> pos1 = new ArrayList<>(this.binaryCutter(inChr,cutter1));
//        List<Integer> pos2= new ArrayList<>(this.binaryCutter(inChr, cutter2));
//        Integer[] a = (Integer[])pos1.toArray();
//        Integer[] b = pos1.toArray(new Integer[pos1.size()]);
//        int index = Arrays.binarySearch(b ,pos2);
//        index = - index -1;
//        System.out.println(index);
//    }
    
    public List<Integer> posSearching(StringBuilder inChr,String cutter1,String cutter2){
        
        List<Integer> aa = new ArrayList();
//        System.out.println(inChr.length());
        aa = (binaryCutter(inChr,cutter1));
        
        List<Integer> bb = new ArrayList<>(this.binaryCutter(inChr, cutter2));
        
        List<Integer> cc = new ArrayList<>();
        
        List<Integer> getpos = new ArrayList<>();
        
//        String[] both = concact(pos1, pos2);    
//        int c[] = this.contact(pos1, pos2);

        cc.addAll(aa);
        cc.addAll(bb);
        
//        Integer[] aa = pos1.toArray(new Integer[pos1.size()]);
//        Integer[] bb = pos2.toArray(new Integer[pos2.size()]);
//        Integer[] cc = pos3.toArray(new Integer[pos3.size()]);
        Collections.sort(cc);
        
        boolean ispos1 = false;
        boolean ispos2 = false;
        int pos1V = 0;
        int pos2V = 0;
        int j = 0;
        
        for(int i = 0; i<cc.size();i++){
            if(j >= aa.size()){
                aa.add(-1);
            }
//            if(j<aa.size()){
                
                if(cc.get(i) == aa.get(j)){
                    ispos1 = true;  
                    pos1V = cc.get(i);
                    j++;
                    if(ispos1 && ispos2){                       
                        getpos.add(pos2V);
                        getpos.add(pos1V+cutter1.length());
                        ispos2 = false;
                    }
                }
                
                
                else{
                    ispos2 = true;
                    pos2V = cc.get(i);
                    if(ispos1 && ispos2){                       
                      
                            getpos.add(pos1V);
                            getpos.add(pos2V+cutter2.length());
                        
                        ispos1 = false;
                    }
                } 
            }
//        }
        return getpos;     
    }
     

    public void getReducedLibrary(String Chr,List<Integer> getpos,StringBuilder inChr,
            String outfileS1 ,String outfileS2){
        try{
        String cutterinChr = null;
        String headcutterinchr = null;
        String headcutterbed = null;
//        BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            for(int i = 0;i < getpos.size();i = i+2){

                cutterinChr = inChr.substring(getpos.get(i),getpos.get(i+1));
                headcutterinchr = ">" + Chr + "_" + getpos.get(i) + "_" + getpos.get(i+1) + "\n";
                headcutterbed = Chr + "\t" + getpos.get(i) + "\t" + getpos.get(i+1) + "\t0";
                getWriteStreamAppend(outfileS1,headcutterinchr);
                getWriteStreamAppend(outfileS1,cutterinChr+"\n");
                getWriteStreamAppend(outfileS2,headcutterbed +"\n");
//                getWriteStreamAppend(outfileS1,cutterinChr+"\n");
//                bw.write(cutterinChr);
//                bw.write(headcutterinchr);
            }
//            bw.flush();  
//            bw.close(); 
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
