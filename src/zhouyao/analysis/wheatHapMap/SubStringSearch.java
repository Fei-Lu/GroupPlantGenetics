/*
 * This was written for getting the sequence of a genome between two REs;
 * cutter1 and cutter2 are two different REs
 * Only support one genome; so this function must be changed before to testing
 * the whole genome
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class SubStringSearch {
    public SubStringSearch(String genome,String Cutter1,String Cutter2){
        this.getIndex(genome,Cutter1,Cutter2);
    }
    public static List<Integer> getPosition(String genome, String Cutter){
        int LC = Cutter.length();
        List<Integer> pos =  new ArrayList();
        for (int i = LC; i < genome.length()+1;i++){
            if(genome.subSequence(i-LC, i).equals(Cutter)){
                pos.add(i-LC);
            }
        }
        return pos;
    }
    public void getIndex(String genome,String Cutter1,String Cutter2){
        List <Integer> pos1 = getPosition(genome, Cutter1);
        List <Integer> pos2 = getPosition(genome, Cutter2);
//        Object a[] = pos1.toArray();
        List <Integer> index = new ArrayList();
        List <Integer> pos = new ArrayList();
        pos.addAll(pos1);
        pos.addAll(pos2);
        Collections.sort(pos);
        int j = 0;
        boolean isPos1 = false, isPos2 = false;
        int p1 = -1, p2 = -1;
        
        for (int i = 0; i < pos.size(); i++){
            if(j > pos1.size()) pos1.add(-1);
            if(pos.get(i)!=pos1.get(j)){
                isPos2 = true;
                p2 = pos.get(i);
                if (isPos1 && isPos2 ){
                    index.add(p1);
                    index.add(p2);
                    isPos1 = false;
                }
            }else{
                isPos1 = true;
                p1 = pos.get(i);
                j++;
                if (isPos1 && isPos2 ){
                    index.add(p2);
                    index.add(p1);
                    isPos2 = false;
                }
            }
//            if (isPos1 && isPos2 ){
//                index.add(p1);
//                index.add(p2);
//                isPos1 = false;
//                isPos2 = false;
//            }
        }
//        Collections.sort(index);
        BufferedWriter bw = YaoIOUtils.getTextWriter("test.txt");
        try {
            bw.write("");
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(SubStringSearch.class.getName()).log(Level.SEVERE, null, ex);
        }

        for (int q = 0; q < index.size(); q++){
//                BufferedWriter bw = YaoIOUtils.getTextWriter("test.txt");
//            System.out.println(index.get(q));
            int aa = index.get(q);
//                bw.write(aa+"\n");
//                bw.flush();
//                bw.close();
            getAppend("test.txt",aa + "\n");
            getAppend("test.txt",aa + "\n");
        }        
    }  
    
    public  static   void  getAppend(String file, String conent) {  
        BufferedWriter out = null ;  
        try  {  
            out = new  BufferedWriter( new  OutputStreamWriter(  
                    new  FileOutputStream(file,  true )));  
            out.write(conent);  
        } catch  (Exception e) {  
            e.printStackTrace();  
        } finally  {  
            try  {  
                out.close();  
            } catch  (IOException e) {  
                e.printStackTrace();  
            }  
        } 
    }
}
