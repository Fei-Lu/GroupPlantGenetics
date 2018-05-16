/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

/**
 *
 * @author yaozhou
 */
public class RemoveDupProtein {
  
    public RemoveDupProtein(String inFile,String outFile){
        this.removeDuplicate(inFile, outFile);
    }

    private void removeDuplicate(String inFile, String outFile) {
        BufferedReader br ;
        BufferedWriter bw,bw1;
        if(inFile.endsWith("gz")){
            br = YaoIOUtils.getTextGzipReader(inFile);
            bw = YaoIOUtils.getTextGzipWriter(outFile);
            bw1 = YaoIOUtils.getTextGzipWriter(outFile+".index");
        }
        else {
            br = YaoIOUtils.getTextReader(inFile);
            bw = YaoIOUtils.getTextWriter(outFile);
            bw1 = YaoIOUtils.getTextWriter(outFile+".index");
        }
        Set set  =   new  HashSet(); 
        String temp = null,temp1 = null,temp3 = null;
        int i = 0;
       
        try {
            while ((temp = br.readLine())!=null){
                i++;
                if(i%2 == 1){
                    temp1 = temp;
                }else{
                    if(set.add(temp)){
//                    System.out.println(temp.substring(temp.length()-104,temp.length()-103));
                        bw.write(temp1+"\n");
                        bw.write(temp+"\n");
                        bw1.write(i+"\n");
                    }
                }
            }
            bw.flush();
            bw1.flush();
            bw.close();
            bw1.close();
        } catch (Exception e) {
            System.out.println("Error in Reading: " + inFile);
        }
    }
    private static void  removeDuplicateWithOrder(List list)   { 
        Set set  =   new  HashSet(); 
        List newList  =   new  ArrayList(); 
        for(Iterator iter  =  list.iterator(); iter.hasNext();)   { 
            Object element  =  iter.next(); 
            if  (set.add(element)) 
                newList.add(element); 
            } 
            list.clear(); 
            list.addAll(newList); 
            System.out.println( " remove duplicate "   +  list); 
        } 
}
