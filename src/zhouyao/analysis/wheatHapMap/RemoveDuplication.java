/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class RemoveDuplication {
    public RemoveDuplication(String inFile,String outFile){
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
        String temp = null,temp3 = null;
        int i=0;
        StringBuilder temp2 = new StringBuilder();
        StringBuilder temp1 = new StringBuilder();
       
        try {
            while ((temp = br.readLine())!=null){
                temp1.delete(0, temp1.length());
                temp2.delete(0, temp2.length());
                String[] tmp = temp.split("\t");
                temp3 = temp.substring(temp.length()-110,temp.length());
                if(i < 1){
                    System.out.println(temp3+"\n");
                }
                if(i>0){
                    for (int j = 0; j<10;j++){
                        temp1.append(tmp[j]+"\t");
                    }
                    temp1.append(tmp[10]);
//                    System.out.println(tmp.length);
                    for(int j = 11; j < tmp.length;j++){
                        if(tmp[j].equals("A")|tmp[j].equals("T")){
                            temp2.append("\t"+tmp[j]);
                        }else if (tmp[j].equals("G")|tmp[j].equals("C")){
                            temp2.append("\t"+tmp[j]);
                        }else{
                            temp2.append("\tN");
                        }
                    }
                }else{
                    temp1.append(temp);
//                    temp2.append("");
                }
                i++;
                if(i==3){
                    int a = 1;
                }
                if(set.add(temp3)){
//                    System.out.println(temp.substring(temp.length()-104,temp.length()-103));
                    bw.write(temp1.toString());
                    bw.write(temp2.toString());
                    bw.write("\n");
                    bw1.write(i+"\n");
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
