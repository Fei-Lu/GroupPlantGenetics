/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author yaozhou
 */
public class SubstractFromFasta {
    public SubstractFromFasta(String inFile,String names,String outFile){
        this.subNames(inFile, names,outFile);
    }

    private void subNames(String inFile, String names,String outFile) {
        BufferedReader br,bn ;
        BufferedWriter bw,bw1;
        if(inFile.endsWith("gz")){
            br = YaoIOUtils.getTextGzipReader(inFile);
            bw = YaoIOUtils.getTextGzipWriter(outFile);
            bw1 = YaoIOUtils.getTextGzipWriter(outFile+".index");
        }
        else {
            br = YaoIOUtils.getTextReader(inFile);
            bw = YaoIOUtils.getTextWriter(outFile);
        }
        bn = YaoIOUtils.getTextReader(names);
        Set set  =   new  HashSet(); 
        String temp = null,temp1 = null,temp3 = null;
        int i = 0;
       
        try {
            while ((temp = bn.readLine())!=null){
                set.add(temp);
            }
            
            while ((temp = br.readLine())!=null){
                i++;
                if(i%2 == 1){
                    temp1 = temp;
                }else{
                    if(!set.add(temp1)){
//                    System.out.println(temp.substring(temp.length()-104,temp.length()-103));
                        bw.write(temp1+"\n");
                        bw.write(temp+"\n");
                    }
                }
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            System.out.println("Error in Reading: " + inFile);
        }
    }
}
