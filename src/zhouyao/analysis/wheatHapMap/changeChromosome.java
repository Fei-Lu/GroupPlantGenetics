/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class changeChromosome {
    public changeChromosome(String inFile, String outFile,String pos){
        this.getChroChanged(inFile,outFile,pos);
    }
    public void getChroChanged(String inFile,String outFile,String pos){
        try {
            Map bed = new HashMap();
            BufferedReader brbed = YaoIOUtils.getTextReader(pos);
            BufferedReader brvcf = YaoIOUtils.getTextReader(inFile);
            BufferedWriter bwvcf = YaoIOUtils.getTextWriter(outFile);
            String temp_bed = null,temp_vcf =null;
            String[] temp= null;
            while((temp_bed = brbed.readLine())!=null){
               temp = temp_bed.split("\t");
               bed.put(temp[0],temp); 
            }
            String[] line = null;
           
            while((temp_vcf = brvcf.readLine())!=null){
                if(temp_vcf.startsWith("#")){
                    bwvcf.write(temp_vcf+"\n");
                }else{
                    temp = temp_vcf.split("\t");
                    line = (String[]) bed.get(temp[0]);
                    temp[1] =  String.valueOf(Integer.parseInt(temp[1]) +  Integer.valueOf(line[4]));
                    temp[0] =  String.valueOf(line[3]);
                    StringBuilder newline = null;
                    for (int i =0;i<temp.length-2;i++){
                        bwvcf.write(temp[i]+"\t");
                    }
                    bwvcf.write(temp[temp.length-1]+"\n");
                    bwvcf.flush();
                }
            }
            bwvcf.close();
        } catch (IOException ex) {
            Logger.getLogger(changeChromosome.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
        
}
