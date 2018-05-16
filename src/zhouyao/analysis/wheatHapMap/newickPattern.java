/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class newickPattern {
   public newickPattern(String inFile){
       this.getPattern(inFile);
   }
   private void getPattern(String inFile){
       BufferedReader br = YaoIOUtils.getTextReader(inFile);
       BufferedWriter bw = YaoIOUtils.getTextWriter(inFile+"_ordered.txt");
       try {
           String a = br.readLine();
           String a1 = a.replace(")", "");
           String a2 = a1.replace("(", "");
           a1 = a2.replace("'", "");
           String[] b = a1.split(",");
           for (int i = 0; i < b.length; i++){
               bw.write(b[i].split(":")[0]+"\n");
           }
           bw.flush();
           bw.close();
       } catch (IOException ex) {
           
       }
       
   }
}
