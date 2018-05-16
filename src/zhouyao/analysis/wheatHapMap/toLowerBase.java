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
public class toLowerBase {
    public toLowerBase(String inFile,String newIn, String outFile){
        this.getLower(inFile,newIn,outFile);
    }
    public void getLower(String inFile,String newIn,String outFile){
        try {
            BufferedReader br1 = YaoIOUtils.getTextReader(inFile);
            BufferedReader br2 = YaoIOUtils.getTextReader(newIn);
            BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
            String temp = null,temp1 =null;
            
            while((temp=br1.readLine())!=null){
                if(temp.startsWith(">")){
                    temp1 = br2.readLine();
                    bw.write(temp);
                    bw.newLine();
                    bw.flush();
                }else{
                    StringBuilder outString = new StringBuilder();
                    temp1 = br2.readLine();
                    for (int i = 0; i< temp1.length();i++){
                        char chr1 = temp.charAt(i);
                        char chr2 = temp1.charAt(i);
                        if(Character.isUpperCase(chr1)){
                            outString.append(chr2);
                        }else if(Character.isLowerCase(chr1)){
                            outString.append(chr1);
                        }
                    }
                    bw.write(outString.toString());
                    bw.newLine();
                    bw.flush();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(toLowerBase.class.getName()).log(Level.SEVERE, null, ex);
        }

    }
}
