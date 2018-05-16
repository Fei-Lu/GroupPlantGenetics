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
public class unmask {
    public unmask(String inFile){
        this.getUnmasked(inFile);
    }
    private void getUnmasked(String inFile){
        try {
            BufferedReader br = YaoIOUtils.getTextGzipReader(inFile);
            String outFile = inFile.replace(".fa.gz", ".unmasked.bed");
            BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
            String temp = null;
            String chr = null;
            int startPos = 0; // using 1-base
            int endPos = 0; //
            boolean newBed = false, writeBed = false;
            int j = 0;
            while ((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    j = 0;
                    chr = temp.replace(">","");
                    newBed = true;
                }else{
                    for (int i = 0; i< temp.length();i++){
                        char chr1 = temp.charAt(i);
                        if(Character.isUpperCase(chr1)){
                            if(chr1!='N'){
                                if(newBed){
                                    startPos = j;
                                    newBed = false;
                                }
                            }else{
                                if(newBed){
                                    startPos = j;
                                    newBed = false;
//                                    if (i==temp.length()-1){
//                                        endPos = j;
//                                        newBed = true;
//                                        bw.write(chr+"\t"+startPos+"\t"+endPos);
//                                        bw.newLine();
//                                        bw.flush();
//                                    }
                                }else{
                                    endPos = j;
                                    newBed = true;
                                    bw.write(chr+"\t"+startPos+"\t"+endPos);
                                    bw.newLine();
                                    bw.flush();  
                                }
                            }
                        }else if(Character.isLowerCase(chr1)){
                           if(!newBed){
                                endPos = j-1;
                                newBed = true;
                                bw.write(chr+"\t"+startPos+"\t"+endPos);
                                bw.newLine();
                                bw.flush();
                            }
                        }
                        j++;
                    }  
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        
    }
}
