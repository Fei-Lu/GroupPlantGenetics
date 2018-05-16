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
public class ReducedStat {
    public ReducedStat(String inFile,String outFile){
        this.getStat(inFile, outFile);
    };
    public void getStat(String inFile, String outFile){
        BufferedReader br;
        if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
        else  br = YaoIOUtils.getTextReader(inFile);
        String temp;
        String[] tem;
        outFile = outFile + ".stat.gz";
        BufferedWriter bw = YaoIOUtils.getTextGzipWriter(outFile);
        Integer len = 0;
        int i = 0;
        try {
            while ((temp = br.readLine()) != null) {
                if(temp.startsWith(">")){
                   i++;
                   len = Integer.parseInt(temp.split("\t")[2]) - Integer.parseInt(temp.split("\t")[1]);
//                   System.out.println(len);
//                   bw.write(i.toString());
//                   bw.write("\t");
                   bw.write(len.toString());
                   bw.newLine();
                   bw.flush();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(ReducedStat.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    
}
