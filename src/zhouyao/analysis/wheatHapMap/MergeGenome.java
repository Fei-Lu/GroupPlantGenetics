/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class MergeGenome {
    MergeGenome(String inFile, String outFile){
        getMerged(inFile,outFile);
    }
    public void getMerged(String inFile, String outFile){
        File InFile = new File (inFile);
        File[] fs = YaoIOUtils.listRecursiveFiles(InFile);
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
        File[] subFs = YaoIOUtils.listFilesEndsWith(fs, ".fna.gz");
        String temp = null;
        BufferedWriter bw = YaoIOUtils.getTextWriter(outFile+"/AllinOne.fa");
        for (File subFile:subFs){
            System.out.println("Processing new genome....");
                BufferedReader br = YaoIOUtils.getTextGzipReader(subFile.toString());
            try {
                while((temp = br.readLine())!=null){
                    bw.write(temp);
                    bw.newLine(); 
                    bw.flush();
                }
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
        try {
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
       
    }
    
}
