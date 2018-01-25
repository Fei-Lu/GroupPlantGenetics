/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author yaozhou
 */
public class GenomeStatistic {
    public GenomeStatistic(String inFile, String outFile){
        this.getStatistic(inFile,outFile);
    }
    public void getStatistic(String inFile,String outFile){
        try{
            BufferedReader var;
            if(inFile.endsWith("gz"))  var = IOUtils.getTextGzipReader(inFile);
            else  var = IOUtils.getTextReader(inFile);
            BufferedWriter bw = IOUtils.getTextGzipWriter(outFile + ".stat");
            String content = null;
            int size = -1;
            while((content = var.readLine()) != null){
                if(content.startsWith(">")){
                    System.out.println(content);
                    if(size > -1){
                        bw.write(size +"\n");
//                        bw.flush();
                    }
                    bw.write(content.split(">")[1]+"\t");
//                    bw.flush();
                    size = 0;
                }else{
                    size = size + content.length();
                }
            }
            bw.write(size +"\n");
            bw.flush();
            bw.close();
        }       
        catch (Exception e){
            e.printStackTrace();
        }  
    }
}
