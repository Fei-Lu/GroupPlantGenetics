/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;

/**
 *
 * @author yaozhou
 */
public class LogReadStatistic {
    public LogReadStatistic(String inFile){
        this.ReadPrint(inFile);
    }
    public void ReadPrint(String inFile){
        try {
            BufferedReader br ;
            if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
            else  br = YaoIOUtils.getTextReader(inFile);
            String temp;
            int read = 0;
            while ((temp = br.readLine()) != null) {
                String[] tem = temp.split(" ");
                for (int i = 0; i < tem.length; i++)
                    if (tem[i].startsWith("reads")){
                        read = Integer.valueOf(tem[i-1]) + read;
                    }
            }
            System.out.println("read number is: " + read); 
            br.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
    
}
