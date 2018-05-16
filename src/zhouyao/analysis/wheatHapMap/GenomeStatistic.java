/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.List;

/**
 *
 * @author yaozhou
 */
public class GenomeStatistic {
    public GenomeStatistic(String inFile, String outFile,String size){
        List ArraySize = this.getStatistic(inFile);
        this.getStatistic(inFile,outFile,ArraySize,size);
    }
    public GenomeStatistic(String inFile, String outFile){
        this.getStatistic(inFile,outFile);
    }
    public void getStatistic(String inFile,String outFile){
        try{
            BufferedReader var;
            if(inFile.endsWith("gz"))  var = YaoIOUtils.getTextGzipReader(inFile);
            else  var = YaoIOUtils.getTextReader(inFile);
            BufferedWriter bw = YaoIOUtils.getTextGzipWriter(outFile + ".stat");
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
    public List getStatistic(String inFile){
        List sizeArray = new ArrayList();
        try{
            BufferedReader var;
            if(inFile.endsWith("gz"))  var = YaoIOUtils.getTextGzipReader(inFile);
            else  var = YaoIOUtils.getTextReader(inFile);
            String content = null;
            int size = -1;
            while((content = var.readLine()) != null){
                if(content.startsWith(">")){
//                    System.out.println(content);
                    if(size > -1){
                        sizeArray.add(size);
//                        bw.flush();
                    }
                    size = 0;
                }else{
                    size = size + content.length();
                }
            }
            sizeArray.add(size);
            
        }       
        catch (Exception e){
            e.printStackTrace();
        } 
        return sizeArray;
    }
    public void getStatistic(String inFile,String outFile, List ArrySize,String size0){
        List sizeArray = new ArrayList();
        try{
            BufferedReader var;
            if(inFile.endsWith("gz"))  var = YaoIOUtils.getTextGzipReader(inFile);
            else  var = YaoIOUtils.getTextReader(inFile);
            BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
            String content = null;
            int size = 0;
            boolean write = false;
            while((content = var.readLine()) != null){
                if(content.startsWith(">")){
                    Integer size1 = Integer.parseInt(ArrySize.get(size).toString());
                    if(size1 > Integer.parseInt(size0)){
                        size++;
                        System.out.println(size1);
                        bw.write(content);
                        bw.newLine();
                        write = true;
                    }else{
                        write= false;
                    }
                }else{
                    if(write){
                        bw.write(content);
                        bw.newLine();
                        bw.flush();
                    }
                }
            }
            bw.flush();
            bw.close();
        }       
        catch (Exception e){
            e.printStackTrace();
        } 
    }
    
}
