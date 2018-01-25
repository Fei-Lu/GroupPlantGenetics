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
public class GenomeDistribute {
    public GenomeDistribute(String inFile,String outFile){
        this.getGenome(inFile,outFile);
    }   
    public void getGenome(String inFile, String outFile){
         try{
            BufferedReader br;
            if(inFile.endsWith("gz"))  br = IOUtils.getTextGzipReader(inFile);
            else  br = IOUtils.getTextReader(inFile);
            String temp = null;
            BufferedWriter bwA = IOUtils.getTextGzipWriter(outFile + "_AA.fa.gz");
            BufferedWriter bwB = IOUtils.getTextGzipWriter(outFile + "_BB.fa.gz");
            BufferedWriter bwAB = IOUtils.getTextGzipWriter(outFile + "_AABB.fa.gz");
            BufferedWriter bwD = IOUtils.getTextGzipWriter(outFile + "_DD.fa.gz");
            BufferedWriter bwABD = IOUtils.getTextGzipWriter(outFile + "_AABBDD.fa.gz");
            boolean aw = false;
            boolean abw = false;
            boolean abdw = false;
            boolean dw = false;
            boolean bw = false;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    aw = false;
                    abw = false;
                    abdw = false;
                    dw = false;
                    bw = false;
                    if(temp.substring(5,6).equals("A")){
                        
                        bwA.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwA.flush();
                        bwAB.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwAB.flush();
                        bwABD.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwABD.flush();
                        aw = true;
                        abw = true;
                        abdw = true;
                    }else if(temp.substring(5,6).equals("B")){
                        bwB.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
                        bwAB.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwAB.flush();
                        bwABD.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwABD.flush();
                        abw = true;
                        abdw = true;
                    }else if(temp.substring(5,6).equals("D")){
                       
                        bwABD.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12) +"\n");
//                        bwABD.flush();
                        bwD.write(">" + temp.substring(4,6)+"_"+temp.substring(11,12)+"\n");
//                        bwD.flush();
                        abdw = true;
                        dw = true;
                    }else{
                        bwA.write(">un" + "\n");
                        bwB.write(">un" + "\n");
//                        bwA.flush();
                        bwAB.write(">un" + "\n");
//                        bwAB.flush();
                        bwABD.write(">un" + "\n");
//                        bwABD.flush();
                        bwD.write(">un" + "\n");
//                        bwD.flush();
                        abdw = true;
                        dw = true;
                        abw = true;
                        aw = true;
                    }
                }else{
                    if(aw){
                        bwA.write(temp+"\n");
//                        bwA.flush();
                    }
                    if(abw){
                        bwAB.write(temp+"\n");
//                        bwAB.flush();
                    }
                    if(abdw){
                        bwABD.write(temp+"\n");
//                        bwABD.flush();
                    }
                    if(dw){
                        bwD.write(temp+"\n");
//                        bwD.flush();
                    }
                    if(bw){
                        bwB.write(temp+"\n");
//                        bwD.flush();
                    }
                }
            }
            bwA.flush();
            bwB.flush();
            bwAB.flush();
            bwABD.flush();
            bwD.flush();
            bwABD.close();
            bwA.close();
            bwB.close();
            bwD.close();
            bwAB.close();
         }
         catch(Exception e){
             e.printStackTrace();
         }
    }
}
