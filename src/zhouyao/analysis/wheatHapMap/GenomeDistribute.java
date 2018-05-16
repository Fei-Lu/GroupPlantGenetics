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
    public GenomeDistribute(String inFile,String outFile,boolean byChr){
        this.getChr(inFile,outFile);
    }   
    public GenomeDistribute(String inFile,String outFile,String chrNUM){
        this.getChrNum(inFile,outFile,chrNUM);
    }
    private void getGenome(String inFile, String outFile){
         try{
            BufferedReader br;
            if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
            else  br = YaoIOUtils.getTextReader(inFile);
            String temp = null;
            BufferedWriter bwA = YaoIOUtils.getTextGzipWriter(outFile + "_AA.fa.gz");
            BufferedWriter bwB = YaoIOUtils.getTextGzipWriter(outFile + "_BB.fa.gz");
            BufferedWriter bwAB = YaoIOUtils.getTextGzipWriter(outFile + "_AABB.fa.gz");
            BufferedWriter bwD = YaoIOUtils.getTextGzipWriter(outFile + "_DD.fa.gz");
            BufferedWriter bwABD = YaoIOUtils.getTextGzipWriter(outFile + "_AABBDD.fa.gz");
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
    private void getChr(String inFile, String outFile){
         try{
            BufferedReader br;
            if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
            else  br = YaoIOUtils.getTextReader(inFile);
            String temp = null;
            int count = 0;
            boolean write = false;
            String name = null;
            BufferedWriter bw = null;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    System.out.println("Processing: " + temp + "\n");
                    if(write){
                        bw.flush();
                        bw.close();
                    }else{
                        write = true;
                    }
                    name = temp.replace(">","");
                    String[] names = name.split(" ");
                    name = names[0];
                    bw = YaoIOUtils.getTextWriter(outFile +"/"+name+".fa");
                    bw.write(temp+"\n");
                }else{
                    bw.write(temp+"\n");
                } 
            }
            bw.flush();
            bw.close();
         }
         catch(Exception e){
             e.printStackTrace();
         }
    }
    private void getChrNum(String inFile, String outFile,String chrNum){
         try{
            BufferedReader br;
            if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
            else  br = YaoIOUtils.getTextReader(inFile);
            String temp = null;
            int count = 0;
            boolean write = false;
            String name = null;
            BufferedWriter bw = null;
            Integer nm = Integer.parseInt(chrNum);
            int i = 0;
            while((temp = br.readLine())!=null){
                if(temp.startsWith(">")){
                    i++;
                    if(i>nm) break;
                    System.out.println("Processing: " + temp + "\n");
                    if(write){
                        bw.flush();
                        bw.close();
                    }else{
                        write = true;
                    }
                    name = temp.replace(">","");
                    String[] names = name.split(" ");
                    name = names[0];
                    bw = YaoIOUtils.getTextWriter(outFile +"/"+name+".fa");
                    bw.write(temp+"\n");
                }else{
                    bw.write(temp+"\n");
                } 
            }
            bw.flush();
            bw.close();
         }
         catch(Exception e){
             e.printStackTrace();
         }
    }
}
