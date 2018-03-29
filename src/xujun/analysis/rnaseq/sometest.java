/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import format.genomeAnnotation.GeneFeature;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author xujun
 */
public class sometest {
    public sometest(){
//        this.test1();
        this.test2();
    }
    public void test1(){//Zm00001d031168 对这个超高表达的基因进行统计 第一列是read的长度 第二列是质量值 第三列是这条read的质量均值
        String inputfile="/Users/xujun/Desktop/notSIDname.text";
        String outputfile="/Users/xujun/Desktop/Zm00001d031168.text"; 
        String news=null;
        String phred=null;
        
        try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                phred=news.split("\t")[2];
                int phred1=0;
                if(news.split("\t")[0].equals("Zm00001d031168")){
                    for(int i=0;i<phred.length();i++){
                        phred1=phred1+(int)phred.charAt(i)-33;
                    }
                    bw.write(news.split("\t")[1]/*+"\t"+phred*/+"\t"+phred1/phred.length()); 
                    bw.newLine();
                }
                
             }
             br1.close();
             bw.flush();bw.close();
             
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }           
        
    }
    public void test2(){//找出大于10000的gene
        String inputfile="/Users/xujun/Desktop/notSDIgenecount.text";
        String outputfile="/Users/xujun/Desktop/high.text";
        String news=null;
        String number=null;
        String length=null;
        String name=null;
        String phred=null;
        try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                number=news.split("\t")[1];               
                if(Integer.parseInt(number)==1){
                    name=news.split("\t")[0];
                    length=news.split("\t")[2];
                    phred=news.split("\t")[3];
                    if(length.length()>=3&Integer.parseInt(phred)>25){
                        bw.write(name+"\t"+length+"\t"+phred);             
                        bw.newLine();
                    }
                    
                } 
                    
                }
             br1.close();
             bw.flush();bw.close();                     
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         }           
        
    }
    public void test3(){//
        String inputfile="/Users/xujun/Desktop/notSDIgenecount.text";
        String outputfile="/Users/xujun/Desktop/high.text";
        String news=null;
        String number=null;
        try{
             BufferedReader br1 = utils.IOUtils.getTextReader(inputfile);
             BufferedWriter bw = utils.IOUtils.getTextWriter(outputfile);
             while((news= br1.readLine()) != null){
                number=news.split("\t")[1]; 
                for(int i=1;i<=10;i++){
                    if(1<=Integer.parseInt(number)&Integer.parseInt(number)<100){
                        
                    }
                }
                if(Integer.parseInt(number)==1){
                                 
                    bw.newLine();
                } 
                    
                }
             br1.close();
             bw.flush();bw.close();                     
         }
         catch(Exception ex){
 //            System.out.println((String) w.get(i)+"\t"+intArray[i]+"\t"+averagelength[location]/intArray[i]);
             ex.printStackTrace();
         } 
    }

    public static void main(String[] args){
        new sometest();
    }
}
