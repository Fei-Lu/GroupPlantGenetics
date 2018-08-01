/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author xuebozhao
 */
public class LabNumChrGFF3 {
    public LabNumChrGFF3(String infileS,String outfileS){
        this.ReplaceGFF3ToNum(infileS,outfileS);
    }
    public void ReplaceGFF3ToNum(String infileS,String outfileS){
        try{
            int chrNum = 0;
            int i = 0;
            String temp = null;
            //String tem = null;  
            BufferedReader br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bw = XueboIOUtils.getTextWriter(outfileS);
            while((temp = br.readLine()) != null){
                ++i;
                if(i % 10000 == 0){
                    System.out.println("It's time to" + i);
                }
                String[] tem = temp.split("\t");
                if(temp.startsWith("chrUn")){
                    temp = temp.replaceAll("chrUn", "0");
                    bw.write(temp + "\n");
                }
                if(temp.startsWith("chr1A"))  {
                    if(Integer.valueOf(tem[4]) < 471304005){
                        temp = temp.replaceAll("chr1A", "1");
                    }else{
                        temp = temp.replaceAll("chr1A", "2");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr1B"))  {
                    if(Integer.valueOf(tem[4]) < 438720154){
                        temp = temp.replaceAll("chr1B", "3");
                    }else{
                        temp = temp.replaceAll("chr1B", "4");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr1D"))  {
                    if(Integer.valueOf(tem[4]) < 452179604){
                        temp = temp.replaceAll("chr1D", "5");
                    }else{
                        temp = temp.replaceAll("chr1D", "6");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr2A"))  {
                    if(Integer.valueOf(tem[4]) < 462376173){
                        temp = temp.replaceAll("chr2A", "7");
                    }else{
                        temp = temp.replaceAll("chr2A", "8");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr2B"))  {
                    if(Integer.valueOf(tem[4]) < 453218924){
                        temp = temp.replaceAll("chr2B", "9");
                    }else{
                        temp = temp.replaceAll("chr2B", "10");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr2D"))  {
                    if(Integer.valueOf(tem[4]) < 462216879){
                        temp = temp.replaceAll("chr2D", "11");
                    }else{
                        temp = temp.replaceAll("chr2D", "12");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr3A"))  {
                    if(Integer.valueOf(tem[4]) < 454103970){
                        temp = temp.replaceAll("chr3A", "13");
                    }else{
                        temp = temp.replaceAll("chr3A", "14");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr3B"))  {
                    if(Integer.valueOf(tem[4]) < 448155269){
                        temp = temp.replaceAll("chr3B", "15");
                    }else{
                        temp = temp.replaceAll("chr3B", "16");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr3D"))  {
                    if(Integer.valueOf(tem[4]) < 476235359){
                        temp = temp.replaceAll("chr3D", "17");
                    }else{
                        temp = temp.replaceAll("chr3D", "18");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr4A"))  {
                    if(Integer.valueOf(tem[4]) < 452555092){
                        temp = temp.replaceAll("chr4A", "19");
                    }else{
                        temp = temp.replaceAll("chr4A", "20");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr4B"))  {
                    if(Integer.valueOf(tem[4]) < 451014251){
                        temp = temp.replaceAll("chr4B", "21");
                    }else{
                        temp = temp.replaceAll("chr4B", "22");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr4D"))  {
                    if(Integer.valueOf(tem[4]) < 451004620){
                        temp = temp.replaceAll("chr4D", "23");
                    }else{
                        temp = temp.replaceAll("chr4D", "24");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr5A"))  {
                    if(Integer.valueOf(tem[4]) < 453230519){
                        temp = temp.replaceAll("chr5A", "25");
                    }else{
                        temp = temp.replaceAll("chr5A", "26");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr5B"))  {
                    if(Integer.valueOf(tem[4]) < 451372872){
                        temp = temp.replaceAll("chr5B", "27");
                    }else{
                        temp = temp.replaceAll("chr5B", "28");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr5D"))  {
                    if(Integer.valueOf(tem[4]) < 451901030){
                        temp = temp.replaceAll("chr5D", "29");
                    }else{
                        temp = temp.replaceAll("chr5D", "30");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr6A"))  {
                    if(Integer.valueOf(tem[4]) < 452440856){
                        temp = temp.replaceAll("chr6A", "31");
                    }else{
                        temp = temp.replaceAll("chr6A", "32");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr6B"))  {
                    if(Integer.valueOf(tem[4]) < 452077197){
                        temp = temp.replaceAll("chr6B", "33");
                    }else{
                        temp = temp.replaceAll("chr6B", "34");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr6D"))  {
                    if(Integer.valueOf(tem[4]) < 450509124){
                        temp = temp.replaceAll("chr6D", "35");
                    }else{
                        temp = temp.replaceAll("chr6D", "36");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr7A"))  {
                    if(Integer.valueOf(tem[4]) < 450046986){
                        temp = temp.replaceAll("chr7A", "37");
                    }else{
                        temp = temp.replaceAll("chr7A", "38");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr7B"))  {
                    if(Integer.valueOf(tem[4]) < 453822637){
                        temp = temp.replaceAll("chr7B", "39");
                    }else{
                        temp = temp.replaceAll("chr7B", "40");
                    }
                    bw.write(temp + "\n");
                } 
                if(temp.startsWith("chr7D"))  {
                    if(Integer.valueOf(tem[4]) < 453812268){
                        temp = temp.replaceAll("chr7D", "41");
                    }else{
                        temp = temp.replaceAll("chr7D", "42");
                    }
                    bw.write(temp + "\n");
                }                 
                //bw.write(temp + "\n");
            }
            bw.flush();
            bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    
}
