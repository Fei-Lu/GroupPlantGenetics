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
public class GenerateHybrid {
    public GenerateHybrid(String inFile, String outFile,boolean hmp){
        if(hmp){
            this.getHybrid_hmp(inFile, outFile);
        }else{
            this.getHybrid(inFile,outFile);
        }
    }

    private void getHybrid(String inFile, String outFile) {
        BufferedReader br = YaoIOUtils.getTextReader(inFile);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
        String tmp = null;
        String[] temp = null;
        StringBuilder to = new StringBuilder();
        try {
            while ((tmp = br.readLine())!=null){
                temp = tmp.split("\t");
                for (int i = 0; i < temp.length-1; i++){
                    for (int j = i+1; j<temp.length;j++){
                        if(temp[i].equals("0") && temp[j].equals("0")){
                            to.append("0\t");
                        }else if(temp[i].equals("0") && temp[j].equals("2")) {
                            to.append("1\t");
                        }else if(temp[i].equals("2") && temp[j].equals("0")) {
                            to.append("1\t");
                        }else if(temp[i].equals("2") && temp[j].equals("2")) {
                            to.append("2\t");
                        }else {
                            System.out.println("Debug: only 0/2 supported");
                        }
                    }
                }
                bw.write(to+"\n");
                to.delete(0,to.length());
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {

        }
    }
    private void getHybrid_hmp(String inFile, String outFile){
        BufferedReader br = YaoIOUtils.getTextReader(inFile);
        BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
        String tmp = null;
        String[] temp = null;
        int line = 0;
        try {
            while ((tmp = br.readLine())!=null){
                line++;
                temp = tmp.split("\t");
                if(line==1){
                    for(int i =0; i <10;i++){
                        bw.write(temp[i]+"\t");
                    }
                    bw.write(temp[10]);
                    for (int j = 11; j < temp.length -1; j++){
                        for(int w = j+1; w<temp.length; w++){
                            bw.write("\t" + temp[j] + "_" + temp[w]);
                        }
                    }
                    bw.write("\n");
                }else{
                   for(int i =0; i <10;i++){
                        bw.write(temp[i]+"\t");
                    }
                    bw.write(temp[10]);
                    for (int j = 11; j < temp.length -1; j++){
                        for(int w = j+1; w<temp.length; w++){
                            bw.write("\t" + getWrite(temp[j], temp[w]));
                        }
                    } 
                    bw.write("\n");
                }
                
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {

        }
    }
    private String getWrite(String a, String b){
        String H = null;
        int A =0, B = 0;
        if("A".equals(a)){
            A = 1;
        }else if("T".equals(a)){
            A = 4;
        }else if("G".equals(a)){
            A = 9;
        }else if("C".equals(a)){
            A = 11;
        }else{
            System.out.println("Only support homozygous");
        }
        
        if("A".equals(b)){
            B = 1;
        }else if("T".equals(b)){
            B = 4;
        }else if("G".equals(b)){
            B = 9;
        }else if("C".equals(b)){
            B = 11;
        }else{
            System.out.println("Only support homozygous");
        }
        
        int C = A + B;
        if(C == 2){
            H = "A";
        }else if (C == 5){
            H = "W";
        }else if (C == 10){
            H = "R";
        }else if (C == 12){
            H = "M";
        }else if (C == 8){
            H = "T";
        }else if (C == 13){
            H = "K";
        }else if (C == 15){
            H = "Y";
        }else if (C == 18){
            H = "G";
        }else if (C == 20){
            H = "S";
        }else if (C == 22){
            H = "C";
        }else{
            System.out.println("Only support homozygous/without missing");
        }
        return H;
    }
}
