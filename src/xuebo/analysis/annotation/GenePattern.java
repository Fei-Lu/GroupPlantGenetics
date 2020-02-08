/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import pgl.infra.utils.IOUtils;
//import utils.IOUtils;

/**
 *
 * @author xuebozhao
 */

public class GenePattern {
    
    Integer CutterSize = 50;
    
    public GenePattern (String infileS,String chr,String inFilePos,Integer CutterSize,String outFileS ) {
        
        this.CutterSize = CutterSize;

        
        HashMap CpScorehmap = this.getAll(infileS, chr);
        
        this.getFullPos(inFilePos,CpScorehmap,CutterSize,outFileS);
        
//        List chrt = this.getPosCpscore(CpScorehmap, 6, 2, 6);
//        for (int i = 0;  i < 4;i++){
//           System.out.println(CpScorehmap.get(3+"_"+Integer.toString(i)));
//        }
//        for (int i = 0;  i < chrt.size();i++){
//           System.out.println(chrt.get(i));
//        }
//        List ta = getCutterMeanCpscore(chrt,4);
//        for (int i = 0;  i < ta.size();i++){
//           System.out.println(ta.get(i));
//        }
        
    }
    
    
  
    private HashMap<String,Double> getAll(String inFileS,String chrinfo){
        
        HashMap<String, Double> cpScoreHmAll = new HashMap<String,Double>();        
        List<Integer> chr = getChr(chrinfo);
        
        for (int i = chr.get(0); i < chr.get(1) + 1;i++){
            
            String chrString = inFileS + "_" + i + ".txt";
//            String chrString = inFileS;

            HashMap<String, Double> cpScoreHm = readFile(chrString,i);
            
            cpScoreHmAll.putAll(cpScoreHm);
//            for (int j = 0;  j < cpScoreHm.size();j++){
//                System.out.println(cpScoreHm.get(chr+"_"+Integer.toString(j)));
//            }
//            for (int j = 0;  j < cpScoreHmAll.size();j++){
//                System.out.println(cpScoreHmAll.get(1+"_"+Integer.toString(j)));
//            }
        }
        return cpScoreHmAll;
    }
    
    private List<Integer> getChr(String chrinfo){
        String[] temp = chrinfo.split("-");
        Integer chrS = Integer.parseInt(temp[0]);
        Integer chrE = Integer.parseInt(temp[1]);
        List<Integer> chrR = new ArrayList<Integer>(2);
        chrR.add(chrS);
        chrR.add(chrE);
        return(chrR);
    }   
    
    
    
    private List<Integer> getPos(String posinfo){
        String[] temp = posinfo.split(";");
        Integer poslen = temp.length;
        List<Integer> pos = new ArrayList<Integer>(poslen*2);
        for (int i = 0; i < poslen; i++){
            String[] tmp = temp[i].split(":");
            pos.add(Integer.valueOf(tmp[0]));
            pos.add(Integer.valueOf(tmp[1]));
        }
        return pos;
    }
    
    
    
   
    public HashMap<String, Double> readFile(String infileS,Integer chr){
        
        try{             
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;  
            HashMap<String, Double> CpScoreHmap = new HashMap<String, Double>();
            String kp = null;
            double value_double = 0;
            int i = -1;
            
            while((temp = br.readLine())!=null){
                
                i++;
                
                if("NA".equals(temp)){
                    //value_double = 0;
                    continue;
                }
                else{
                    value_double = Double.parseDouble(temp);
                }              
                
                kp = chr +"_"+ Integer.toString(i);
                CpScoreHmap.put(kp,value_double);
            }
//            for (int j = 0;  j < CpScoreHmap.size();j++){
//                System.out.println(CpScoreHmap.get(chr+"_"+Integer.toString(j)));
//            }
//            br.flush();
//            br.close(); 
            return CpScoreHmap; 
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
    }
 
    
      
    public List<Double> getPosCpscore(HashMap<Integer,Double> CpScorehmap,Integer chr,Integer pos1,Integer pos2,Integer strand){
//        Integer tmp;
//        if(pos1 > pos2){
//            tmp = pos1;
//            pos1 = pos2;
//            pos2 = tmp;
//        }
        Integer pos_start = pos1;
        Integer pos_end = pos2;
        Integer pos_strand = strand;
        
        List<Double> cpScoreArray = new ArrayList<Double>(pos2-pos1+1);
       
        for(int j = pos_start-1; j < pos_end; j++){
            
            cpScoreArray.add(CpScorehmap.get(chr + "_" +j));
            
        }
        
        if (pos_strand == 0){
            
            Collections.reverse(cpScoreArray);
         
        }
        
        return cpScoreArray;
    } 
    
    
    
    public List<Double> getCutterMeanCpscore(List<Double>cpScoreArray,Integer CutterSize){
        
//        int aa = cpScoreArray.size();
        List<Double> meancpScoreArray = new ArrayList<Double>(CutterSize);
        
        double bb = 0;
        int m = 1;
        int meanSize = 0;
        
//        if (cpScoreArray.size()/CutterSize == 0){
//            meanSize = cpScoreArray.size() / CutterSize;
//        }
//        else{
//             meanSize = cpScoreArray.size() / CutterSize + 1;
//        }
        meanSize = cpScoreArray.size()/CutterSize;
       
        double sum = 0;
        if ( cpScoreArray.size() > CutterSize ){ 
            for( int i = 1;i < meanSize * (CutterSize - 1)+1;i++){
//                int j = i / meanSize;
                sum = sum + cpScoreArray.get(i-1);
                if(i % meanSize == 0){
                    sum = sum/meanSize;
                    meancpScoreArray.add(sum);
                    sum  = 0 ;
                }
            }
            for(int j = meanSize * (CutterSize - 1);j < cpScoreArray.size();j++ ){
               sum = sum + cpScoreArray.get(j);
            }
            sum = sum /(cpScoreArray.size() - meanSize * (CutterSize - 1));
            meancpScoreArray.add(sum);
  
        }
        else{
            for(int k = 0; k < cpScoreArray.size();k++){
                bb = bb + cpScoreArray.get(k);
            }
            meancpScoreArray.add(bb/cpScoreArray.size());                  
        }
        return meancpScoreArray;  
    }
    
    private void getFullPos(String inFilePos,HashMap<Integer,Double> CpScorehmap,Integer CutterSize,String outFileS) {
        
        BufferedReader br = IOUtils.getTextReader(inFilePos);
        BufferedWriter bw = IOUtils.getTextWriter(outFileS);
        String tmp = null;
        String[] temp = null;
        Integer Chr = 0;
        Integer strand  = 0;
        
        List<Integer> pos = new ArrayList<Integer>();
        
        
        List<Double> meanCpScore = new ArrayList<Double>();
        
        try {
            int j = 0;
            while((tmp = br.readLine())!= null){
               List<Double> CpScore = new ArrayList<Double>();
                temp = tmp.split("\t");
                
                Chr = Integer.valueOf(temp[0]);                
                pos = getPos(temp[1]);
                strand = Integer.valueOf(temp[2]);
                
//                 j++;
//                if(j>48){
//                    System.out.printf("biubiubiu");
//                    for(int c =0; c<pos.size();c++){
//                        System.out.println(pos.get(c));
//                    }
//                }
                
                for (int i = 0; i < pos.size(); i = i + 2){
                    CpScore.addAll(getPosCpscore(CpScorehmap,Chr,pos.get(i),pos.get(i+1),strand)) ;
                }
                
                meanCpScore = getCutterMeanCpscore(CpScore,CutterSize);
//                for(int i= 0; i<10;i++){
//                    System.out.println(meanCpScore.get(i));
//                }
                for (int i = 0; i < meanCpScore.size() -1;i++){
                    bw.write(meanCpScore.get(i)+"\t");
                }
                bw.write(meanCpScore.get(meanCpScore.size() -1)+"\n");
                bw.flush();
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            System.out.println("Died!");
        }
             
    }

   
}

