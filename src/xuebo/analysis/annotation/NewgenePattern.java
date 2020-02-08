/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOUtils;
//import utils.IOUtils;

/**
 *
 * @author xuebozhao
 */

public class NewgenePattern{
    
    Integer CutterSize = 50;
    
    public NewgenePattern (String infileS1,String infileS2,Integer CutterSize,String outfileS ) {
        
        this.CutterSize = CutterSize;
        List<Double> CpScoreHmap = this.readFile(infileS1);
        this.getFullPos(infileS2,CpScoreHmap,CutterSize,outfileS);
        
    }
    
    
    public static void main (String[] args) {
         String infileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/test111.txt";
          String infileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/testPos.txt";
          String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/test222.txt";
          new BinarySearch(infileS1,infileS2,100,outfileS);
    }
    
    
//    public void  test(){
//        
//        int size = 10000;
//        int[] selectedPos = new int[size];
//        String[] result = new String[size];
//        Arrays.sort(selectedPos);
//        String infileS1 = null;
//        String infileS2 = null;
//        String outileS = null;
//        try {
//            BufferedReader br1 = IOUtils.getTextReader(infileS1);
//            BufferedReader br2 = IOUtils.getTextReader(infileS2);
//            BufferedWriter bw = IOUtils.getTextWriter(outileS);
//            String temp = null;
//            int pos = 0;
//            int count = 0;
//            while ((temp = br2.readLine()) != null) {
//                pos++;
//                int index = Arrays.binarySearch(selectedPos, pos);
//                if (index < 0) continue;
//
//                            //bw.write(temp); or
//
//                            result[count] = temp;
//                            count ++;
//            }
//            bw.flush();
//                        bw.close();
//                        br.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//    
//    } 
//    
    
    private List<Integer> getPos(String temp2){
        String[] temp = temp2.split(";");
        Integer poslen = temp.length;
        List<Integer> pos = new ArrayList<Integer>(poslen*2);
        for (int i = 0; i < poslen; i++){
            String[] tmp = temp[i].split(":");
            pos.add(Integer.valueOf(tmp[0]));
            pos.add(Integer.valueOf(tmp[1]));
        }
//        for(int t = 0;t<pos.size();t++){
//            System.out.println(pos.get(t));
//        }
        return pos;
    }
    
    
    public List<Double> readFile(String infileS1){
        
        try{             
            BufferedReader br = IOUtils.getTextReader(infileS1);
            String temp = null;  
            List<Double> CpScoreHmap = new ArrayList<Double>();
            double value_double = 0;
                     
            while((temp = br.readLine())!=null){
                
                if("NA".equals(temp)){
                    value_double = 0;
                    //continue;
                }
                else{
                    value_double = Double.parseDouble(temp);
                }              
                
               CpScoreHmap.add(value_double);
               //System.out.println(CpScoreHmap.size());           
            }          
            return CpScoreHmap;
            
        }
        catch (Exception e) {
            e.printStackTrace();
            return null;
        }
        
    }
    
    
    public List<Double> getPosCpscore(List<Double> CpScorehmap,Integer pos1,Integer pos2,Integer strand){

        Integer pos_start = pos1;
        Integer pos_end = pos2;
        Integer pos_strand = strand;
        
        List<Double> cpScoreArray = new ArrayList<Double>(pos2-pos1+1);
        //List<Double> cpScoreArraytemp = new ArrayList<Double>(pos2-pos1+1);
       
        for(int j = pos_start-1; j < pos_end; j++){            
            cpScoreArray.add(CpScorehmap.get(j));                       
        }    
        //int j = pos_end-1; j < pos_start; j--
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
    
    private void getFullPos(String infileS2,List<Double> CpScorehmap,Integer CutterSize,String outfileS) {
        
        BufferedReader br = IOUtils.getTextReader(infileS2);
        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
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
               List<Double> CpScoreremp = new ArrayList<Double>();
                temp = tmp.split("\t");
                
                //Chr = Integer.valueOf(temp[0]);                
                pos = getPos(temp[1]);
                strand = Integer.valueOf(temp[2]);
                
                for (int i = 0; i < pos.size(); i = i + 2){
                    CpScore.addAll(getPosCpscore(CpScorehmap,pos.get(i),pos.get(i+1),strand));
//                    CpScoreremp = getPosCpscore(CpScorehmap,pos.get(i),pos.get(i+1),strand);
//                    for(int r =0;r < CpScoreremp.size();r++){
//                        //System.out.println(CpScoreremp.get(r));
//                        CpScore.add(CpScoreremp.get(r));
//                    }
                }
                
//                for(int t = 0;t<CpScore.size();t++){
//                    System.out.println(CpScore.get(t));
//                }
                //System.out.println(CpScore.size());
                meanCpScore = getCutterMeanCpscore(CpScore,CutterSize);
//                System.out.println(meanCpScore);
                for (int m = 0; m < meanCpScore.size() -1;m++){
                    bw.write(meanCpScore.get(m)+"\t");
                    //System.out.println(meanCpScore);
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


