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
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class BinarySearch {
    
    Integer CutterSize = 50;
    
    public BinarySearch (String infileS1,String infileS2,Integer CutterSize,String outfileS ) {
        
        this.CutterSize = CutterSize;
        List<Integer> posArray = this.getSitePos(infileS2);
        List<Integer> strandArray = this.getSiteStrand(infileS2);
        this.readFileScore(infileS1,posArray,strandArray,CutterSize,outfileS);
        
    }
    
    
    public static void main (String[] args) {
//         String infileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/SitemethyTT10CHH_10.txt";
//          String infileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/GeneFeaturePosCDS.txt";
//          String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/Methy_10/PosCDSCHH.txt";
//          new BinarySearch(infileS1,infileS2,100,outfileS);

            String infileS1 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/test111.txt";
            String infileS2 = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/testPos.txt";
            String outfileS = "/Users/xuebozhao/Documents/LuLab/cpScore/MaizeGeneFeature/MethySingleGeneFeaturePos/test222.txt";
            new BinarySearch(infileS1,infileS2,1,outfileS);
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
    
    private List<Integer> getSitePos(String infileS2) {
        
        BufferedReader br = IOUtils.getTextReader(infileS2);
        String tmp = null;
        String[] temp = null;
        Integer Chr = 0;
        Integer strand  = 0;       
        List<Integer> posArray = new ArrayList<Integer>();   
        List<Integer> strandArray = new ArrayList<Integer>();
        try {
            int j = 0;
            while((tmp = br.readLine())!= null){
                temp = tmp.split("\t");             
                posArray.addAll(getPos(temp[1]));
                //strand = Integer.valueOf(temp[2]);
                //strandArray.add(strand);
                } 
           
            br.close();
        } catch (IOException ex) {
            System.out.println("Died!");
        }
        for(int t = 0;t<posArray.size();t++){
            System.out.println(posArray.get(t));
        }
        return posArray;   
        //return strandArray;
    }
    
    
     private List<Integer> getSiteStrand(String infileS2) {
        
        BufferedReader br = IOUtils.getTextReader(infileS2);
        String tmp = null;
        String[] temp = null;
        Integer Chr = 0;
        Integer strand  = 0;       
        //List<Integer> posArray = new ArrayList<Integer>();   
        List<Integer> strandArray = new ArrayList<Integer>();
        try {
            int j = 0;
            while((tmp = br.readLine())!= null){
                temp = tmp.split("\t");             
                //posArray = getPos(temp[1]);
                strand = Integer.valueOf(temp[2]);
                strandArray.add(strand);
                } 
           
            br.close();
        } catch (IOException ex) {
            System.out.println("Died!");
        }
        for(int t = 0;t<strandArray.size();t++){
            System.out.println(strandArray.get(t));
        }
        //return posArray;   
        return strandArray;
    }
    
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
    
    
    public  void readFileScore(String infileS1,List<Integer>posArray,List<Integer>strandArray,Integer CutterSize,String outfileS){
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS1);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            List<Double> CpScoretemp = new ArrayList<Double>();
            List<Double> CpScoreNew = new ArrayList<Double>();
            List<Double> meanCpScore = new ArrayList<Double>();
            String temp = null;
            int count = 0;
            int arraycount = 0;
            double valuedouble = 0;
            boolean yes = false;
            
            while ((temp = br.readLine()) != null){
                
                valuedouble = Double.parseDouble(temp);
                //arraycount = posArray.size()-2;
                
                //System.out.println(posArray.get(arraycount));
                    if(count==599186){
                        System.out.println();
                    }
                    if( count == posArray.get(arraycount)){  
                        CpScoretemp = new ArrayList<Double>();
                        CpScoretemp.add(valuedouble);
                        yes = true;
                    }               
                    else if(yes){
                        if(count < posArray.get(arraycount+1)){
                            CpScoretemp.add(valuedouble);
                        }
                        else{
                            if((strandArray.get(arraycount/2)) == 0){
                                Collections.reverse(CpScoretemp); 
                            }
//                            CpScoreNew.addAll(CpScoretemp);
                            
                            meanCpScore = getCutterMeanCpscore(CpScoretemp,CutterSize);
                            for (int m = 0; m < meanCpScore.size() -1;m++){
                                bw.write(meanCpScore.get(m)+"\t");
                                    //System.out.println(meanCpScore);
                            }
                            bw.write(meanCpScore.get(meanCpScore.size() -1)+"\n");
                            bw.flush();
//                            CpScoretemp = ;
                            yes = false;
                            arraycount = arraycount + 2;
                            if(arraycount == posArray.size()) {
                                break;
                            }
                        }                       
                    }
                    count = count + 1;
                }
                bw.flush();
                bw.close();                 
            }
        catch(Exception e){
            e.printStackTrace();
        }
        
    }    
             
    
}
