/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;

import java.io.File;



/**
 *
 * @author feilu
 */
public class GLoadEntrance {
    public GLoadEntrance() {
        //this.dataOrginazed();
        this.analysisPipeline();
    }
    
    public GLoadEntrance(String a) {
        new PopGenPara(a);
    }
    
    public void analysisPipeline(){
        //new VariantSummary();
        //new Recombination();
        //new TaxaInfo();
        //new Expression();
        //new popGenGroup();
        new PopGenPara();
       // new VCFprocessor();
        //new ItolTree();
        
        
    }

    private void dataOrginazed() {
        new DataOrginazed();
        
    } 
    
    public static void main (String[] args){
        //String inFileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/Parameters_hapScanner.txt";
        //new HapScanner (inFileS);
        System.out.println("this is the entrance of maize GLoadEntrance");
        
        //new GLoadEntrance(args[0]);
        new GLoadEntrance();
        
    }
       
        
       
        
    
}
