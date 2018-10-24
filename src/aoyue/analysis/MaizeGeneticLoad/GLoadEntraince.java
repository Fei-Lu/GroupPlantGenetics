/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.MaizeGeneticLoad;



/**
 *
 * @author feilu
 */
public class GLoadEntraince {
    public GLoadEntraince() {
        //this.dataOrginazed();
        this.analysisPipeline();;
        
        
        
        
    }
    public void analysisPipeline(){
        //new VariantSummary();
        //new Recombination();
        //new TaxaInfo();
        new Expression();
        
    }

    private void dataOrginazed() {
        new DataOrginazed();
        
    } 
    
    public static void main (String[] args){
        //String inFileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/Parameters_hapScanner.txt";
        //new HapScanner (inFileS);
        new GLoadEntraince();
       //System.out.println("2222");
        
       
        
    }  
}
