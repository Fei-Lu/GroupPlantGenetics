/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.sv;

import analysis.maize2k.HapScanner;

/**
 *
 * @author feilu
 */
public class SVEntraince {
    public SVEntraince() {
        this.dataOrginazed();
        //new TaxaDiversitycp();
        
        //String inFileS = "/Users/Aoyue/Documents/Data/pipeline/hapScanner/Parameters_hapScanner.txt";
        //new HapScanner (inFileS);
        
    }

    private void dataOrginazed() {
        new DataOrginazed();
        //new HapMapTaxaProcessorcp ();
    } 
    
    public static void main (String[] args){
        //new SVEntraince();
        new FastqQualitycp();
        
    }  
}
