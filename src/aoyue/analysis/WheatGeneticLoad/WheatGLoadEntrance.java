/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import format.table.RowTable;
import utils.IOFileFormat;

/**
 *
 * @author Aoyue
 */
public class WheatGLoadEntrance {
    public WheatGLoadEntrance(){
        this.firstProcess();
        
    }
    public void firstProcess(){
        //new MapMake();
        //new Wheat120cleandataProcessor();
        //new Wheat120bamProcessor();
        //new Wheat200cleanDataProcessor();
        new WheatBamDatabase();
        
    }
    
    
    public static void main (String[] args){
        System.out.println("Here is the entrance of wheatGload !");
        System.out.println("I made some revise on itellij ");
        //new WheatGLoadEntrance();

       
        
    }  
    
}
