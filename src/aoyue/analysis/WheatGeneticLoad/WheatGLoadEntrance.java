/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;

/**
 *
 * @author Aoyue
 */
public class WheatGLoadEntrance {
    public WheatGLoadEntrance(){
        //this.firstProcess();
        
        
    }
    public void firstProcess(){
        //new MapMake();
        //new Wheat120cleandataProcessor();
        //new Wheat120bamProcessor();
        //new Wheat200cleanDataProcessor();
        //new WheatBamDatabase();
        new WheatABandDcleandataProcessor();
        //new WheatABDvcfProcessor();
        
    }
    
    
    public static void main (String[] args){
        System.out.println("Here is the entrance of wheatGload !");
        //System.out.println("I made some revise on itellij ");
        new WheatGLoadEntrance();
        
     
        //new WheatGLoadEntrance().test("s", "s");
        
    }  
    
    public void test(String a, String b){
        //a = "aaa";
        //b = "bbb";
        System.out.println("a is " + a);
        System.out.println("b is " + b);
    }
}
