/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import pgl.infra.dna.FastaByte;

/**
 *
 * @author xuebozhao
 */
public class ArgsGetWheatpseudoGenome {
    /**
     * @param args the command line arguments
     * -i: input file name 1
     * -j: input file name 2
     * -o: output file name 1
     */
    public static void main(String[] args) {
        
        int len = args.length;
        
        String infileS = "";     
        FastaByte outfileS1 = null;
        String outfileS2 = "";
        
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS = args[i+1];
                    //infileS2 = args[i+2];
                    i = i + 1;
                    break;
                    
                case "-o":
                    outfileS1 = new FastaByte(args[i+1]);
                    //outfileS2 = args[i+2];
                    i = i + 1 ;
                    break;
                    
                case "-p":
                    outfileS2 = args[i+1];
                    i++;
                    break;
                    
                default:
                    break;
            }
        }
 
        if(outfileS2.equals("")){
            String[] temp = infileS.split("/");
            outfileS2 = temp[temp.length -1].split("\\.")[0] + ".pseudogenome.fa";
        }
        
        new WheatGeneFeature(infileS , outfileS1,outfileS2 ); 
        
    }
    
}
