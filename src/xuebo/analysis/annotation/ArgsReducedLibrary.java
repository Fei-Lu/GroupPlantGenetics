/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

/**
 *
 * @author xuebozhao
 */
public class ArgsReducedLibrary {
    /**
     * @param args the command line arguments
     * -i: input file name 
     * -p: input RE cutter 1
     * -q: input RE cutter 1
     * -o: output file name 
     */
    public static void main(String[] args) {
        
        int len = args.length;
        
        String infileS = "";
        
        String cutter1 = "";        
        String cutter2 = "";
        
        String outfileS = "";
           
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS = args[i+1];
                    i = i + 4 ;
                    break;
                    
                case "-p":
                    cutter1 = args[i+1];
                    i = i + 4 ;
                    break; 
                    
                case "-q":
                    cutter2 = args[i+1];
                    i = i + 4 ;
                    break; 
                    
                case "-o":
                    outfileS = args[i+1];
                    i = i + 4;
                    break;
                    
                default:
                    break;
            }
        }
        
        
        if(outfileS.equals("")){
            String[] temp = infileS.split("/");
            outfileS = temp[temp.length -1].split("\\.")[0] + ".reducedlibrary.fa.gz";
        }
        
        
        new Containing(infileS ,cutter1 , cutter2, outfileS); 
        
    }
}