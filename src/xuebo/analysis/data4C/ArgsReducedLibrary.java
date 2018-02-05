/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.data4C;

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
        
        String outfileS1 = "";
        String outfileS2 = "";
           
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS = args[i+1];
                    i++;
                    break;
                    
                case "-p":
                    cutter1 = args[i+1];
                    i++;
                    break; 
                    
                case "-q":
                    cutter2 = args[i+1];
                    i++;
                    break; 
                    
                case "-o":
                    outfileS1 = args[i+1];
                    outfileS2 = args[i+2];
                    i = i +2 ;
                    break;
                    
                default:
                    break;
            }
        }
        
        
        if(outfileS1.equals("")){
            String[] temp = infileS.split("/");
            outfileS1 = temp[temp.length -1].split("\\.")[0] + ".reducedlibrary.fas";
        }
        
        if(outfileS2.equals("")){
            String[] temp = infileS.split("/");
            outfileS2 = temp[temp.length -1].split("\\.")[0] + ".libraryAll.bed";
        }
        new ReducedLibrary(infileS ,cutter1 , cutter2, outfileS1,outfileS2); 
        
    }
}