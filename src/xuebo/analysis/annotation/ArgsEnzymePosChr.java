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
public class ArgsEnzymePosChr {
    
    public static void main(String[] args) {
        
        int len = args.length;
        
        String infileS = "";       
        String cutter = "";        
        String outfileS = "";
           
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS = args[i+1];
                    i++;
                    break;
                    
                case "-p":
                    cutter = args[i+1];
                    i++;
                    break; 
                    
                case "-o":
                    outfileS = args[i+1];
                    i++ ;
                    break;
                    
                default:
                    break;
            }
        }
        
        
        if(outfileS.equals("")){
            String[] temp = infileS.split("/");
            outfileS = temp[temp.length -1].split("\\.")[0] + ".EnzymePosChr.bed";
        }
        
        new EnzymePosChr(infileS ,cutter, outfileS); 
        
    }
    
}
