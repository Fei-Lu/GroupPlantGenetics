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
public class ArgsContaining {
    /**
     * @param args the command line arguments
     * -i: input file name 1
     * -j: input file name 2
     * -o: output file name 1
     * -p: output file name 2
     */
    public static void main(String[] args) {
        
        int len = args.length;
        
        String infileS1 = "";
        String infileS2 = "";
        
        String outfileS1 = "";
        String outfileS2 = "";
        
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS1 = args[i+1];
                    infileS2 = args[i+2];
                    i = i +2 ;
                    break;
                    
//                case "-j":
//                    infileS2 = args[i+1];
//                    i++;
//                    break; 
                    
                case "-o":
                    outfileS1 = args[i+1];
                    outfileS2 = args[i+2];
                    i = i +2 ;
                    break;
                    
//                case "-p":
//                    outfileS2 = args[i+1];
//                    i++;
//                    break;
                    
                default:
                    break;
            }
        }
        
        
        if(outfileS1.equals("")){
            String[] temp = infileS1.split("/");
            outfileS1 = temp[temp.length -1].split("\\.")[0] + ".contained.fq.gz";
        }
        
        
        if(outfileS2.equals("")){
            String[] temp = infileS2.split("/");
            outfileS2 = temp[temp.length -1].split("\\.")[0] + ".contained.fq.gz";
        }
        
        new Containing(infileS1 ,infileS2 , outfileS1,outfileS2); 
        
    }
    
}
