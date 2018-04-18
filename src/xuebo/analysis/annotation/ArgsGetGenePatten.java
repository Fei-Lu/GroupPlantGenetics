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
public class ArgsGetGenePatten {
     /**
     * @param args the command line arguments
     * -i: input file name 1
     * -j: input file name 2
     * -o: output file name 
     */
    public static void main(String[] args) {
        
        int len = args.length;
        
        String infileS = "";     
        String inFilePos = "";
        String chr = "10-10";
        Integer CutterSize = 50; 

        String outfileS = "";
           
        for (int i = 0; i < len; i++){
            
            if (null != args[i])switch (args[i]) {
                case "-i":
                    infileS = args[i+1];
                    i++;
                    break;
                    
                case "-j":
                    inFilePos = args[i+1];
                    i++;
                    break;  
                                                    
                case "-p":
                    CutterSize = Integer.valueOf(args[i+1]);
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
        
        
//        if(outfileS.equals("")){
//            String[] temp = infileS.split("/");
//            outfileS = temp[temp.length -1].split("\\.")[0] + ".EnzymePosChr.txt";
//        }
        
        new GenePattern(infileS,chr,inFilePos,CutterSize,outfileS); 
        
    }
    
    
}
