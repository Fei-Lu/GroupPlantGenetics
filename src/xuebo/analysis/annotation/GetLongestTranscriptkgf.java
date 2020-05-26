/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author xuebozhao
 */
public class GetLongestTranscriptkgf {
    public GetLongestTranscriptkgf (String infileS,String outfileS ) {
//        this.readFile(infileS);
      this.TTEST(infileS,outfileS);
//      this.writeFile(outfileS);
    }
  public void TTEST(String infileS,String outfileS) {
      try {
        BufferedReader br = IOUtils.getTextReader(infileS);
        BufferedWriter bw = IOUtils.getTextWriter(outfileS);
        String temp = null;
        String temp2 = null;
        String tempbbb = null;
        String[] tembbb = null;
        int num = 0;
        int i = 0;
        temp = br.readLine();
        String tempaa1 = temp;
        String tempaa2 = null;
        String tempaa3 = null;
        String tempaa4 = null;
        String tempaa5 = null;
        String tempaa6 = null;
        
        
        while (( temp = br.readLine()) != null) {
            
            String[] tem = temp.split("\t");
//            num = appearNumber(tem[2], "Transcript");
            if(tem[0].equals("TranscriptNumber")){
                if (tem[1].equals("1") | tem[2].equals("0") ){
                    tempaa1 = temp;
                    tempaa2 = br.readLine();
                    tempaa3 = br.readLine();
                    tempaa4 = br.readLine();
                    tempaa5 = br.readLine();
                    tempaa6 = br.readLine();
                }
                else{
                    int index = Integer.valueOf(tem[2]); 
                    int indexnumber = Integer.valueOf(tem[1]);
                    tempaa1 = temp;
                    int countTranscript = 1;
                    for(i =0 ; countTranscript < index + 1 ;i++){
                        tempbbb  = br.readLine();
                        tembbb = tempbbb.split("\t");
                        if(tembbb[0].equals("Transcript")){
//                            countTranscript++;
                            if(countTranscript == index){
                                tempaa2 = tempbbb;
                                tempaa3 = br.readLine();
                                tempaa4 = br.readLine();
                                tempaa5 = br.readLine();
                                tempaa6 = br.readLine();
                                countTranscript = countTranscript + 1 ;
                            }
                            else{
                                countTranscript = countTranscript + 1 ;
                                
                            }
                            
                        }
                    }
                }
                
                bw.write(tempaa1 + "\n");
                bw.write(tempaa2 + "\n");
                bw.write(tempaa3 + "\n");
                bw.write(tempaa4 + "\n");
                bw.write(tempaa5 + "\n");
                bw.write(tempaa6 + "\n");
            } 
            
        }
        
       bw.flush(); 
       bw.close();   
      }
       catch (Exception e) {
            e.printStackTrace();
        }
   }
    public static int appearNumber(String srcText, String findText) {
        int count = 0;
        Pattern p = Pattern.compile(findText);
        Matcher m = p.matcher(srcText);
        while (m.find()) {
            count++;
        }
        return count;
    }

}
