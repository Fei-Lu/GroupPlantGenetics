/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.io.BufferedReader;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;
import utils.IOUtils;

/**
 *
 * @author kanglipeng
 */
public class repeatMaskDepth {
    String maskedGenomeFile=null;

    public repeatMaskDepth(String infiles){
        this.parseParameters(infiles);
    
    }
    public repeatMaskDepth(String maskedGenomeFile,String geneName){
    this.getMaskedPercent(maskedGenomeFile);
    }
    
    public void getMaskedPercent(String maskedGenomeFile){
          try {
              File maskedGenome=new File(maskedGenomeFile);
              String genomeName = maskedGenome.getName();
           BufferedReader br = IOUtils.getTextReader(maskedGenomeFile);
           String temp = null;
           int maskedBaseNum=0;
           int allBaseNum=0;
            while ((temp = br.readLine()) != null) {
                if(temp.startsWith(">"))continue;
                for(int i=0;i<=temp.length()-1; i++){
                    char a=temp.charAt(i);
                    allBaseNum++;
                    if (!Character.isUpperCase(a)){
                    maskedBaseNum++;
                    }

                }
                
            }
            DecimalFormat df=new DecimalFormat("0.00");
            System.out.println(100*(float)maskedBaseNum/allBaseNum +"of"+genomeName+""+"is soft-masked");
          }
            catch (Exception e) {
            System.out.println("Error in getting maskedGenomeDepth " );
        }
    
    
    }
    private void parseParameters(String infileS) {
        List<String> pLineList = new ArrayList<>();
        try {
            BufferedReader br = IOUtils.getTextReader(infileS);
            String temp = null;
            boolean ifOut = false;
            if (!(temp = br.readLine()).equals("calling maf cover gene")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Author: clip")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Email: lipenkang@163.com")) {
                ifOut = true;
            }
            if (!(temp = br.readLine()).equals("Homepage: http://plantgeneticslab.weebly.com/")) {
                ifOut = true;
            }
            if (ifOut) {
                System.out.println("Thanks for using TEP.");
                System.out.println("Please keep the authorship in the parameter file. Program stops.");
                System.exit(0);
            }
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith("#")) {
                    continue;
                }
                if (temp.isEmpty()) {
                    continue;
                }
                pLineList.add(temp);
            }
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
        this.maskedGenomeFile = pLineList.get(0);
       }
    public static void main(String[] args) {
        String geneName=null;
      repeatMaskDepth le=  new repeatMaskDepth(args[0]);
      new repeatMaskDepth(le.maskedGenomeFile, geneName);
        
    }
    
    
    
}
