/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author xuebozhao
 */
public class SplitABD {
    public SplitABD(String infileS,String infileS1,String infileS2,String infileS3,String outfileS1,String outfileS2,String outfileS3){
        //this.GetSplitABD(infileS,infileS1,infileS2,infileS3,outfileS1,outfileS2,outfileS3);
        this.getBlastFileForMcscanX(infileS,infileS1,infileS2,infileS3,outfileS1,outfileS2,outfileS3);
          }       
    public void GetSplitABD(String infileS,String infileS1,String infileS2,String infileS3,String outfileS1,String outfileS2,String outfileS3){
        try{
            int chrNum = 0;
            int i = 0;
            String temp = null;
            String tempAll = null;
            boolean tempA = false;
            boolean tempB = false;
            boolean tempD = false;
            BufferedReader brA = XueboIOUtils.getTextReader(infileS1);
            BufferedReader brB = XueboIOUtils.getTextReader(infileS2);
            BufferedReader brD = XueboIOUtils.getTextReader(infileS3);
            BufferedReader br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bwA = XueboIOUtils.getTextWriter(outfileS1);
            BufferedWriter bwB = XueboIOUtils.getTextWriter(outfileS2);
            BufferedWriter bwD = XueboIOUtils.getTextWriter(outfileS3);
            Set A = new HashSet();
            Set B = new HashSet();
            Set D = new HashSet();
            String chr = null;
            String chrtemp = null;
            while((temp = brA.readLine())!= null){
                A.add(temp);
            }
            while((temp = brB.readLine())!= null){
                B.add(temp);
            }
            while((temp = brD.readLine())!= null){
                D.add(temp);
            }        
            boolean wa = false;
            boolean wb = false;
            boolean wd = false;
            while((tempAll = br.readLine()) != null){
                ++i;
                if(i % 10000 == 0){
                    System.out.println("It's time to" + i);
                }
                if(tempAll.startsWith(">")){
                      chrtemp = tempAll.split(("_"))[0];
                      //System.out.println(chrtemp);
                      chr = chrtemp.substring(1,chrtemp.length());
                      //System.out.println(chr);
                      if(!A.add(chr)){
                        bwA.write(tempAll);
                        bwA.newLine();
                        wa = true;
                        wb = false;
                        wd = false;
                      }else if(!B.add(chr)){
                         A.remove(chr);
                         bwB.write(tempAll);
                         bwB.newLine();
                         wa = false;
                         wb = true;                     
                         wd = false;
                      }else if(!D.add(chr)){
                          A.remove(chr);
                          B.remove(chr);
                          bwD.write(tempAll);
                          bwD.newLine();
                          wa = false;
                          wb = false;
                          wd = true;
                      }else{
                          A.remove(chr);
                          B.remove(chr);
                          D.remove(chr);
                          System.out.println("chromosome unknown!");
                      }
                }else{
                    if(wa){
                        bwA.write(tempAll);
                        bwA.newLine();
                    }else if (wb){
                        bwB.write(tempAll);
                        bwB.newLine();
                    }else if(wd){
                        bwD.write(tempAll);
                        bwD.newLine();
                    }else{
                        continue;
                    }
                }
                  
            }
            bwA.flush();
            bwB.flush();
            bwD.flush();
            bwA.close();
            bwB.close();
            bwD.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
    
    
    
    public void getBlastFileForMcscanX(String infileS,String infileS1,String infileS2,String infileS3,String outfileS1,String outfileS2,String outfileS3){
        try{
            int chrNum = 0;
            int i = 0;
            String temp = null;
            String tempAll = null;
            boolean tempA = false;
            boolean tempB = false;
            boolean tempD = false;
            BufferedReader brA = XueboIOUtils.getTextReader(infileS1);
            BufferedReader brB = XueboIOUtils.getTextReader(infileS2);
            BufferedReader brD = XueboIOUtils.getTextReader(infileS3);
            BufferedReader br = XueboIOUtils.getTextReader(infileS);
            BufferedWriter bwA = XueboIOUtils.getTextWriter(outfileS1);
            BufferedWriter bwB = XueboIOUtils.getTextWriter(outfileS2);
            BufferedWriter bwD = XueboIOUtils.getTextWriter(outfileS3);
            Set A = new HashSet();
            Set B = new HashSet();
            Set D = new HashSet();
            String chr = null;
            String chrtemp = null;
            while((temp = brA.readLine())!= null){
                A.add(temp);
            }
            while((temp = brB.readLine())!= null){
                B.add(temp);
            }
            while((temp = brD.readLine())!= null){
                D.add(temp);
            }        
            boolean wa = false;
            boolean wb = false;
            boolean wd = false;
            while((tempAll = br.readLine()) != null){
                ++i;
                if(i % 10000 == 0){
                    System.out.println("It's time to" + i);
                }
                String[] tem = tempAll.split("\t");
                chr = tem[2];
                if(!A.add(chr)){
                    bwA.write("SubA" + "\t" + tem[2] + "_" + tem[3] + "_" + tem[4] + "_" + tem[1] + "\t" + tem[3] + "\t" + tem[4]);
                    bwA.newLine();
                }else if(!B.add(chr)){
                    A.remove(chr);
                    bwB.write("SubB" + "\t" + tem[2] + "_" + tem[3] + "_" + tem[4] + "_" + tem[1] + "\t" + tem[3] + "\t" + tem[4]);
                    bwB.newLine();
                }else if(!D.add(chr)){
                    A.remove(chr);
                    B.remove(chr);
                    bwD.write("SubD" + "\t" + tem[2] + "_" + tem[3] + "_" + tem[4] + "_" + tem[1] + "\t" + tem[3] + "\t" + tem[4]);
                    bwD.newLine();    
                }else{
                    A.remove(chr);
                    B.remove(chr);
                    D.remove(chr);
                    System.out.println("chromosome unknown!");
                }                 
            }
            bwA.flush();
            bwB.flush();
            bwD.flush();
            bwA.close();
            bwB.close();
            bwD.close();
        }
        catch(Exception e){
            e.printStackTrace();
        }
    }
   


}
