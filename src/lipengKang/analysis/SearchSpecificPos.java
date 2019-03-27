/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package lipengKang.analysis;

import java.io.BufferedReader;
import zhouyao.analysis.wheatHapMap.YaoIOUtils;

/**
 *
 * @author kanglipeng
 /逢10无法count*/
public class SearchSpecificPos {

    public SearchSpecificPos(String genomeDir) {
        this.getSpecificSegment(genomeDir);
    }

    public void getSpecificSegment(String genomeDir) {
        BufferedReader br;

        if (genomeDir.endsWith("gz")) {
            br = YaoIOUtils.getTextGzipReader(genomeDir);

        } else {
            br = YaoIOUtils.getTextReader(genomeDir);

        }
        String temp = null;
        String[] tem = null;
        int pos = 0;
        int m=70557;//pay attention this number is 1 based
        
        try {
            while ((temp = br.readLine()) != null) {
                if (temp.startsWith(">")) {
                    continue;
                }
                if (pos==(int)(60*Math.floor(m/60))) {
                    //58508	58768//s wheatA.chr1    148052 2251 + 594102056 gcaaggattgacaggct-----
                    
                    
                   /* AGCAAGGATTGACAGGCTCTCTTTCTTGA
CAACACGGGGAAACTTACCAGGTCTAGACATAGCAAGGATTGACAGGCTCTCTTTCTTGA
TTCTATGGGTGGTGGTGCATGGCCGTTCTTAGTTGGTGGAGCGATTTGTCTGGcTTAATTC
148080
60*/
                   /*agcaaggattgacaggctctctttcttga
caacacggggaaacttaccaggtctagacatagcaaggattgacaggctctctttcttga
ttctatgggtggtggtgcatggccgttcttagttggtggagcgatttgtctggttaattc
148080
60*/
                    System.out.println(temp.substring((int) (m-60*Math.floor(m/60)-1)));
                    System.out.println(temp);
                    temp = br.readLine();
                    pos = pos + temp.length();
                    if (pos >= m) {
                        
                        System.out.println(temp);
                        System.out.println(pos);
                        temp=br.readLine();
                        System.out.println(temp);
                        System.out.println(pos+60);
                         temp=br.readLine();
                        System.out.println(temp);
                        System.out.println(pos+120); temp=br.readLine();
                        System.out.println(temp);
                        System.out.println(pos+180); temp=br.readLine();
                        System.out.println(temp);
                        System.out.println(pos+240); temp=br.readLine();
                        System.out.println(temp);
                        System.out.println(pos+300);
                        System.out.println(temp.length());
                        break;
                    }
                }else{pos=pos+60;}
            }
            //123 4 5 6
            //789101112
            //6- (12-10)-1=3
           br.close();
        } catch (Exception e) {
            e.printStackTrace();
            System.exit(1);
        }
    }
    public static void main(String[] args){
    new SearchSpecificPos(args[0]);
    }

}

