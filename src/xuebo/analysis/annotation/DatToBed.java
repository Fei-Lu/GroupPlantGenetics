/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import pgl.infra.utils.IOUtils;


/**
 *
 * @author xuebozhao
 */
public class DatToBed {
    public DatToBed (String infileS,String outfileS ) {
//        this.readFile(infileS);
      this.readFile(infileS,outfileS);
//      this.writeFile(outfileS);
    }
    
    public void readFile (String infileS,String outfileS) {
        try {
            BufferedReader br;
            br = IOUtils.getTextReader(infileS);
            String temp = null;
//          StringBuilder sb = new StringBuilder(temp);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            while (( temp = br.readLine()) != null) {
                temp = "chr"+ temp;
                bw.write(temp+"\n");
                bw.flush();
            }
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }
//    public void writeFile (String outfileS){
//        try {
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            bw.flush();
//            bw.close();
//        }
//        catch (Exception e) {
//            e.printStackTrace();
//        }
//    }
}

    
