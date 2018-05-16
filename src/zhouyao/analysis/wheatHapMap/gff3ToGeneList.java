/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class gff3ToGeneList {
    public gff3ToGeneList(String InFile, String OutFile){
        this.getGeneList(InFile,OutFile);
    }

    private void getGeneList(String InFile, String OutFile) {
        BufferedReader br = YaoIOUtils.getTextReader(InFile);
        BufferedWriter bw = YaoIOUtils.getTextWriter(OutFile);
        String tmp = null;
        String[] tmps = null;
        String chr = null,pos1 = null,pos2 = null;
        boolean s = false;
        String gene = null;
        try {
            while((tmp = br.readLine())!=null){
                tmps = tmp.split("\t");
//                if(tmps.length > 8){
//                    if(tmps[1].equals("wareLab")){
//                        chr = tmps[3];
//                    }
//                }
                if(s){
//                    System.out.println(tmps[tmps.length-1]);
                    pos1 = tmps[3];
                    pos2 = tmps[4];
                    gene = tmps[tmps.length-1].split(";")[0].split(":")[1];
                    bw.write(tmps[0] + "\t"+ pos1+ "\t"+pos2+"\t"+ gene +"\t"+ tmps[6]+"\n");
                    s = false;
                }
                if(tmp.startsWith("###")) s = true;
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(gff3ToGeneList.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
}
