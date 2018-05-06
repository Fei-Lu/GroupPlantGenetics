/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import utils.IOUtils;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author xujun
 */
public class GBS {
    private void barcode () {
        String sampleInfor="";
        RowTable rt=new RowTable(sampleInfor);
        HashMap hm1=new HashMap();
        HashMap hm2=new HashMap();
        List sample=new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            hm1.put(rt.getCell(i, 3), rt.getCell(i, 2));
            hm2.put(rt.getCell(i, 4), rt.getCell(i, 2));
            sample.add(rt.getCell(i, 2));
        }
        String infileDirS = "";
        String outfileDirS = "";
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        int[] cnt=new int[sample.size()];
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_1.clean.fq").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_2.clean.fq").getAbsolutePath();
            String outfile = new File (outfileDirS, "count.txt").getAbsolutePath();
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;           
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);
                BufferedReader br2 = utils.IOUtils.getTextReader(infile2);   
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while ((temp1 = br1.readLine()) != null) {
                    seq1=br1.readLine();
                    if(hm1.get(seq1.substring(0,6)).equals(null)){
                        br2.readLine();br2.readLine();br2.readLine();br2.readLine();
                    }else{
                        temp2=br2.readLine();
                        seq2=br2.readLine();
                        if(hm2.get(seq2.substring(0,6)).equals(hm1.get(seq1.substring(0,6)))){
                            int index=sample.indexOf(hm2.get(seq2.substring(0,6)));
                            cnt[index]++;
                        }
                        br2.readLine();br2.readLine();
                    }
                    br1.readLine();br1.readLine();
                }
                for(int i=0;i<cnt.length;i++){
                    bw.write(sample.get(i)+"\t"+cnt[i]);
                    bw.newLine();
                }
                br1.close();br2.close();
                bw.flush();bw.close();
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
}
