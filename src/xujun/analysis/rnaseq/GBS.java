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
    public GBS(){
        this.barcode();
    }
    private void barcode () {
        String sampleInfor="/Users/xujun/Desktop/barcodenews.fq";
        RowTable rt=new RowTable(sampleInfor);
        HashMap hm1=new HashMap();
        HashMap hm2=new HashMap();
        List sample=new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            hm1.put(rt.getCell(i, 2), rt.getCell(i, 3));
            hm2.put(rt.getCell(i, 4), rt.getCell(i, 2));
            sample.add(rt.getCell(i, 2));
        }
        String infileDirS = "/Users/xujun/Desktop/Cleandata/20180427-GBS";
        String outfileDirS = "/Users/xujun/Desktop/Cleandata";
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("-")[0]);
        }
        int[] cnt=new int[sample.size()];
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"-GBS_R1.fq").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"-GBS_R2.fq").getAbsolutePath();
            String outfile = new File (outfileDirS, "count.txt").getAbsolutePath();
            int count=0;int count1=0;
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;           
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);
                BufferedReader br2 = utils.IOUtils.getTextReader(infile2);   
                BufferedWriter bw = utils.IOUtils.getTextWriter(outfile);
                while ((temp1 = br2.readLine()) != null) {
                    count1++;
                    seq1=br2.readLine();
                    if(hm2.get(seq1.substring(0,6))==null){
                        br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                        count++;
                    }else{
                        temp2=br1.readLine();
                        seq2=br1.readLine();
                        if(hm1.get(hm2.get(seq1.substring(0,6)))!=null){
                            int index=sample.indexOf(hm2.get(seq1.substring(0,6)));
                            cnt[index]++;
                        }else{
                            count++;
                        }
                        br1.readLine();br1.readLine();
                    }
                    br2.readLine();br2.readLine();
                }
                for(int i=0;i<cnt.length;i++){
                    bw.write(sample.get(i)+"\t"+cnt[i]);
                    bw.newLine();
                }
                System.out.println(count);
                System.out.println(count1);
                br1.close();br2.close();
                bw.flush();bw.close();
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

        });
    }
}
