/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;
import com.koloboke.collect.map.hash.HashIntIntMap;
import com.koloboke.collect.map.hash.HashIntIntMaps;
//import utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import pgl.infra.utils.IOUtils;
/**
 *
 * @author guoyafei
 */
public class MpileupMerge {
    public void merge1() throws IOException{
        String inputfsK="/data2/yafei/Bam_Pileup/Pileup.out/";
        File[] fsK = new File(inputfsK).listFiles();
        List<File> fListK = new ArrayList(Arrays.asList());
        fsK = IOUtils.listFilesEndsWith(fsK, ".pileup");
        fListK=Arrays.asList(fsK);
        HashIntIntMap[] posNum=new HashIntIntMap[42];
        for(int i=0;i<posNum.length;i++){
            posNum[i] = HashIntIntMaps.getDefaultFactory().withDefaultValue(-1).newMutableMap();
        }
        BufferedWriter bw =IOUtils.getTextWriter(new File("/data2/yafei/Bam_Pileup/Pileup.out/merge1.txt").getAbsolutePath());
        fListK.stream().forEach(f -> {
            BufferedReader br =IOUtils.getTextReader(f.getAbsolutePath());
            String temp =null; 
            try{
                while((temp=br.readLine())!=null){
                    int chr=Integer.valueOf(temp.split(" ")[0])-1;
                    int pos = Integer.valueOf(temp.split(" ")[1]);
                    if(posNum[chr].get(pos)!=-1){
                        int num=posNum[chr].get(pos);
                        posNum[chr].put(pos, num+Integer.valueOf(temp.split(" ")[2]));
                    }else{
                        int num=Integer.valueOf(temp.split(" ")[2]);
                        posNum[chr].put(pos,num);
                    }
                }
                br.close();
            }
            catch (Exception e) {
               e.printStackTrace();
            } 
        });
            for(int i=0;i<posNum.length;i++){
                for(int a: posNum[i].keySet()){
                    bw.write(i+1+"\t"+a+"\t"+posNum[i].get(a));
                    bw.newLine();
                }
            } 
    }
}
