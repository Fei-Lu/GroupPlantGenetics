/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author yxh
 */
public class mergeMappedRate {
    
    public mergeMappedRate(){
        this.merge();
       
    }
    
    public void merge(){
        String inputpath = "/Users/yxh/Documents/RNA-seq/test/mergerate";
        String outputpath = "/Users/yxh/Documents/RNA-seq/test/mergerate/result";
        File[] fs = new File(inputpath).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".samLog.final.out");
        List<File> fList = Arrays.asList(fs); 
        List<String> fileList=new ArrayList<>();
        for(int i=0;i<fs.length;i++){
            fileList.add(fList.get(i).getName().replace(".samLog.final.out", ""));
        }
        fileList.stream().forEach((String p) ->{
            try{
                String infileS = new File (inputpath,p+".samLog.final.out").getAbsolutePath();
                BufferedReader br = IOUtils.getTextReader(infileS);
                BufferedWriter [] bw = new BufferedWriter [2];
                for (int i = 0;i<2;i ++){
                    bw [i] = IOUtils.getTextWriter(new File(outputpath,p+".rate").getAbsolutePath());
                }
                int pos = -1;
                String temp = null;
                String rate = null;
                String seq1 = null;String seq2 = null;String seq3 = null;
                String seq4 = null;String seq5 = null;String seq6 = null;
                String seq7 = null;String seq8 = null;String seq9 = null;
                String seq10 = null;
                    while((temp = br.readLine())!=null){
                        seq1 = br.readLine();
                        seq2 = br.readLine();
                        seq3 = br.readLine();
                        seq4 = br.readLine();
                        seq5 = br.readLine();
                        seq6 = br.readLine();
                        seq7 = br.readLine();
                        seq8 = br.readLine();
                        seq9 = br.readLine();
                        seq10 = br.readLine();
                        rate=seq10.split("|")[1];
                        pos =1;
                        bw[pos].write(rate);
                    }
                for(int i=0;i<2;i++){
                      bw[i].flush();bw[i].close();
                }
                br.close();
            }
            catch (Exception e) {
                e.printStackTrace();
                System.exit(1);
            }
        });
        
        
    }
    
    public static void main (String[] args){
        new mergeMappedRate();
    }
}
