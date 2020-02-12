/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import java.io.BufferedWriter;
import java.io.File;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author Aoyue
 */
public class Statistic {
    public Statistic(String inDirS, String outS, String suffix, String unit){
        this.fileSize(inDirS, outS, suffix, unit);
    }
    
    /**
     * @param 
     */
    public void fileSize(String inDirS, String outS, String suffix,String unit){
//        String infileS = "/Users/Aoyue/Documents/maizeGeneticLoad/oriData";
//        String outfileS = "/Users/Aoyue/Documents/statisticTest.txt";
        String infileS = inDirS;
        String outfileS = outS;
        String sizeunit = unit;
        File[] fs = new File(infileS).listFiles();
        //String suffix = ".vcf";
        fs = IOUtils.listFilesEndsWith(fs, suffix);
        ArrayList<String> l = new ArrayList<>();
        for (int i = 0; i < fs.length;i++){
            String name = fs[i].getName().replaceFirst(suffix, "");
            l.add(name);
        }
        Collections.sort(l);
//        for(String a : l){
//            System.out.println(a);
//        }
        
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("ID\t");bw.write(suffix);bw.write("_Size(");bw.write(unit);bw.write(")");bw.newLine();
            DecimalFormat df = new DecimalFormat ("0.00");
            double summary = 0.00;
            double d =0.00;
            for(int i = 0; i < l.size(); i++){
                String name = l.get(i);
                File f = new File(infileS, name + suffix);
                if(sizeunit.equals("TB")){
                    d = (double)f.length()/1024/1024/1024/1024;
                }
                if(sizeunit.equals("GB")){
                    d = (double)f.length()/1024/1024/1024;
                }
                if(sizeunit.equals("MB")){
                    d = (double)f.length()/1024/1024;
                }
                if(sizeunit.equals("KB")){
                    d = (double)f.length()/1024;
                }
                summary = summary + d;
                String size = df.format(d);
                bw.write(name);bw.write("\t");bw.write(size);bw.newLine();
            }
            bw.write("Summary\t");bw.write(String.format("%.2f", summary));bw.newLine();
            bw.flush();bw.close();
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        System.out.println("Calculating the size done !");
    }
}
