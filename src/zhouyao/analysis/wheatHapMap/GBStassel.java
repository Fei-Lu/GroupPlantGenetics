/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class GBStassel {
    String endString ="R1.fastq.gz";
    String keyFile = null;
    public GBStassel(String inFileDir, String KeyFileName){
        File inFile = new File(inFileDir);
//        String inFileDir = inFile.getAbsolutePath();
        this.checkGBSDir(inFileDir);
        File outFileDir = new File(inFileDir+"/GBStassel") ;
        if(!outFileDir.exists()) outFileDir.mkdir();
//        File outFileDir2 = new File(inFileDir+"/GBStassel_R2") ;
        this.keyFile = inFileDir + "/" + KeyFileName;
        this.checkKeyFile(keyFile);
//        this.getKeyFile(inFile, KeyFileName);
//        this.fullName = this.getFullName(inFile);
        this.referAllFile(inFile, endString,outFileDir);
    }
    
    public void referAllFile(File inDir,String endString,File outDir) {
        File[] files = inDir.listFiles();
        Set fileSet = new HashSet();
        if (files == null) {  
          return;  
        }  
        for (File subfile : files) {  
            if (subfile.isDirectory()) {  
                referAllFile(subfile,endString,outDir);  
            } else {  
                if (subfile.getName().endsWith(endString)) {
                    String inFileDir = inDir.getAbsolutePath();
                    if(fileSet.add(inFileDir)){
                        String[] sampleInfo = getInfo(inFileDir);
                        String[] inName = inFileDir.split("/");
                        String fullNamePre = inName[inName.length-2] + inName[inName.length-1];
                        String fullName = fullNamePre.replace("_", "");
                        sampleInfo[3] = fullName;
//                        System.out.println(System.currentTimeMillis());
                        updateFastQ(new File(inFileDir),outDir,sampleInfo);
//                        System.out.println(System.currentTimeMillis());
                        getKeyfile(keyFile,sampleInfo);
                    }
                }  
            }  
        }  
    }  
    
    private void checkGBSDir(String inFile){
        File GBStassel1 = new File(inFile + "_GBStassel_R1");
        if(!GBStassel1.exists()){
            GBStassel1.mkdir();
        }
        File GBStassel2 = new File(inFile + "_GBStassel_R2");
        if(!GBStassel2.exists()){
            GBStassel2.mkdir();
        }
    }
    private void getKeyfile(String keyFile,String[] Info){
        try {
            BufferedWriter bw = new BufferedWriter(new FileWriter(keyFile, true));
            bw.write(Info[1]+"\t"+Info[2]+"\t"+Info[0]+"\t"+Info[3]+"\n");
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            Logger.getLogger(GBStassel.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
    private void checkKeyFile(String keyFile){
        File file=new File(keyFile);
        if(!file.exists()) {      
            try {      
                file.createNewFile();
                BufferedWriter bw = new BufferedWriter(new FileWriter(keyFile, true));
                bw.write("Flowcell\tLane\tBarcode\tFullSampleName\n");
                bw.flush();
                bw.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        } 
    }
    public void updateFastQ(File inFileDir,File outFileDir,String[] sampleInfo){
        File[] files = inFileDir.listFiles();
        String barcode = sampleInfo[0];
        StringBuilder barC = null;
        StringBuilder barS = null;
        StringBuilder out = null;
        try {
            BufferedWriter bw = YaoIOUtils.getTextGzipWriter(outFileDir+"/"+sampleInfo[1]+"_"+sampleInfo[3]+"_"+sampleInfo[2]+"_fastq.gz");
            String temp = null;
            for ( File subfiles : files){
                BufferedReader br = YaoIOUtils.getTextGzipReader(subfiles.toString());
                while((temp =br.readLine())!=null){
                    String[] tem = temp.split(":");
                    barcode = tem[tem.length-1];
                    barC = new StringBuilder(barcode);
                    barS = getCopy(">",barC.length());
                    out = new StringBuilder();
                    bw.write(temp+"\n");
                    temp = br.readLine();
                    barC.append(temp);
                    barC.append("\n");
                    bw.write(barC.toString());
                    temp = br.readLine();
                    out.append(temp);
                    out.append("\n");
                    bw.write(out.toString());
                    temp = br.readLine();
                    barS.append(temp);
                    barS.append("\n");
                    bw.write(barS.toString());
                    bw.flush();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }  
    }
    private String[] getInfo(String inFileDir){
       BufferedReader br = YaoIOUtils.getTextGzipReader(inFileDir+"/R1.fastq.gz");
        String[] out = new String[4];
        String barcode = null;
        try {
            String temp = null;
            while((temp =br.readLine())!=null){
                String[] tem = temp.split(":");
                barcode = tem[tem.length-1];
                if(barcode.contains("N")){
                    temp = br.readLine();
                    temp = br.readLine();
                    temp = br.readLine();
                }else{
                    out[0] = barcode;
                    out[1] = tem[2];//Followcell
                    out[2] = tem[3]; // lane
                    break;
                }
            };
        } catch (IOException ex) {
            ex.printStackTrace();
        }  
        return out;
    }
    public StringBuilder getCopy(String c, int count){
        StringBuilder sb = new StringBuilder(count);
        while(--count>0) sb.append(c);
        return sb;
    }
    public void getCopy(File inDir,File outDir,String R,String fc,String lane,String s){
        String cmd = "cp " + inDir + "/"+R+".fastq.gz "+ outDir + "/"+fc+"_"+s+"_"+lane+"_fastq.gz";
        try {
            Runtime.getRuntime().exec(cmd);
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }
   
}
