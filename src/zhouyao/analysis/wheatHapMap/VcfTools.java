/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class VcfTools {
    public VcfTools(String inFile,String outFile){
        this.getABD(inFile,outFile);
    }
    public VcfTools(String inFile,String outFile,boolean depth){
        this.getDepth(inFile,outFile);
    }
    public VcfTools(String inFile,String outFile,String model,String chr){
        int SNPnum = 0; // count the SNP number;
        int SNPnonbi = 0; // count the the non-bialletic SNPs
            
    }
    private void getDepth(String inFile,String outFile){
        try {
            System.out.println("Now analyzing SNP depth...");
            BufferedReader vcf;
            BufferedWriter VcfDepth;
            if(inFile.endsWith("gz"))  vcf = YaoIOUtils.getTextGzipReader(inFile);
            else  vcf = YaoIOUtils.getTextReader(inFile);
            VcfDepth = YaoIOUtils.getTextWriter(outFile+".all");
            String temp = null;
            String[] tem = null;
//            Set depth = new HashSet(); // store the depth catalog
            Integer d = 0;
            int depthArray[][] = new int[101][300];// store the depth for each catalog
            int snp = 0;
            int sampleNum = 0;
            int depthAll[] = new int[501];
            while ((temp = vcf.readLine())!=null){
                if(!temp.startsWith("#")){
                    snp++;
                    if(snp%1000000 == 0) System.out.println("Anlyzing " + snp +"...");
                    tem = temp.split("\t");
                    sampleNum = tem.length - 9;
                    if(tem[4].length()==1){
//                        for (int i = 9 ; i < tem.length; i++){
//                            String re = tem[i];
//                            if(re.split(":").length < 3){
////                                System.out.println(temp);
//                                continue;
//                            }
//                            d = Integer.parseInt(re.split(":")[2]);
//                            if (d > 100) d = 100;
//                            depthArray[d][i-9]++;
//                        } 
                        int dp = Integer.parseInt(tem[7].split("DP=")[1].split(";")[0]);
                        if (dp > 500) dp = 500;
                        depthAll[dp]++;
                   }
                }
            }
//            for (int i = 0; i < 101;i++){
//                VcfDepth.write(Integer.toString(i));
//                for (int j = 0; j < sampleNum;j++){
//                    VcfDepth.write("\t" + depthArray[i][j]);
//                    VcfDepth.flush();
//                }
//                VcfDepth.newLine();
//            }
//            VcfDepth.close();
            for (int i = 0; i < 501; i++){
                VcfDepth.write(Integer.toString(i));
                VcfDepth.write("\t" + depthAll[i]);
                VcfDepth.newLine();
            }
            VcfDepth.flush();
            VcfDepth.close();
        } catch (IOException ex) {
           ex.printStackTrace();
        }
    }
    private void getABD(String inFile, String outFile) {
        BufferedReader vcf;
        BufferedWriter nvcf,bed,stat,Avcf;
        if(inFile.endsWith("gz"))  vcf = YaoIOUtils.getTextGzipReader(inFile);
        else  vcf = YaoIOUtils.getTextReader(inFile);
        nvcf = YaoIOUtils.getTextWriter(outFile + ".all.vcf");
        Avcf = YaoIOUtils.getTextGzipWriter(outFile + ".A.vcf.gz");
        bed = YaoIOUtils.getTextWriter(outFile + ".bed");
        stat = YaoIOUtils.getTextWriter(outFile + ".stat");
        String temp = null;
        String[] temps = null;
        int a = 0, b = 0, d = 0;
        String chr = "0";
        boolean wrt = false;
        try {
            while((temp = vcf.readLine())!=null){
                if(temp.startsWith("#")){
                    nvcf.write(temp + "\n");
                    Avcf.write(temp+"\n");
                }else{
                    temps = temp.split("\t");
                    if(temps[0].substring(1,2).equals("A")){
//                        a++;
                        nvcf.write(temp + '\n');
                        Avcf.write(temp+"\n");
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }else if(temps[0].substring(1,2).equals("B")){
//                        b++;
                        nvcf.write(temp + '\n');
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }else if(temps[0].substring(1,2).equals("D")){
//                        d++;
                        nvcf.write(temp + '\n');
                        bed.write(temps[0]+"\t" + temps[1]+"\t"+temps[1]+"\n");
                        wrt = true;
                    }
                    if(wrt){
                        if(chr.equals("0")){
                            chr = temps[0];
                            a++;
                        }else{
                            if(chr.equals(temps[0])){
                                a++;
                            }else{
                                stat.write(chr+"\t"+a+"\n");
                                a = 1;
                                chr = temps[0];
                            }
                        }
                    }
                }
            }
            stat.write(chr+"\t"+a);
            stat.flush();
            stat.close();
            bed.flush();
            bed.close();
            nvcf.flush();
            nvcf.close();
            Avcf.flush();
            Avcf.close();
        } catch (IOException ex) {
            System.out.println("Reading vcf file failed!");
        }
    }
}
