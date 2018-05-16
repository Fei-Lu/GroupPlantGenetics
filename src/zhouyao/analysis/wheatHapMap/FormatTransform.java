/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 * Transform data  among plink, blink(numeric) and hapmap format
 * --file: input file name without suffix
 * --make-plink: transform blink format to plink format
 *           *.plk.ped and *.plk.map file will be generated
 * --make-blink: transform hmp, plink or vcf format to blink format
 *           *.blk.dat and *.blk.map file will be generated
 * --out: outFile name without suffix, suffix will be added automatically:
 *        
 * @author yaozhou (yao.zhou@genetics.ac.cn)
 */
public class FormatTransform {
    
    public FormatTransform(String file, boolean plink, boolean blink,String outFile, boolean gzip){
       this.getDecode(file,plink,blink,outFile,gzip); 
    }

    private void getDecode(String file, boolean plink, boolean blink, String outFile,boolean gzip) {
        if(plink){
            BufferedReader br;
            String plkg = outFile + ".plk.ped";
            String plkm = outFile + ".plk.map";
            try{
                if(gzip) br = YaoIOUtils.getTextGzipReader(file);
                else br = YaoIOUtils.getTextReader(file);
                BufferedWriter bwgeno = YaoIOUtils.getTextWriter(plkg);
                BufferedWriter bwmap = YaoIOUtils.getTextWriter(plkm);
            }catch (Exception e){
                System.out.println("Only blink format supported to plink!");
            }
        }else if(blink){
            String blkg = outFile + ".blk.dat";
            String blkm = outFile + ".blk.map";
            BufferedReader br;
            try{
                if(gzip) br = YaoIOUtils.getTextGzipReader(file);
                else br = YaoIOUtils.getTextReader(file);
                BufferedWriter bwgeno = YaoIOUtils.getTextWriter(blkg);
                BufferedWriter bwmap = YaoIOUtils.getTextWriter(blkm);
                String o = "NA";
                String temp = null;
                while((temp = br.readLine())!= null){
                    String[] tem = temp.split("\t");
                    if(tem[0].startsWith("rs")){
                        System.out.println("Converting hapmap to blink...");
                        bwmap.write("rs\tChr\tpos\n");
                    }else{
                        String a = tem[1].split("/")[0];
                        String b = tem[1].split("/")[1];
                        for(int i = 11;i < tem.length;i++){
                            if (tem[i].equals(a)){
                                o = "0";
                            }else if (tem[i].equals(b)){
                                o = "2";
                            }else if(tem[i].equals("N")){
                                o = "NA";
                            }else{
                                o = "1";
                            }
                            if( i == 11){
                                bwgeno.write(o);
                            }else{
                                bwgeno.write("\t" + o);
                            }
                        }
                        bwgeno.write("\n");
                        bwmap.write(tem[0]+"\t"+tem[2]+"\t"+tem[3]+"\n");
                    }
                }
                bwgeno.flush();
                bwmap.flush();
                bwmap.close();
                bwgeno.close();
            }catch (Exception e){
                System.out.println("Can't find input file");
            }
        }else{
            System.out.println("No output format!");
        }
    }
}
