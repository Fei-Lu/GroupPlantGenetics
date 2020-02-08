/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.Annotationtwo;
import java.io.BufferedReader;
import pgl.infra.utils.IOUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;

/**
 *
 * @author Xuebo Zhao
 */
public class MeanScore {
    public static void main(String args[]) { 
//        String outfileS = "D:\\zxb\\UniquenessScore\\temp\\Extract_chr10\\maizeAGOv4chr010_CpScore.txt";
//        try {
//            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
//            bw.write("D:\\zxb\\UniquenessScore\\temp\\Extract_chr10\\maize.txt");
//            bw.flush();
//            bw.close();
//        }
//        catch(Exception e) {
//            e.printStackTrace(); 
//        }  

            String infileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/temp_end6_no0.txt";
            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCpGaddpvalue.txt";
            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCHHaddpvalue.txt";
            String outfileS3 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCHGaddpvalue.txt";           
            new BinomialTestMethylation (infileS , outfileS1, outfileS2 , outfileS3 );
            
            
//            String infileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforestMehty_10/siteScoreCpGaddpvalue3X10.txt";
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforestMehty_10/AllsiteScoreCpGaddpvalue3X10.txt";
//            new BinomialTestMethylation (infileS,outfileS);

    }  
}
