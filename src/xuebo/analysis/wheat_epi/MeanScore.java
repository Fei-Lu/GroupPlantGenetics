/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;
import pgl.infra.dna.FastaByte;
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

//            String infileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/temp_end6_no0.txt";
//            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCpGaddpvalue.txt";
//            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCHHaddpvalue.txt";
//            String outfileS3 = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/V4/siteScoreCHGaddpvalue.txt";           
//            new BinomialTestMethylation (infileS , outfileS1, outfileS2 , outfileS3 );
            
//            String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/testwheat.gff3";
//            FastaByte outfileS1 = new FastaByte("/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/testFasta50001.fa");
//            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/123.txt";
//            new WheatGeneFeature (infileS , outfileS1,outfileS2);


//            String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/testwheatNumNoUn.gff3";
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/AAAtestwheatNumNoUn.gff3";
//            new WheatGeneFeature (infileS ,outfileS);

//            String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatRNA_seq/CountAnalyze/countsRNA_outNo0.txt";
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatRNA_seq/CountAnalyze/countsRNA_ABD.txt";
//            new CountRNAseq (infileS ,outfileS);
//           
//            
                    
//            String infileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/testFasta5000.fa";
//            String infileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/testKGF.txt";
//            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/Z_test/123";
//            new GetWheatpseudoGenome (infileS1,infileS2,outfileS1);     
                    
//            String infileS = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/NumNoUnKGF.gff3";
//            String infileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/SplitSubGenomeA.txt";
//            String infileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/SplitSubGenomeB.txt";
//            String infileS3 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/SplitSubGenomeD.txt";
//            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/BlastFileForMcscanXA.fa";
//            String outfileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/BlastFileForMcscanXB.fa";
//            String outfileS3 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatgenome/BlastFileForMcscanXD.fa";
//            new SplitABD(infileS,infileS1,infileS2,infileS3,outfileS1,outfileS2,outfileS3);
                    
               
            String infileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/filternextAB.txt";
            String infileS2 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/filternextAD.txt";
            String infileS3 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/filternextBD.txt";
            String outfileS1 = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/ABDReciprocalBestD.txt";
            new TheReciprocalBestBLASTHit(infileS1,infileS2,infileS3,outfileS1);
//            
            
//            String infile = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/filterAB.txt";
//            String outfile = "/Users/xuebozhao/Documents/LuLab/WheatEpigenome/wheatHomology/TheReciprocalBestBLAST/filternextAB.txt";
//            new TheReciprocalBestBLASTHit(infile,outfile);



//            String infileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforestMehty_10/siteScoreCpGaddpvalue3X10.txt";
//            String outfileS = "/Users/xuebozhao/Documents/LuLab/CpScore/data/methylation/randomforestMehty_10/AllsiteScoreCpGaddpvalue3X10.txt";
//            new BinomialTestMethylation (infileS,outfileS);

    }  
}
