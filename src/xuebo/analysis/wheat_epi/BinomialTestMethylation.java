/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.wheat_epi;

//import xuebo.analysis.annotation.*;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

/**
 *
 * @author xuebozhao
 */
public class  BinomialTestMethylation{
	//采样100次
//	public static void main(String[] args) {
////		for (int i = 0; i < 100; i++) {
//                BinomialTest a = new BinomialTest();
//                AlternativeHypothesis b = AlternativeHypothesis.TWO_SIDED;
//		System.out.println(a.binomialTest(0,10,0.8,b));
////		}
//	}
//////	//二项分布采样
//////	public static double binomialsampler(int trials, double p){
//////		BinomialDistribution binomial=new BinomialDistribution(trials,p);
//////                System.out.println(p);
//////		return binomial.sample();
//////	}
    public BinomialTestMethylation(String infileS, String outfileS1, String outfileS2,String outfileS3) {

        this.getMethylationiBnomialAnnotation(infileS,outfileS1, outfileS2,outfileS3);
    }
        
        
    public void getMethylationiBnomialAnnotation(String infileS, String outfileS1, String outfileS2,String outfileS3) {

        try {
            
            BufferedReader br;
            
            if (infileS.endsWith("gz")) {

                br = XueboIOUtils.getTextGzipReader(infileS);

            } else {

                br = XueboIOUtils.getTextReader(infileS);
            }
            
            String temp = null;
            int i = 0;
            double methylationScore = 0;   
            double pvalue = 0;
            String LCpG = null;
            String LCHH = null;
            String LCHG = null;     
            int Pcategorical = 0;
            BufferedWriter bw1 = XueboIOUtils.getTextWriter(outfileS1);
            BufferedWriter bw2 = XueboIOUtils.getTextWriter(outfileS2);
            BufferedWriter bw3 = XueboIOUtils.getTextWriter(outfileS3);            
            BinomialTest a = new BinomialTest();
            AlternativeHypothesis b = AlternativeHypothesis.TWO_SIDED;
            
            while (( temp = br.readLine()) != null) {
                
                    ++i;
                    
                    if (i % 10000000 == 0) {

                    System.out.println("MethylationAnnotation" + i + "....");

                    }
                    
                    String[] tem = temp.split("\t"); 
                    
                    boolean outWrite = false;
                    
                    int methylationC = Integer.valueOf(tem[5]);                   
                    int methylationT = Integer.valueOf(tem[5])+ Integer.valueOf(tem[6]);
                    methylationScore = methylationC / methylationT;                    
                    pvalue = a.binomialTest(methylationT,methylationC,0.3075,b);
                    
                    if( pvalue < 0.001 ){
                        Pcategorical = 0;
                    }
                    else{
                        Pcategorical = 1;
                    }
                    
                    if(tem[3].equals("CG")){
                        
                        LCpG = tem[0] + "\t" + tem[1] + "\t"  + methylationScore  + "\t"  + Pcategorical + "\t"  + pvalue ;
                        bw1.write(LCpG + "\n"); 
                    }
                    
                    if(tem[3].equals("CHH")){
                        
                        LCHH = tem[0] + "\t" + tem[1] + "\t"  + methylationScore + "\t"  + Pcategorical + "\t"  + pvalue ;
                        bw2.write(LCHH + "\n");
                    }
                    
                    if(tem[3].equals("CHG")){
                        
                        LCHG = tem[0] + "\t" + tem[1] + "\t"  + methylationScore  + "\t"  + Pcategorical + "\t"  + pvalue ;
                        bw3.write(LCHG + "\n");
                    }
                    
//                    if (outWrite) {
//
//                    bw1.write(LCpG + "\n");
//                    bw2.write(LCHH + "\n");
//                    bw3.write(LCHG + "\n");
//                    
//                    bw1.flush();
//                    bw2.flush();
//                    bw3.flush();
//                    }
            }
            bw1.close();
            bw2.close();
            bw3.close();
        }catch (Exception e) {
            e.printStackTrace();
        }
    
    }
}
