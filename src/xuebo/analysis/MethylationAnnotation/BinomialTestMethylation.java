/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.MethylationAnnotation;

//import xuebo.analysis.annotation.*;
import cern.jet.random.Binomial;
import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;

/**
 *
 * @author xuebozhao
 */
public class  BinomialTestMethylation{
	//采样100次
	public static void main(String[] args) {
//		for (int i = 0; i < 100; i++) {
                BinomialTest a = new BinomialTest();
                AlternativeHypothesis b = AlternativeHypothesis.TWO_SIDED;
		System.out.println(a.binomialTest(0,10,0.8,b));
//		}
	}
//	//二项分布采样
//	public static double binomialsampler(int trials, double p){
//		BinomialDistribution binomial=new BinomialDistribution(trials,p);
//                System.out.println(p);
//		return binomial.sample();
//	}
      
}
