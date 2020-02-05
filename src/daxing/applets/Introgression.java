package daxing.applets;

import gnu.trove.list.array.TDoubleArrayList;
import org.apache.commons.math3.util.MathUtils;

import java.util.Arrays;

public class Introgression {

    public static double caculate_D(double[] p1DerivedArray, double[] p2DerivedArray, double[] p3DerivedArray ){
        int len=p1DerivedArray.length;
        if (len!=p2DerivedArray.length || len!= p3DerivedArray.length){
            System.out.println("check parameter, program quit");
            System.exit(1);
        }
        double[] abba=new double[len];
        double[] baba=new double[len];
        for (int i = 0; i < len; i++) {
            abba[i]=(1-p1DerivedArray[i])*(p2DerivedArray[i])*(p3DerivedArray[i]);
            baba[i]=(p1DerivedArray[i])*(1-p2DerivedArray[i])*(p3DerivedArray[i]);
        }
        double sumABBA= Arrays.stream(abba).sum();
        double sumBABA=Arrays.stream(baba).sum();
        return (sumABBA-sumBABA)/(sumABBA+sumBABA);
    }

    public static double caculate_fd(double[] p1DerivedArray, double[] p2DerivedArray, double[] p3aDerivedArray,
                                     double[] p3bDerivedArray){
        int len=p1DerivedArray.length;
        if (len!=p2DerivedArray.length || len != p3aDerivedArray.length || len!=p3bDerivedArray.length){
            System.out.println("check parameter, program quit");
            System.exit(1);
        }
        double[] abba_numerator=new double[len];
        double[] baba_numerator=new double[len];
        double[] abba_denominator=new double[len];
        double[] baba_denominator=new double[len];
        for (int i = 0; i < len; i++) {
            abba_numerator[i]=(1-p1DerivedArray[i])*p2DerivedArray[i]*p3aDerivedArray[i];
            baba_numerator[i]=p1DerivedArray[i]*(1-p2DerivedArray[i])*p3aDerivedArray[i];
            abba_denominator[i]=(1-p1DerivedArray[i])*p3bDerivedArray[i]*p3aDerivedArray[i];
            baba_denominator[i]=p1DerivedArray[i]*(1-p3bDerivedArray[i])*p3aDerivedArray[i];
        }
        double sum_abba_numerator= Arrays.stream(abba_numerator).sum();
        double sum_baba_numerator=Arrays.stream(baba_numerator).sum();
        double sum_abba_denominator=Arrays.stream(abba_denominator).sum();
        double sum_baba_denominator=Arrays.stream(baba_denominator).sum();
        return (sum_abba_numerator-sum_baba_numerator)/(sum_abba_denominator-sum_baba_denominator);
    }


}
