package daxing.load.complementary;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.util.FastMath;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * An modified implementation of the Wilcoxon signed-rank test.
 * It can perform Wilcoxon signed-rank test and Wilcxon rank-sum test(also called Mann–Whitney–Wilcoxon) when sample
 * size greater or equal than 20
 */
public class WilcoxonRank {

    /**
     * Calculates y[i] - x[i] for all i
     * @param x first sample
     * @param y second sample
     * @return z = y - x
     */
    private static double[] calculateDifferences(double[] x,double[] y) {
        double[] z = new double[x.length];
        for (int i = 0; i < x.length; ++i) {
            z[i] = y[i] - x[i];
        }
        return z;
    }

    /**
     * Calculates |z[i]| for all i
     * @param z sample
     * @return |z|
     */
    private static double[] calculateAbsoluteDifferences(double[] z){
        double[] zAbs = new double[z.length];
        for (int i = 0; i < z.length; ++i) {
            zAbs[i] = FastMath.abs(z[i]);
        }
        return zAbs;
    }

    /**
     * calculate ties variance of wilcoxon signed-rank
     * @param ranks rank result
     * @return variance of ties
     */
    private static double getSignedRankVarianceOfTies(double[] ranks){
        TIntArrayList tiesList=WilcoxonRank.getTies(ranks);
        double varianceSum=0;
        for (int i = 0; i < tiesList.size(); i++) {
            varianceSum+=(Math.pow(tiesList.get(i), 3)-tiesList.get(i))/48;
        }
        return  varianceSum;
    }

    /**
     * calculate ties variance of wilcoxon rank-sum
     * @param ranks rank result
     * @return variance of ties
     */
    private static double getRankSumVarianceOfTies(double[] ranks){
        TIntArrayList tiesList=WilcoxonRank.getTies(ranks);
        double varianceSum=0;
        double n=ranks.length;
        for (int i = 0; i < tiesList.size(); i++) {
            varianceSum+=(Math.pow(tiesList.get(i), 3)-tiesList.get(i))/(n*(n-1));
        }
        return varianceSum;
    }

    /**
     *
     * @param ranks rank result
     * @return the number of occurrences (>1) of each element in the ranks
     */
    private static TIntArrayList getTies(double[] ranks){
        Map<Double,Long> res= Arrays.stream(ranks).boxed().parallel().collect(Collectors.groupingBy(i -> i,Collectors.counting()));
        TIntArrayList tiesList=new TIntArrayList();
        for (Map.Entry<Double,Long> entry: res.entrySet()){
            if (entry.getValue() < 2) continue;
            tiesList.add(entry.getValue().intValue());
        }
        return tiesList;
    }

    /**
     * z[i] = y[i]-x[i], x.length must be equal y.length
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return [W+,normalizedW+]
     */
    private static double[] wPlusAndNormalizedWilcoxonSignedRank(double[] x, double[] y){
        if (x.length != y.length){
            System.out.println("WilcoxonSignedRank: x.length must be equal y.length");
            System.exit(1);
        }
        double[] res=new double[2];
        Arrays.fill(res, Double.NaN);
        double[] z = WilcoxonRank.calculateDifferences(x, y);
        double[] zAbs = WilcoxonRank.calculateAbsoluteDifferences(z);
        // remove x[i]=y[i] and adjust N
        TDoubleArrayList zList = new TDoubleArrayList();
        TDoubleArrayList zAbsList=new TDoubleArrayList();
        for (int i = 0; i < zAbs.length; i++) {
            if (zAbs[i]==0) continue;
            zList.add(z[i]);
            zAbsList.add(zAbs[i]);
        }
        double n=zList.size();
        // if n less than 20, exact distribution needs to be used
        // return NaN
        if (n < 20) return res;
        NaturalRanking naturalRanking=new NaturalRanking();
        double[] ranks = naturalRanking.rank(zAbsList.toArray());
        double wPlus = 0;
        for (int i = 0; i < zList.size(); ++i) {
            if (zList.get(i) > 0) {
                wPlus += ranks[i];
            }
        }
        double mean=(n*(n+1))/4;
        double variance=(mean*(2*n+1))/6;
        double varianceOfTir=WilcoxonRank.getSignedRankVarianceOfTies(ranks);
        res[0]=wPlus;
        // - 0.5 is a continuity correction
        res[1]=(wPlus-mean-0.5)/Math.sqrt(variance-varianceOfTir);
        return res;
    }

    /**
     *
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return WilcoxonSignedRank statics
     */
    public static double getWilcoxonSignedRank(double[] x, double[] y){
        return WilcoxonRank.wPlusAndNormalizedWilcoxonSignedRank(x, y)[0];
    }

    /**
     *
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return Normalized WilcoxonSignedRank statics
     */
    public static double getNormalizedWilcoxonSignedRank(double[] x, double[] y){
        return WilcoxonRank.wPlusAndNormalizedWilcoxonSignedRank(x, y)[1];
    }

    /**
     *
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return [Uy,normalizedUy]
     */
    private static double[] u1AndNormalizedWilcoxonRankSum(double[] x, double[] y){
        double[] res=new double[2];
        Arrays.fill(res, Double.NaN);
        // if n less than 20, exact distribution needs to be used
        // return NaN
        if ((x.length+y.length) < 20) return res;
        double[] mergedSample=new double[x.length+y.length];
        System.arraycopy(x, 0, mergedSample, 0, x.length);
        System.arraycopy(y, 0, mergedSample, x.length, y.length);
        NaturalRanking naturalRanking=new NaturalRanking();
        double[] mergedRank=naturalRanking.rank(mergedSample);
        double n1=x.length;
        double n2=y.length;
        double r2=0;
        for (int i = x.length; i < mergedRank.length; i++) {
            r2+=mergedRank[i];
        }
        double u2=r2-(n2*(n2+1))/2;
        double mean=n1*n2/2;
        double varianceOfTie=WilcoxonRank.getRankSumVarianceOfTies(mergedRank);
        double sd=Math.sqrt((mean/6)*((n1+n2+1)-varianceOfTie));
        res[0]=u2;
        // - 0.5 is a continuity correction
        res[1]=(u2-mean-0.5)/sd;
        return res;
    }

    /**
     *
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return WilcoxonRankSum statics
     */
    public static double getWilcoxonRankSum(double[] x, double[] y){
        return WilcoxonRank.u1AndNormalizedWilcoxonRankSum(x, y)[0];
    }

    /**
     *
     * @param x sample x (pseudo)
     * @param y sample y (hexaploid)
     * @return Normalized WilcoxonRankSum statics
     */
    public static double getNormalizedWilcoxonRankSum(double[] x, double[] y){
        return WilcoxonRank.u1AndNormalizedWilcoxonRankSum(x, y)[1];
    }
}
