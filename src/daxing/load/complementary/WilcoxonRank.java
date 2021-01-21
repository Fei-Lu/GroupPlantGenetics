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
     *
     * @param ranks rank result
     * @return variance of ties
     */
    private static double getVarianceOfTies(double[] ranks){
        TIntArrayList tiesList=WilcoxonRank.getTies(ranks);
        double varianceSum=0;
        for (int i = 0; i < tiesList.size(); i++) {
            varianceSum+=(Math.pow(tiesList.get(i), 3)-tiesList.get(i))/48;
        }
        return  varianceSum;
    }

    public static TIntArrayList getTies(double[] ranks){
        Map<Double,Long> res= Arrays.stream(ranks).boxed().parallel().collect(Collectors.groupingBy(i -> i,Collectors.counting()));
        TIntArrayList tiesList=new TIntArrayList();
        for (Map.Entry<Double,Long> entry: res.entrySet()){
            if (entry.getValue() < 2) continue;
            tiesList.add(entry.getValue().intValue());
        }
        return tiesList;
    }

    /**
     * z[i] = y[i]-x[i]
     * @param x pseudo
     * @param y hexaploid
     * @return [W+,normalizedW+]
     */
    public static double[] wPlusAndNormalizedWilcoxonSignedRank(double[] x, double[] y){
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
        double varianceOfTir=WilcoxonRank.getVarianceOfTies(ranks);
        res[0]=wPlus;
        res[1]=(wPlus-mean)/Math.sqrt(variance-varianceOfTir);
        return res;
    }

    public static double getWilcoxonSignedRank(double[] x, double[] y){
        return WilcoxonRank.wPlusAndNormalizedWilcoxonSignedRank(x, y)[0];
    }

    public static double getNormalizedWilcoxonSignedRank(double[] x, double[] y){
        return WilcoxonRank.wPlusAndNormalizedWilcoxonSignedRank(x, y)[1];
    }
}
