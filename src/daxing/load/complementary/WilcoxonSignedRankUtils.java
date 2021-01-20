package daxing.load.complementary;

import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.exception.DimensionMismatchException;
import org.apache.commons.math3.exception.NoDataException;
import org.apache.commons.math3.exception.NullArgumentException;
import org.apache.commons.math3.stat.ranking.NaNStrategy;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.commons.math3.util.FastMath;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * An modified implementation of the Wilcoxon signed-rank test.
 */
public class WilcoxonSignedRankUtils {

    /** Ranking algorithm. */
    private NaturalRanking naturalRanking;

    /**
     * Create a test instance where NaN's are left in place and ties get
     * the average of applicable ranks. Use this unless you are very sure
     * of what you are doing.
     */
    public WilcoxonSignedRankUtils() {
        naturalRanking = new NaturalRanking(NaNStrategy.FIXED,
                TiesStrategy.AVERAGE);
    }

    /**
     * Create a test instance using the given strategies for NaN's and ties.
     * Only use this if you are sure of what you are doing.
     *
     * @param nanStrategy
     *            specifies the strategy that should be used for Double.NaN's
     * @param tiesStrategy
     *            specifies the strategy that should be used for ties
     */
    public WilcoxonSignedRankUtils(final NaNStrategy nanStrategy,
                                  final TiesStrategy tiesStrategy) {
        naturalRanking = new NaturalRanking(nanStrategy, tiesStrategy);
    }

    /**
     * Ensures that the provided arrays fulfills the assumptions.
     *
     * @param x first sample
     * @param y second sample
     * @throws NullArgumentException if {@code x} or {@code y} are {@code null}.
     * @throws NoDataException if {@code x} or {@code y} are zero-length.
     * @throws DimensionMismatchException if {@code x} and {@code y} do not
     * have the same length.
     */
    private void ensureDataConformance(final double[] x, final double[] y)
            throws NullArgumentException, NoDataException, DimensionMismatchException {

        if (x == null ||
                y == null) {
            throw new NullArgumentException();
        }
        if (x.length == 0 ||
                y.length == 0) {
            throw new NoDataException();
        }
        if (y.length != x.length) {
            throw new DimensionMismatchException(y.length, x.length);
        }
    }

    /**
     * Calculates y[i] - x[i] for all i
     *
     * @param x first sample
     * @param y second sample
     * @return z = y - x
     */
    private double[] calculateDifferences(final double[] x, final double[] y) {

        final double[] z = new double[x.length];

        for (int i = 0; i < x.length; ++i) {
            z[i] = y[i] - x[i];
        }

        return z;
    }

    /**
     * Calculates |z[i]| for all i
     *
     * @param z sample
     * @return |z|
     * @throws NullArgumentException if {@code z} is {@code null}
     * @throws NoDataException if {@code z} is zero-length.
     */
    private double[] calculateAbsoluteDifferences(final double[] z)
            throws NullArgumentException, NoDataException {

        if (z == null) {
            throw new NullArgumentException();
        }

        if (z.length == 0) {
            throw new NoDataException();
        }

        final double[] zAbs = new double[z.length];

        for (int i = 0; i < z.length; ++i) {
            zAbs[i] = FastMath.abs(z[i]);
        }

        return zAbs;
    }

    /**
     *
     * @param ranks
     * @return
     */
    private double getVarianceOfTies(double[] ranks){
        Map<Double,Long> res= Arrays.stream(ranks).boxed().parallel().collect(Collectors.groupingBy(i -> i,Collectors.counting()));
        TIntArrayList tiesList=new TIntArrayList();
        for (Map.Entry<Double,Long> entry: res.entrySet()){
            if (entry.getValue() < 2) continue;
            tiesList.add(entry.getValue().intValue());
        }
        double varianceSum=0;
        for (int i = 0; i < tiesList.size(); i++) {
            varianceSum+=(Math.pow(tiesList.get(i), 3)-tiesList.get(i))/48;
        }
        return  varianceSum;
    }

    /**
     * z[i] = y[i]-x[i]
     * @param x pseudo
     * @param y hexaploid
     * @return [W+,normalizedW+]
     */
    public double[] wPlusAndNormalizedWilcoxonSignedRank(double[] x, double[] y){
        ensureDataConformance(x, y);

        double[] res=new double[2];
        Arrays.fill(res, Double.NaN);

        // throws IllegalArgumentException if x and y are not correctly
        // specified
        final double[] z = calculateDifferences(x, y);
        final double[] zAbs = calculateAbsoluteDifferences(z);

        // remove x[i]=y[i] and adjust N
        TDoubleArrayList zList = new TDoubleArrayList();
        TDoubleArrayList zAbsList=new TDoubleArrayList();
        for (int i = 0; i < zAbs.length; i++) {
            if (zAbs[i]==0) continue;
            zList.add(z[i]);
            zAbsList.add(zAbs[i]);
        }

        double n=zList.size();
        if (n*(n+1)/2 < 21) return res;

        final double[] ranks = naturalRanking.rank(zAbsList.toArray());

        double wPlus = 0;

        for (int i = 0; i < zList.size(); ++i) {
            if (zList.get(i) > 0) {
                wPlus += ranks[i];
            }
        }

        double mean=(n*(n+1))/4;
        double variance=(mean*(2*n+1))/6;
        double varianceOfTir=getVarianceOfTies(ranks);
        res[0]=wPlus;
        res[1]=(wPlus-mean)/Math.sqrt(variance-varianceOfTir);
        return res;
    }

    public double getWilcoxonSignedRank(double[] x, double[] y){
        return this.wPlusAndNormalizedWilcoxonSignedRank(x, y)[0];
    }

    public double getNormalizedWilcoxonSignedRank(double[] x, double[] y){
        return this.wPlusAndNormalizedWilcoxonSignedRank(x, y)[1];
    }
}
