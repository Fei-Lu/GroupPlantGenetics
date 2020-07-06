package daxing.common;

/**
 * @author Daxing Xu
 * coordinate transformation between maf file (0-based, half-open, [a, b)) and reference genome (1-based,
 * fully-closed, [a, b])
 */
public class CoordinateTransform {

    /**
     *
     * @param oneBasedFullyClosedRange 1-start, fully-closed [a, b]
     * @return zeroBasedHalfOpenRange 0-start, half-open   [a, b)
     */
    public static int[] getZeroBasedHalfOpenRange(int[] oneBasedFullyClosedRange){
        if (oneBasedFullyClosedRange.length != 2){
            System.out.println("error, the length of oneBasedFullyClosedRange must be 2!!!");
            System.exit(1);
        }
        int[] zeroBasedHalfOpen=new int[2];
        zeroBasedHalfOpen[0]=oneBasedFullyClosedRange[0]-1;
        zeroBasedHalfOpen[1]=oneBasedFullyClosedRange[1];
        return zeroBasedHalfOpen;
    }

    /**
     *
     * @param zeroBasedHalfOpenRange 0-start, half-open   [a, b)
     * @return oneBasedFullyClosedRange 1-start, fully-closed” [a, b]
     */
    public static int[] getOneBasedFullyClosedRange(int[] zeroBasedHalfOpenRange){
        if (zeroBasedHalfOpenRange.length != 2){
            System.out.println("error, the length of zeroBasedHalfOpenRange must be 2!!!");
            System.exit(1);
        }
        int[] oneBasedFullyClosedRange=new int[2];
        oneBasedFullyClosedRange[0]=zeroBasedHalfOpenRange[0]+1;
        oneBasedFullyClosedRange[1]=zeroBasedHalfOpenRange[1];
        return oneBasedFullyClosedRange;
    }

    /**
     *
     * @param oneBasedFullyClosedRange 1-start, fully-closed [a, b]
     * @param chrSize
     * @return reverseStrandZeroBasedHalfOpenRange 0-start, half-open [revStart, revEnd)
     */
    public static int[] getReverseStrandZeroBasedHalfOpenRange(int[] oneBasedFullyClosedRange, int chrSize){
        if (oneBasedFullyClosedRange.length != 2){
            System.out.println("error, the length of oneBasedFullyClosedRange must be 2!!!");
            System.exit(1);
        }
        int revStart = chrSize - oneBasedFullyClosedRange[1];
        int revEnd = chrSize - (oneBasedFullyClosedRange[0] - 1);
        int[] reverseStrandZeroBasedHalfOpenRange=new int[2];
        reverseStrandZeroBasedHalfOpenRange[0]=revStart;
        reverseStrandZeroBasedHalfOpenRange[1]=revEnd;
        return reverseStrandZeroBasedHalfOpenRange;
    }

    /**
     *
     * @param reverseStrandZeroBasedHalfOpenRange reverse strand 0-start, half-open [revStart, revEnd)
     * @param chrSize
     * @return forward strand 1-start, fully-closed [a, b]
     */
    public static int[] getForwardStrandOneBasedFullyClosedRange(int[] reverseStrandZeroBasedHalfOpenRange,
                                                                 int chrSize){
        if (reverseStrandZeroBasedHalfOpenRange.length != 2){
            System.out.println("error, the length of reverseStrandZeroBasedHalfOpenRange must be 2!!!");
            System.exit(1);
        }
        int oneStart = chrSize - reverseStrandZeroBasedHalfOpenRange[1] + 1;
        int oneEnd = chrSize - reverseStrandZeroBasedHalfOpenRange[0];
        int[] oneBasedFullyClosedRange=new int[2];
        oneBasedFullyClosedRange[0]=oneStart;
        oneBasedFullyClosedRange[1]=oneEnd;
        return oneBasedFullyClosedRange;
    }


}
