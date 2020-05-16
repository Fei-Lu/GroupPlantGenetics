package daxing.load.ancestralSite;

import java.util.Arrays;
import java.util.Comparator;
import java.util.stream.IntStream;

public class Standardization {

    private static final double[] POINTA=initialPoint(0,0);
    private static final double[] POINTB=initialPoint(3, 3*Math.sqrt(3));
    private static final double[] POINTD=initialPoint(6,0);
    private static final double[] POINTE=initialPoint(1.5, Math.sqrt(6.75));
    private static final double[] POINTF=initialPoint(4.5, Math.sqrt(6.75));
    private static final double[] POINTG=initialPoint(3, 0);
    private static final double[] POINTO=initialPoint(3,Math.sqrt(3));

    public enum PointIndex{

//        POINTA(0, "A_dominant"),
//        POINTB(1, "B_dominant"),
//        POINTD(2, "D_dominant"),
//        POINTE(3, "D_suppressed"),
//        POINTF(4, "A_suppressed"),
//        POINTG(5, "B_suppressed"),
//        POINTO(6, "Balanced"),
//        POINTM000(7, "M000");

        POINTA(0, "M100"),
        POINTB(1, "M010"),
        POINTD(2, "M001"),
        POINTE(3, "M110"),
        POINTF(4, "M011"),
        POINTG(5, "M101"),
        POINTO(6, "M111"),
        POINTM000(7, "M000");

        int index;
        String region;

        PointIndex(int index, String region) {
            this.index=index;
            this.region=region;
        }

        public String getRegion(){
            return region;
        }
    }

    private static double[] initialPoint(double x, double y){
        double[] point=new double[2];
        double num=(double)100/6;
        point[0]=x*num;
        point[1]=y*num;
        return point;
    }

    /**
     *
     * @param x1 坐标（a,b）
     * @param x2 坐标 (c,d)
     * @return
     */
    private static double calculateDistance(double[] x1, double[] x2){
        return Math.sqrt(Math.pow(x1[0]-x2[0],2)+Math.pow(x1[1]-x2[1],2));
    }

    /**
     *
     * @param x
     * @return
     */
    public static PointIndex getNearestPointIndex(double[] x){
        if (x[0]==0 && x[1]==0 && x[2] ==0) return PointIndex.POINTM000;
        double[] transform100= transform100(x);
        boolean[] ifZero=new boolean[3];
        for (int i = 0; i < transform100.length; i++) {
            ifZero[i] = transform100[i]==0 ? true : false;
        }
        if (ifZero[1] && ifZero[2]) return PointIndex.POINTA;
        if (ifZero[0] && ifZero[2]) return PointIndex.POINTB;
        if (ifZero[0] && ifZero[1]) return PointIndex.POINTD;
        if (transform100[0]==transform100[1] && transform100[1]==transform100[2]) return PointIndex.POINTO;
        double[] coordinate=getCoordinate(transform100);
        double[] nearestDistance=new double[7];
        nearestDistance[0]=calculateDistance(POINTA, coordinate);
        nearestDistance[1]=calculateDistance(POINTB, coordinate);
        nearestDistance[2]=calculateDistance(POINTD, coordinate);
        nearestDistance[3]=calculateDistance(POINTE, coordinate);
        nearestDistance[4]=calculateDistance(POINTF, coordinate);
        nearestDistance[5]=calculateDistance(POINTG, coordinate);
        nearestDistance[6]=calculateDistance(POINTO, coordinate);
        int minIndex= IntStream.range(0, nearestDistance.length).boxed().min(Comparator.comparingDouble(index->nearestDistance[index])).get();
        return PointIndex.values()[minIndex];
    }

    /**
     *
     * @param x The three elements of the array cannot be 0 at the same time
     * @return
     */
    private static double[] transform100(double[] x){
        if (x[0]==0 && x[1]== 0 && x[2]==0){
            System.out.println("The three elements of the array cannot be 0 at the same time");
            System.exit(1);
        }
        double[] res=new double[x.length];
        for (int i = 0; i < res.length; i++) {
            res[i]=-1;
        }
        double sum= Arrays.stream(x).sum();
        for (int i = 0; i < x.length; i++) {
            res[i]=100*x[i]/sum;
        }
        return res;
    }

    /**
     * 将A B D转换为坐标中的点
     * @param transform100 为transform100后对应的值
     * @return
     */
    private static double[] getCoordinate(double[] transform100){
        double[] coordinate=new double[2];
        coordinate[0]=0.5*transform100[1]+transform100[2];
        coordinate[1]=Math.sqrt(0.75)*transform100[1];
        return coordinate;
    }

    /**
     *
     * @param point1 double[3] 为transform100前对应的值
     * @param point2 double[3] 为transform100前对应的值
     * @return
     */
    public static double getDistance(double[] point1, double[] point2){
        double[] transform100_point1=transform100(point1);
        double[] transform100_point2=transform100(point2);
        double[] pointA=getCoordinate(transform100_point1);
        double[] pointB=getCoordinate(transform100_point2);
        return calculateDistance(pointA, pointB);
    }


}
