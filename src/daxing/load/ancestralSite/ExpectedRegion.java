package daxing.load.ancestralSite;

public class ExpectedRegion {

    public enum Expected{
        H, L, N
    }

    /**
     *
     * @param p_valueOfABD 六个p_value A B D各两个
     * @return
     */
    public static String getExpectedRegion(double[] p_valueOfABD){
        StringBuilder sb=new StringBuilder();
        sb.setLength(0);
        for (int i = 0; i < p_valueOfABD.length; i=i+2) {
            if (p_valueOfABD[i]<0.05){
                sb.append(Expected.L);
            }else if (p_valueOfABD[i+1]<0.05){
                sb.append(Expected.H);
            }else {
                sb.append(Expected.N);
            }
        }
        return sb.toString();
    }
}
