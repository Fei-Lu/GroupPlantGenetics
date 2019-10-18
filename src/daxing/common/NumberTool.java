package daxing.common;

import java.math.BigDecimal;
import java.text.DecimalFormat;
import java.text.NumberFormat;
import java.text.ParseException;

public class NumberTool {

    /**
     * 四舍五入
     * @param num
     * @param digit 小数点后保留位数
     * @return
     */
    public static double format(double num, int digit){
        BigDecimal d=new BigDecimal(num);
        return d.setScale(digit, BigDecimal.ROUND_HALF_UP).doubleValue();
    }

    /**
     *
     * @param num  以千位分隔符分隔的字符串数字，如"123,256";
     * @return
     */
    public static int intValue(String num){
        int res=Integer.MIN_VALUE;
        try {
            NumberFormat format=NumberFormat.getInstance();
            Number number=format.parse(num);
            res=number.intValue();
        } catch (ParseException e) {
            e.printStackTrace();
        }
        return res;
    }

    /**
     *
     * @param num  以千位分隔符分隔的字符串数字，如"123,453,234.12";
     * @return
     */
    public static double doubleValue(String num){
        double res=Integer.MIN_VALUE;
        try {
            NumberFormat format=NumberFormat.getInstance();
            Number number=format.parse(num);
            res=number.doubleValue();
        } catch (ParseException e) {
            e.printStackTrace();
        }
        return res;
    }

    /**
     *
     * @param a
     * @return
     */
    public static String parse(int a){
        DecimalFormat df = new DecimalFormat();
        return df.format(a);
    }

    /**
     *
     * @param a
     * @return
     */
    public static String parse(double a){
        DecimalFormat df = new DecimalFormat();
        return df.format(a);
    }


}
