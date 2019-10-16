package daxing.common;

import java.math.BigDecimal;

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
}
