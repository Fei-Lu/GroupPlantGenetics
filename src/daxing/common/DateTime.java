package daxing.common;

import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;

public class DateTime {
/**
 * @author Daxing Xu
 */

    /**
     * 2019-10-26 18:03:08
     * @return the date and time of the current system
     */
    public static String getDateTimeOfNow(){
        LocalDateTime localDateTime=LocalDateTime.now();
        return localDateTime.format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss"));
    }

}
