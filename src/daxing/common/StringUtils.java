package daxing.common;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class StringUtils {

    public static int getNumFromString(String str){
        Pattern p=Pattern.compile("\\d+");
        Matcher m=p.matcher(str);
        if(m.find()){
            return Integer.valueOf(m.group());
        }else {
            System.out.println(str+":"+"\t"+"don not contains any number");
            System.exit(1);
        }
        return Integer.MIN_VALUE;
    }
}
