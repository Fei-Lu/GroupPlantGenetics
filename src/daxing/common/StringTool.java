package daxing.common;

import java.util.ArrayList;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class StringTool {

    /**
     * 返回一个字符串中首次出现的数字,如果没有则返回Integer最小值
     * @param str 字符串
     * @return 返回一个字符串中首次出现的数字
     */
    public static int getNumFromString(String str){
        Pattern p=Pattern.compile("\\d+");
        Matcher m=p.matcher(str);
        if(m.find()){
            return Integer.valueOf(m.group());
        }else {
            System.out.println(str+":"+"\t"+"don not contains any number");
            return Integer.MIN_VALUE;
        }
    }

    /**
     * 返回一个字符串中子字符串出现的index数组
     * @param str 字符串
     * @param subStr 子字符串
     * @return index数组
     */
    public static int[] getFrequencyOfSubStr(String str, String subStr){
        List<Integer> indexes =new ArrayList<>();
        int index=0, wordlen=0;
        while (index!=-1){
            index=str.indexOf(subStr, index+wordlen);
            if (index!=-1){
                indexes.add(index);
            }
            wordlen=subStr.length();
        }
        return indexes.stream().mapToInt(Integer::intValue).toArray();
    }
}
