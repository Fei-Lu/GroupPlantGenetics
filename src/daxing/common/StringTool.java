package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.StringUtils;
import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;
import org.jsoup.select.Elements;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * @author Daxing Xu
 */
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
            return Integer.parseInt(m.group());
        }else {
            System.out.println(str+":"+"\t"+"don not contains any number");
            return Integer.MIN_VALUE;
        }
    }

    /**
     * if the given string is numeric (eg. -1.2 0.3 12 et al.), return true
     * @param str
     * @return true, if the given string is all numeric
     */
    public static boolean isNumeric(String str){
        byte[] a=str.getBytes();
        if (a[0]<45 || a[0]>57) return false;
        if (a[0]==46 || a[0]==47) return false;
        if (a[0]==45 && a[1]==46) return false;
        if (a[0]==45 && a[1]==47) return false;
        for (int i = 1; i < a.length; i++) {
            if (a[i]<46 || a[i]>57) return false;
            if (a[i]==47) return false;
        }
        if (a[a.length-1]==46) return false;
        return true;
    }

    /**
     * 返回一个字符串中子字符串出现的所有索引列表
     * @param str 字符串
     * @param subStr 子字符串
     * @return index数组
     */
    public static TIntArrayList getIndexOfSubStr(String str, String subStr){
        TIntArrayList indexes =new TIntArrayList();
        int index=0, wordlen=0;
        while (index!=-1){
            index=str.indexOf(subStr, index+wordlen);
            if (index!=-1){
                indexes.add(index);
            }
            wordlen=subStr.length();
        }
        return indexes;
    }

    /**
     * case-sensitive 大小写敏感
     * @param str
     * @return 子字符串在指定字符串中出现的次数
     */
    public static int getCountOfSubStr(String str, String subStr){
        return StringUtils.countMatches(str, subStr);
    }

    /**
     * 解析HTNL格式的行 如"<td class="numeric" bgcolor="#00fe00"> 2,354 </td>"
     * @param strOfHtml
     * @param tag "body", "td", et al.
     * @return
     */
    public static List<String> parseLineInHtmlFormat(String strOfHtml, String tag){
        Document document= Jsoup.parse(strOfHtml);
        Elements elements=document.getElementsByTag(tag);
        return elements.eachText();
    }

    /**
     * 解析HTNL格式的行 如"<td class="numeric" bgcolor="#00fe00"> 2,354 </td>", 默认的tag为"body"
     * @param strOfHtml
     * @return
     */
    public static List<String> parseLineInHtmlFormat(String strOfHtml){
        return StringTool.parseLineInHtmlFormat(strOfHtml, "body");
    }

    /**
     *
     * @param num  以千位分隔符分隔的字符串数字，如"123,256";
     * @return
     */
    public static int parseInt(String num){
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


}
