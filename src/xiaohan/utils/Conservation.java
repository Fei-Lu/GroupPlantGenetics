package xiaohan.utils;

import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Conservation {

    public Conservation(Object arg) {
        this.getGerpValueMap((String) arg);
    }

    public static HashMap<String, String> getGerpValueMap(String infile) {
        HashMap<String, String> SiteGerpMap = new HashMap<>();
        BufferedReader br = IOUtils.getTextGzipReader(infile);
        String temp = null;
        String[] temps = null;
        List<String> tList = new ArrayList();
        try{
            while((temp = br.readLine())!=null){
                tList = PStringUtils.fastSplit(temp);
                temps = tList.toArray(new String[tList.size()]);
//                String site = temps[0] + "_" + temps[2];
                String site = temps[2];
                String Gerp = temps[4];
                SiteGerpMap.put(site,Gerp);
            }
        }catch (Exception e){
            e.printStackTrace();
        }
        return SiteGerpMap;
    }

    public static void main(String[] args) {
        new Conservation(args[0]);
    }
}
