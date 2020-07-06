/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package yafei;

/**
 *
 * @author guoyafei
 */

import java.io.*;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

/**
 * @author Yafei Guo
 * @create 2020-03-05 10:24 AM
 */
public class ExtractInfo {
    public static Map getInfo() throws IOException {
        String filename = "/Users/guoyafei/Desktop/TREE/NEW/admixture";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        Map map = new HashMap();
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            StringBuilder sb = new StringBuilder();
            for (int i = 1; i < list.size(); i++) {
                sb.append(list.get(i)+"\t");
            }
            map.put(list.get(0),sb.toString().trim());
        }
        reader.close();
        return map;
    }
    public static void main(String[] args) throws IOException {
        File writename = new File("/Users/guoyafei/Desktop/TREE/NEW/labels");
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        Map info = getInfo();
        String filename = "/Users/guoyafei/Desktop/TREE/NEW/names";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        out.write("acceNumb\tcontinents\tsub_continents\n");
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("_"));
            if(info.keySet().contains(list.get(0))){
                out.write(String.valueOf(info.get(list.get(0)))+"\n");
            }
        }
        reader.close();
        out.close();
    }
}

