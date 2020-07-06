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
 * @create 2020-06-12 8:33 PM
 */
public class GetLoc {
    public static Map getLoc() throws IOException {
        String filename = "/Users/guoyafei/Desktop/1640_clim.csv";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        Map map = new HashMap();
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split(","));
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
        File writename = new File("/Users/guoyafei/Desktop/traitall.txt");
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        Map bio = getLoc();
        String filename = "/Users/guoyafei/Desktop/names";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line = null;
        String line2;
        out.write("acceNumb\tlon\tlat\n");
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("_"));
            // out.write(line + "\t" + String.valueOf(bio.get(list.get(0)))+"\t" + String.valueOf(map.get(list.get(0)))+"\n");
            line2 = String.valueOf(bio.get(list.get(0)));
            //System.out.println(line + "\t" + String.valueOf(bio.get(list.get(0))));
            if(line2 != null){
                out.write(line + "\t" + String.valueOf(bio.get(list.get(0)))+ "\n");
            }
        }
        reader.close();
        out.close();
    }
}

