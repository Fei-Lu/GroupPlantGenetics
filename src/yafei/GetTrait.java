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
 * @create 2020-02-25 10:27 AM
 */
public class GetTrait {
    public static Map getElevation() throws IOException {
        String filename = "/data1/home/yafei/WorldClim/1640_GBS_location.txt";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        Map map = new HashMap();
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("\t"));
            if(list.size() == 19){
                map.put(list.get(0), list.get(18));
            }
        }
        reader.close();
        return map;
    }
    public static Map getBio() throws IOException {
        String filename = "/data1/home/yafei/WorldClim/1640_clim.csv";
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
        File writename = new File("/data1/home/yafei/WorldClim/traitAll.txt");
        BufferedWriter out = new BufferedWriter(new FileWriter(writename));
        Map map = getElevation();
        Map bio = getBio();
        String filename = "/data1/home/yafei/WorldClim/names";
        BufferedReader reader = new BufferedReader(new FileReader(new File(filename)));
        String line;
        out.write("acceNumb\tlon\tlat\tone\ttwo\tthree\tfour\tfive\tsix\tseven\teight\tnine\tten\televen\ttwelve\tthirteen\tfourteen\tfifteen\tsixteen\tseventeen\teighteen\tnineteen\televation\tsrad\n");
        while ((line = reader.readLine()) != null) {
            List<String> list = Arrays.asList(line.split("_"));
           // out.write(line + "\t" + String.valueOf(bio.get(list.get(0)))+"\t" + String.valueOf(map.get(list.get(0)))+"\n");
            out.write(line + "\t" + String.valueOf(bio.get(list.get(0))));
        }
        reader.close();
        out.close();
    }
}
