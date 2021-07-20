package xiaohan.utils;

import xujun.analysis.rnaseq.GeneFeature;

import java.io.*;
import java.net.HttpURLConnection;
import java.net.URL;
import java.net.URLConnection;
import java.util.ArrayList;
import java.util.Arrays;


public class SNPmapping {
    //search from the start. while loop
    public static Integer[] search(int[][] array, int target) {
        if(array==null||array.length<=0||target > array[array.length-1][array[array.length-1].length-1]){
            Integer[] temp = new Integer[]{-1};
            return temp;
        }
        int row = 0;
        int i = 0;
        int result = 0;
        ArrayList<Integer> res = new ArrayList<>();
        while(row <= array.length-1){
            int column = array[row].length-1;
            if (target > array[row][column]){
                row++;i++;
            }
            else {
                result = row;break;
            }
        }
        if (Arrays.asList(array[result]).contains(target)) {
            res.add(result-1);
        }
        res.add(result);
        Integer[] f = res.toArray(new Integer[res.size()]);
        return f;
    }


    //bSearch Without Recursion
    public static int[] bSearch(int[][] a, int key) {
        if(a==null||a.length<=0||key > a[a.length-1][a[a.length-1].length-1]){
            return new int[]{-1};
        }
        int low = 0;
        int high = a.length - 1;
        int mid = 0;
        ArrayList<Integer> res = new ArrayList<>();
        while (low <= high) {
            mid = low + (high - low) / 2;
            if (a[mid][0] > key){
                high = mid -1 ;
            }
            else if (a[mid][1] < key) {
                low = mid +1;
            }
            else {break;}
        }

        if (a[mid][0] > key) { res.add(mid-1);}
        res.add(mid);
        if (a[mid][1] < key) { res.add(mid+1);}

        int[] ints = new int[res.size()];
        for (int i = 0; i < res.size(); i++) { ints[i] = res.get(i); }
        return ints;
//        return new int[]{-1};
    }

    public static int[] bSearch(ArrayList<int[]> a, int key) {
        if(a==null||a.size()<=0||key > a.get(a.size()-1)[1]) {
            return new int[]{-1};
        }
        int low = 0;
        int high = a.size() - 1;
        int mid = 0;
        ArrayList<Integer> res = new ArrayList<>();
        while (low <= high) {
            mid = low + (high - low) / 2;
            if (a.get(mid)[0] > key){
                high = mid -1 ;
            }
            else if (a.get(mid)[1] < key) {
                low = mid +1;
            }
            else {break;}
        }
        if (a.get(mid)[0] > key) { res.add(mid-1);}
        res.add(mid);
        if (a.get(mid)[1] < key) { res.add(mid+1);}

        int[] ints = new int[res.size()];
        for (int i = 0; i < res.size(); i++) { ints[i] = res.get(i); }
        return ints;
//        return new int[]{-1};
    }


    public static String[] mapping(String chr, int pos){
        GeneFeature gf = new GeneFeature("H:/Nature/2020PreExperiment/JNW/Supplement/wheat_v1.1_Lulab.gff3");
        ArrayList<String> subGeneName = new ArrayList<>();
        ArrayList<int[]> subGeneRange = new ArrayList<>();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String temp = gf.getGeneName(i);
            if (temp.contains(chr)){
                subGeneName.add(temp);
                int[] tt = new int[]{gf.getGeneStart(i), gf.getGeneEnd(i)};
                subGeneRange.add(tt);
            }
        }
        int[] temp = bSearch(subGeneRange, pos);
        String[] result = new String[temp.length];
        for (int i = 0; i < temp.length; i++) {
            result[i] = subGeneName.get(temp[i]);
        }
        return result;
    }

    public static String[] mapping(String chr, int pos, int count){
        GeneFeature gf = new GeneFeature("H:/Nature/2020PreExperiment/JNW/Supplement/wheat_v1.1_Lulab.gff3");
        ArrayList<String> subGeneName = new ArrayList<>();
        ArrayList<int[]> subGeneRange = new ArrayList<>();
        for (int i = 0; i < gf.getGeneNumber(); i++) {
            String temp = gf.getGeneName(i);
            if (temp.contains(chr)){
                subGeneName.add(temp);
                int[] tt = new int[]{gf.getGeneStart(i), gf.getGeneEnd(i)};
                subGeneRange.add(tt);
            }
        }
        int[] temp = bSearch(subGeneRange, pos);
        ArrayList<String> tt = new ArrayList<>();

        int start;
        int end;

        if (temp.length == 2){
            if (temp[0] > count -1){
                start = temp[0] - count;
            }else{
                start = 0;
            }
            if (temp[1] < subGeneRange.get(subGeneRange.size() -1)[1] - count){
                end = temp[1] + count;
            }else {
                end = subGeneRange.get(subGeneRange.size())[1];
            }
        }else {
            if (temp[0] < count -1 ){
                start = 0;
                end = temp[0] +count;
            }else if (temp[0] > subGeneRange.get(subGeneRange.size())[1] - count){
                start = temp[0];
                end = subGeneRange.get(subGeneRange.size())[1];
            }else {
                start = temp[0] - count;
                end = temp[0] +count;
            }
        }
        for (int i = start; i <= end; i++) {
            tt.add(subGeneName.get(i));
        }
        String[] result = tt.toArray(new String[tt.size()]);
        return result;
    }

    public static ArrayList<String[]> getTriads(String[] genes) throws IOException{
        BufferedReader br = new BufferedReader(new FileReader("H:/Nature/2020PreExperiment/JNW/Supplement/gene_v1.1_triads.txt"));
        String temp;
        ArrayList<String[]> triads = new ArrayList<>();
        br.readLine();
        ArrayList<String[]> result = new ArrayList<>();

        while ((temp = br.readLine()) != null){
            String[] arr = temp.split("\t");
            triads.add(arr);
        }

        for (String gene:genes) {
            for (int i = 0; i < triads.size(); i++) {
                if (Arrays.asList(triads.get(i)).contains(gene)){
                    result.add(triads.get(i));
                }
            }
        }

        return result;
    }

    public static ArrayList<String[]> flanking(String chr, int pos) throws IOException {
        return getTriads(mapping(chr, pos));
    }

//    public Integer[] judgeIfOut(int target, int [][] array){
//        if(array==null||array.length<=0||target > array[array.length-1][array[array.length-1].length-1]){
//            Integer[] temp = new Integer[]{-1};
//            return temp;
//        }
//    }


//    public static String[] readHomologs(){
//
//    }

    public static String[][] openFile(String filePath) {
        int HttpResult; // 服务器返回的状态
        String[][] ee = new String[][]{};
        try
        {
            URL url =new URL(filePath); // 创建URL
            URLConnection urlconn = url.openConnection(); // 试图连接并取得返回状态码
            urlconn.connect();
            HttpURLConnection httpconn =(HttpURLConnection)urlconn;
            HttpResult = httpconn.getResponseCode();
            if(HttpResult != HttpURLConnection.HTTP_OK) {
                System.out.print("无法连接到");
            } else {
                int filesize = urlconn.getContentLength(); // 取数据长度
                InputStreamReader isReader = new InputStreamReader(urlconn.getInputStream(),"UTF-8");
                BufferedReader reader = new BufferedReader(isReader);
                String str; // 用来保存每行读取的内容
                int i;
                i = 0;
                while ((str = reader.readLine()) != null) {
                    String[] arr = str.split("\t");
                    ee[i] = arr;
                    i++;
                }
            }
        }
        catch (FileNotFoundException e) {
            e.printStackTrace();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return  ee;
    }
}