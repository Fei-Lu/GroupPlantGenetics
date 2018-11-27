/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import utils.IOFileFormat;
import utils.IOUtils;
import utils.PStringUtils;

/**
 *
 * @author Aoyue
 */
public class BWAJiao {
    public BWAJiao(){
        //this.getName();
        //this.sort();
        //this.checkDupli();
        this.newCheckSecondData();
        
    }
    
    /**
     * 查找第二批数据是否和上次所得的 Jiao_nameList_NOdata.txt 信息一致，如果一致，说明数据正确无误。
     * 把 dbfileS 当成库，转化为数组，然后将infileS读进去，一一搜索。
     */
    public void newCheckSecondData(){
        String infileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60taxaname_later.txt";
        String dbfileS = "/Users/Aoyue/Documents/Jiao/001_db/Jiao_nameList_NOdata.txt";
        
        RowTable<String> t = new RowTable<>(dbfileS);
        List<String> l = t.getColumn(0);
        String[] db = l.toArray(new String[l.size()]);
        Arrays.sort(db);  /*千万不能忘了排序啊！！！*/
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            String header = br.readLine();
            String temp = null;
            int cnt =0;
            while((temp = br.readLine()) != null){
                String query = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(db, query);
                if(index <0){
                    System.out.println(query);
                    cnt++;  
                }  
            }
            System.out.println(cnt);
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
        /**
         * 最后结果，一一对应，说明数据没有问题。
         */
    }
    
    public void checkDupli(){
        String infileS = "/Users/Aoyue/Documents/Jiao120taxa_country.txt";
        String Jiao60fileS = "/Users/Aoyue/Documents/Jiao60taxaname.txt";
        String outfileS = "/Users/Aoyue/Documents/Jiao_nameList_NOdata.txt";
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsText(0);
        List<String> l = t.getColumn(0);
        Set s = new HashSet(l);
        System.out.println(l.size());
        System.out.println(s.size());
        
        t = new RowTable<String>(Jiao60fileS);
        l = t.getColumn(0);
        String[] taxa60 = l.toArray(new String[l.size()]);
        Arrays.sort(taxa60);
        
        try{
            BufferedReader br = IOUtils.getTextReader(infileS);
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            String header = br.readLine();
            bw.write(header);bw.newLine();
            String temp = null;
            while((temp = br.readLine()) != null){
                String query = PStringUtils.fastSplit(temp).get(0);
                int index = Arrays.binarySearch(taxa60, query);
                if (index < 0) {
                    bw.write(temp);bw.newLine();
                }
            }
            bw.flush();bw.close();br.close();
        }
        catch(Exception e){
            //System.out.println(temp);
            e.printStackTrace();
            System.exit(1);
            
        }
        
        t = new RowTable<String>(outfileS);
        List<String> country = t.getColumn(1);
        Set<String> sco = new HashSet<>(country);
        for(String a : sco){
            System.out.println(a + "    " + Collections.frequency(country, a));
        }

    }
    
    public void sort(){
//        String infileS = "/Users/Aoyue/Documents/Jiao60name.txt";
//        String outfileS = "/Users/Aoyue/Documents/Jiao60taxaname.txt";
        
        String infileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60_later.txt";
        String outfileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60taxaname_later.txt";
        
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsText(0);
        t.writeTextTable(outfileS, IOFileFormat.Text);
        
        
    }
    
    public void getName(){
//        String infileDirS = "/Volumes/WheatAABBDD/ABD001";
//        String outfileS = "/Users/Aoyue/Documents/Jiao1.txt";
        
//        String infileDirS = "/Volumes/ABD002/ABD002";
//        String outfileS = "/Users/Aoyue/Documents/Jiao2.txt";
        
//        String infileDirS = "/Volumes/AABBDD003/ABD003";
//        String outfileS = "/Users/Aoyue/Documents/Jiao3.txt";
        
        String infileDirS = "/Volumes/Seagate Backup Plus Drive/NHT151096_60s/release_1";
        String outfileS = "/Users/Aoyue/Documents/Jiao/001_db/60sample_later/Jiao60_later.txt";
        
        
        Set<String> s = new HashSet<>();
        File[] fs = new File(infileDirS).listFiles();
        fs = IOUtils.listFilesEndsWith(fs, ".gz");
        int cnt = 0;
        for(int i=0; i< fs.length;i++){
            String name = fs[i].getName().split("_")[0];
            System.out.println(name);
            s.add(name);
            cnt++;
        }
        System.out.println(cnt + "\tfq.gz 文件数");
        System.out.println(s.size() + "\tname文件数");
        
        String[] names = s.toArray(new String[s.size()]);
        Arrays.sort(names);
        
        try{
            BufferedWriter bw = IOUtils.getTextWriter(outfileS);
            bw.write("Name");bw.newLine();
            for(int i=0; i< names.length;i++){
                bw.write(names[i]);bw.newLine();
            }
            bw.flush();bw.close();
            
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
            
        }
    }
    
    
    public static void main (String[] args){
        new BWAJiao();
    }  
}
