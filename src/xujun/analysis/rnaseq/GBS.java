/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

import pgl.infra.range.RangeValStr;
import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import pgl.infra.utils.IOUtils;
import xuebo.analysis.annotation.FStringUtils;

/**
 *
 * @author xujun
 */
public class GBS {
    public GBS(){
//        this.barcode();
        this.DisIndex();
//        this.barcode();
    }
    private void DisIndex(){
        String sampleInfor="/Users/xujun/Desktop/indexNews.txt";
        RowTable rt=new RowTable(sampleInfor); 
        HashMap<String, List> indexBarcodeListMap = new HashMap<>();
        HashMap barcodeSampleMap=new HashMap();
        List library=new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            barcodeSampleMap.put( rt.getCellAsString(i, 3),rt.getCellAsString(i, 4));
            List barcodePool = new ArrayList();
//            indexBarcodeListMap.put(rt.getCellAsString(i, 2), rt.getCell(i, 3));
            if(!library.contains(rt.getCell(i, 0))){
                library.add(rt.getCell(i, 0));
            }
        }
        String infileDirS = "/Users/xujun/Desktop/20180705mix1_R1.fq.gz";
        String outfileDirS = "/Users/xujun/Desktop/DisLibrary.txt";
        int count [] =new int [library.size()];
            try {
                BufferedReader br = IOUtils.getTextGzipReader(infileDirS);
                BufferedWriter bw = IOUtils.getTextWriter(outfileDirS);
                String temp1=null;String index=null;
                while((temp1 = br.readLine()) != null){
                    index=temp1.split(" ")[1].split(":")[3].substring(1, 6);
//                    if(hm.get(index)!=null){
//                        int pos=library.indexOf(hm.get(index));
//                        count[pos]++;
//                    }
                    br.readLine();br.readLine();br.readLine();
                    
                }
                for(int i=0;i< count.length;i++){
                    bw.write(library.get(i)+"\t"+count[i]);
                    bw.newLine();
                }
                br.close();
                bw.flush();bw.close();
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

    }
    private void DisIndex1(){
        String sampleInfor="/Users/xujun/Desktop/indexNews.txt";
        RowTable rt=new RowTable(sampleInfor);
        HashMap hm=new HashMap();
        HashMap<String, BufferedWriter> libraryWriterMap1 = new HashMap<>(); 
        HashMap<String, BufferedWriter> libraryWriterMap2 = new HashMap<>();
        List library=new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            hm.put( rt.getCellAsString(i, 2).substring(1),rt.getCell(i, 0));
            if(!library.contains(rt.getCell(i, 0))){
                library.add(rt.getCell(i, 0));                
            }
        }
        BufferedWriter[] bws1 = new BufferedWriter[library.size()];
        BufferedWriter[] bws2 = new BufferedWriter[library.size()];
        for(int i=0;i<library.size();i++){
            bws1[i] = IOUtils.getTextWriter(library.get(i).toString()+"R1");
            libraryWriterMap1.put(library.get(i).toString(), bws1[i]);
            bws2[i] = IOUtils.getTextWriter(library.get(i).toString()+"R2");
            libraryWriterMap2.put(library.get(i).toString(), bws2[i]);
        }
        
        String infileDirS = "/Users/xujun/Desktop";
        String outfileDirS = "/Users/xujun/Desktop";
        int count [] =new int [library.size()];
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"R1.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"R2.fq.gz").getAbsolutePath();
            String outfile = new File (outfileDirS).getAbsolutePath();
            try {            
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);
                BufferedWriter bw1 = IOUtils.getTextWriter(outfileDirS);
                BufferedWriter bw2 = IOUtils.getTextWriter(outfileDirS);
                BufferedWriter tw1 = null;BufferedWriter tw2 = null;
                String temp1=null;String temp2=null;String index=null;
                while((temp1 = br1.readLine()) != null){
                    temp2=br2.readLine();
                    index=temp1.split(" ")[1].split(":")[3].substring(1, 6);
                    if(hm.get(index)!=null){
                        tw1 = libraryWriterMap1.get(hm.get(index));
                        tw1.write(temp1);tw1.newLine();
                        tw1.write(br1.readLine());tw1.newLine();tw1.write(br1.readLine());tw1.newLine();tw1.write(br1.readLine());tw1.newLine();   
                        tw2 =libraryWriterMap2.get(hm.get(index));
                        tw2.write(temp1);tw2.newLine();
                        tw2.write(br2.readLine());tw2.newLine();tw2.write(br2.readLine());tw2.newLine();tw2.write(br2.readLine());tw2.newLine();   
                    }else{
                        br1.readLine();br1.readLine();br1.readLine();
                        br2.readLine();br2.readLine();br2.readLine();
                    }
                    
                    
                }
                br1.close();br2.close();
                for (int i = 0; i < library.size(); i++) {
                    bws1[i].flush();bws1[i].close();
                    bws2[i].flush();bws2[i].close();
                }
            }
            
            catch (Exception e) {
                e.printStackTrace();
            }

        });
            
    }
    private void barcode () {
        String sampleInfor="/Users/xujun/Desktop/barcodeNews.txt";
        RowTable rt=new RowTable(sampleInfor);
        HashMap hm1=new HashMap();
        HashMap hm2=new HashMap();
        List sample =new ArrayList();
        for(int i=0;i<rt.getRowNumber();i++){
            hm1.put(rt.getCell(i, 1), rt.getCellAsString(i, 3).substring(1));
            hm2.put(rt.getCell(i, 1), rt.getCellAsString(i, 4).substring(1));
            sample.add(rt.getCell(i, 1));
        }
        String infileDirS = "/Users/xujun/Desktop/s20180622-LHY2_TKD180602248";
        String outfileDirS = "/Users/xujun/Desktop";
        File[] fs = new File(infileDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        int[] cnt=new int[sample.size()];
        nameSet.parallelStream().forEach(name -> {
            String infile1 = new File (infileDirS, name+"_TKD180602248_1.clean.fq.gz").getAbsolutePath();
            String infile2 = new File (infileDirS, name+"_TKD180602248_2.clean.fq.gz").getAbsolutePath();
            String outfile = new File (outfileDirS, "s20180622-LHY2_TKD180602248New.txt").getAbsolutePath();
            int count=0;int count1=0;
            String temp1=null;String seq1=null;
            String temp2=null;String seq2=null;  
            List<String> h1= new ArrayList();
            List<String> h2= new ArrayList();
            String sameName=null;int index=0;
            try {
                BufferedReader br1 = IOUtils.getTextGzipReader(infile1);
                BufferedReader br2 = IOUtils.getTextGzipReader(infile2);   
                BufferedWriter bw = IOUtils.getTextWriter(outfile);
                while ((temp1 = br2.readLine()) != null) {
                    count1++;
                    seq1=br2.readLine();
                    for(int i=4;i<8;i++){
//                        if(hm2.get(seq1.substring(0,i))==null){
                        if(getKeyList(hm2,seq1.substring(1,i)).isEmpty()){
                            if(i==8){
                                br1.readLine();br1.readLine();br1.readLine();br1.readLine();
                                count++;
                                break;
                            }    
                        }else{
                            h2=getKeyList(hm2,seq1.substring(1,i));
                            temp2=br1.readLine();
                            seq2=br1.readLine();
//                            if(hm1.get(hm2.get(seq1.substring(0,i)))!=null){
                            for(int a=4;a<8;a++){
                                if(getKeyList(hm1,seq2.substring(1,a)).isEmpty()){
                                    if(a==8){
                                        br1.readLine();br1.readLine();
                                        count++;
                                        break;
                                    } 
                                }else{
                                    h1=getKeyList(hm1,seq2.substring(1,a));
                                    if(!getTheSameSection(h1,h2).isEmpty()){
                                        sameName=getTheSameSection(h1,h2).toString().replace("[","").replace("]", "");
                                        index=sample.indexOf(sameName);
                                        cnt[index]++;
                                        break;
                                    }else{
                                        break;
                                    }   
                                }
                            }
                            br1.readLine();br1.readLine();
                            break;
                        }
                    }                    
                    br2.readLine();br2.readLine();
                }
                for(int i=0;i<cnt.length;i++){
                    bw.write(sample.get(i)+"\t"+cnt[i]);
                    bw.newLine();
                }
                System.out.println(count);
                System.out.println(count1);
                br1.close();br2.close();
                bw.flush();bw.close();
            }
            
            catch (Exception e) {
                System.out.println(sameName);
                System.out.println(index);
                e.printStackTrace();
            }

        });
    }
    public static List<String> getKeyList(HashMap<String,String> map,String value){//the method to get key from value
         List<String> keyList = new ArrayList();
         for(String getKey: map.keySet()){
             if(map.get(getKey).equals(value)){
                 keyList.add(getKey);
             }
         }
         return keyList;
     }
    public List getTheSameSection(List list1,List list2) {
       List resultList = new ArrayList();
       for (Object item : list2) {//遍历list1
            if (list1.contains(item)) {//如果存在这个数
            resultList.add(item);//放进一个list里面，这个list就是交集
        }
       }
       return resultList;
    }
}
