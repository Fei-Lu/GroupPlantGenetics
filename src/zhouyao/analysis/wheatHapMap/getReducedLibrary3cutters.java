/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.zip.GZIPOutputStream;

/**
 *
 * @author yaozhou
 */
public class getReducedLibrary3cutters {
    getReducedLibrary3cutters(String inFile,String cutter1,String cutter2,String cutter3, String outFile) {

        this.SecondFiltering(inFile,cutter1,cutter2,cutter3,outFile);
//        this.outputFiles(outfileS1, outfileS2);
    }

    public void SecondFiltering(String inFile, String cutter1,String cutter2,String cutter3,String outFile) {
        try{
            BufferedReader br ;
            if(inFile.endsWith("gz"))  br = YaoIOUtils.getTextGzipReader(inFile);
            else  br = YaoIOUtils.getTextReader(inFile);
            String outFileg = outFile + ".reduced.gz";
            String temp = null;
//            String temp1 = null;
//            String temp2 = "";
//            String temp3 = null;
            StringBuilder temp2 = new StringBuilder();
//            temp2.append("");
            int i = 0;
           
            try{
                BufferedWriter bw = YaoIOUtils.getTextGzipWriter(outFileg);
                bw.write("");
                bw.flush();  
                bw.close(); 
            }
            catch (Exception e){
                e.printStackTrace();
            }
            String outFiles = outFile +".stat.gz";
            try{
                BufferedWriter bws = YaoIOUtils.getTextGzipWriter(outFiles);
                bws.write("");
                bws.flush();  
                bws.close(); 
            }
            catch (Exception e){
                e.printStackTrace();
            }
            String chr ="";
            List <Integer> getPos = null;
//            long startTime,endTime;
            while (( temp = br.readLine()) != null) {
                if(temp.startsWith(">")){
                    System.out.println("Processing chromosome " + i + "...");
                    if(i > 0){
//                       String chr = Integer.toString(Integer.valueOf(temp.substring(1,2)) - 1);       
                        chr = Integer.toString(i);  
                        getPos = posSearching(temp2,cutter1,cutter2);
//                        startTime = System.currentTimeMillis();
                        getReducedLibrary(chr,getPos,temp2,outFiles,outFileg);
//                        endTime = System.currentTimeMillis();
//                        System.out.println("Testing time is "+(endTime-startTime)/1000+ " seconds");
//                       getReducedLibraryi = 0;
                    }
                    i++;
//                    String[] tem = temp.split(" ");
//                    temp1 = tem[0];                   
//                    StringBuilder temp2 = new StringBuilder();
//                    bw.write(tem[0] + "\n");
                    
                    temp2.delete(0, temp2.length());
                    
                }
                else {
                    temp2.append(temp);
                }
                
            }
            chr = Integer.toString(i);
            getPos = posSearching(temp2,cutter1,cutter2);
            getReducedLibrary(chr,getPos,temp2,outFiles,outFileg);
//            bw.write(temp2+"\n");
//            bw.flush();  
//            bw.close(); 
        }
        catch (Exception e){
            e.printStackTrace();
        }
    }
       
    public List<Integer> binaryCutter(StringBuilder inChr,String cutter){
        String query = null;
        List<Integer> pos = new ArrayList<>();
        for(int i = 0 ; i < inChr.length()-cutter.length() ; i++){
            if (i % 100000000 == 0){
                System.out.println("Processing "+i+" bp");
            }
            query = inChr.substring(i, i + cutter.length());
            if(query.equals(cutter)){
                pos.add(i);
            }       
        }
        return pos;
    }
    private List<Integer> posSearching(StringBuilder inChr,String cutter1,String cutter2){
        
        List<Integer> aa = new ArrayList();
//        System.out.println(inChr.length());
        aa = (binaryCutter(inChr,cutter1));
        
        List<Integer> bb = new ArrayList<>(this.binaryCutter(inChr, cutter2));
        
        List<Integer> cc = new ArrayList<>();
        
        List<Integer> getpos = new ArrayList<>();
        
//        String[] both = concact(pos1, pos2);    
//        int c[] = this.contact(pos1, pos2);

        cc.addAll(aa);
        cc.addAll(bb);
//        for(int i = 0; i < aa.size();i++){
//            System.out.println(aa.get(i));
//        }
//        
//        for(int i = 0; i < cc.size();i++){
//            System.out.println(cc.get(i));
//        }
//        Integer[] aa = pos1.toArray(new Integer[pos1.size()]);
//        Integer[] bb = pos2.toArray(new Integer[pos2.size()]);
//        Integer[] cc = pos3.toArray(new Integer[pos3.size()]);
        Collections.sort(cc);
       
//        for(int i = 0; i < cc.size();i++){
//            System.out.println(cc.get(i));
//        }
        boolean ispos1 = false;
        boolean ispos2 = false;
        int pos1V = 0;
        int pos2V = 0;
        int j = 0;
        
        for(int i = 0; i<cc.size();i++){
            if(j >= aa.size()){
                aa.add(-1);
            }
//            if(j<aa.size()){
                if(cc.get(i).equals(aa.get(j)) ){
                    ispos1 = true;  
                    pos1V = cc.get(i).intValue();
                    j++;
                    if(ispos1 && ispos2){  
                        if ((pos1V - pos2V)< 600 && (pos1V - pos2V)> 200){
//                            System.out.print(pos2V);
//                            System.out.print(pos1V);
                            getpos.add(pos2V);
                            getpos.add(pos1V+cutter1.length());
                            ispos2 = false; 
                        }
                    }
                }else{
                    ispos2 = true;
                    pos2V = cc.get(i).intValue();
                    if(ispos1 && ispos2){
                        if ((pos1V - pos2V)< 600 && (pos1V - pos2V)> 100){
                            getpos.add(pos1V);
                            getpos.add(pos2V+cutter2.length());
                            ispos1 = false; 
                        }
                    }
                } 
            }
//        }
//        for(int i = 0; i<getpos.size();i++){
//            System.out.println(getpos.get(i));
//        }
        return getpos;     
    }
     

    public void getReducedLibrary(String Chr,List<Integer> getpos,StringBuilder inChr,
            String outFiles,String outFileg){
        try{
        String cutterinChr = null;
        String headcutterinchr =null;
        String bed = "";
//        BufferedWriter bw = YaoIOUtils.getTextWriter(outfileS);
            for(int i = 0;i<getpos.size();i = i+2){
                cutterinChr = inChr.substring(getpos.get(i),getpos.get(i+1));
                bed = "Chr"+ Chr + "\t" + getpos.get(i) + "\t" + getpos.get(i+1);
                headcutterinchr = ">" + bed + "\n";
                getWriteStreamAppend(outFileg,headcutterinchr);
//                getWriteStreamAppend(outFileg,Chr);
//                getWriteStreamAppend(outFileg,"\t");
//                getWriteStreamAppend(outFileg,String.valueOf(getpos.get(i)));
//                getWriteStreamAppend(outFileg,"\t");
//                getWriteStreamAppend(outFileg,String.valueOf(getpos.get(i+1)));
//                getWriteStreamAppend(outFileg,"\n");
                getWriteStreamAppend(outFileg,cutterinChr+"\n");
//                getWriteStreamAppend(outFileg,"\n");
//                bw.write(cutterinChr);
//                bw.write(headcutterinchr);
               
                getWriteStreamAppend(outFiles,bed + "\t"+String.valueOf(getpos.get(i+1)-getpos.get(i))+"\n");
//                getWriteStreamAppend(outFiles,"\n");
            }
//            bw.flush();  
//            bw.close(); 
        }
        catch (Exception e){
            e.printStackTrace();
        }
        
        
    }
    
   
    public static void getWriteStreamAppend(String file, String conent) {                         
        BufferedWriter out = null;                                                   
        try {                                                                        
            out = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(
                    new FileOutputStream(file, true),65536)),65536);                              
            out.write(conent);                                                      
        } 
        catch (Exception e) {                                                     
            e.printStackTrace();                                                    
        } 
        finally {                                                                 
            try {                                                                    
                out.close();                                                        
            } catch (Exception e) {                                               
                e.printStackTrace();                                                
            }                                                                       
        }                                                                           
    }
}
