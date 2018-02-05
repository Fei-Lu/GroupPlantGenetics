/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xujun.analysis.rnaseq;

/**
 *
 * @author Jun Xu
 */
import com.koloboke.collect.map.hash.HashByteByteMap;
import format.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import jxl.read.biff.BiffException;
import utils.IOUtils;
public class DemoSample{
    public DemoSample() throws IOException, FileNotFoundException, BiffException {

          this.test();
    }
public void test () throws IOException {
        String barcodeFileS = "E:\\experimental data\\RNA_seq\\3'RNAseq barcode.txt";
        String inputDirS ="E:\\experimental data\\RNA_seq\\twice\\clean_data\\";
        String outputDirS = "E:\\experimental data\\RNA_seq\\twice\\small";        
        RowTable<String> rt = new RowTable<>(barcodeFileS);
        int rowNumber = rt.getRowNumber();
        int columnNumber = rt.getColumnNumber();
        HashMap<String, Integer> barcodeIndexMap = new HashMap<>();
        List<String> nameList = new ArrayList<>();
        for (int i = 0; i < rt.getRowNumber(); i++) {
            barcodeIndexMap.put(rt.getCell(i, 1), i);
            nameList.add(rt.getCell(i, 0));
        }
        File[] fs = new File(inputDirS).listFiles();
        HashSet<String> nameSet = new HashSet();
        for (int i = 0; i < fs.length; i++) {
            if (fs[i].isHidden()) continue;
            nameSet.add(fs[i].getName().split("_")[0]);
        }
        
        nameSet.stream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_1.clean.fq").getAbsolutePath();
            String seq = null;
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);                               
                BufferedWriter[] bws = new BufferedWriter[rowNumber];//这里只是创建了一个bufferedreader类型的数组 并没有对里面的各个new
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)).getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                String temp=null;
                String seq1=null;
                //String seq=null;
                int cnt=0;
                while ((temp = br1.readLine()) != null){
                    cnt++;                   
                    seq=br1.readLine(); 
                    Integer index = barcodeIndexMap.get(seq.substring(0, 8));
                    if (index == null) {
                        br1.readLine();br1.readLine();
                        continue;
                    }
/*                    int i=0;//用charat的方法 运行31s
                    for(int a=9;a<seq.length();a++){
                        if(seq.charAt(a)=='T'){
                            i++;
                        }else{
                            if(i<8){
                                i++;
                            }else{
                                seq=seq.substring(i+8);
                                break;
                            }
                        }
                    }*/
//                    System.out.println(i);
//                   System.out.println(seq);
/*                    seq=seq.substring(8);//每一次都取4个碱基变成byte类型 速度很慢 两次都是4min
                    for(int i=0;i+4<seq.length();i=i+4){
                        byte[] oriByte = seq.substring(i,i+4).getBytes();
                        HashByteByteMap ascMap = format.dna.BaseEncoder.getAscIIByteMap();
                        byte[] tranByte = new byte[oriByte.length];
                        for (int a = 0; a < oriByte.length; a++) {
                            tranByte[a] = ascMap.get(oriByte[a]);
                        }
                        int v=0;
                        for (int b = 0; b < tranByte.length; b++) {            
                            v = (v << 2) + tranByte[b];
                        }
                        if(v!=255){
                            seq=seq.substring(i);
                            break;
                        }
                    }*/
                    seq=seq.substring(8);//除barcode以外的其他碱基全都变成byte 每次取8个bite 速度要快很多 1min以内
                    byte[] oriByte = seq.substring(0).getBytes();
                    HashByteByteMap ascMap = format.dna.BaseEncoder.getAscIIByteMap();
                    byte[] tranByte = new byte[oriByte.length];
                    for (int a = 0; a < oriByte.length; a++) {
                        tranByte[a] = ascMap.get(oriByte[a]);
                    }
                    for(int i=0;i+2<seq.length();i+=4){                        
                        int v=0;
                        for (int b = i; b < i+4; b++) {            
                            v = (v << 2) + tranByte[b];
                        }
                        if(v!=255){
                                seq1=seq.substring(i);
                                break;                           
                        }
                    }
                    if(seq==seq1){
                        seq1=" ";
                    }

                    bws[index].write(temp + "\n");
                    bws[index].write(seq1 + "\n");
                    bws[index].write(br1.readLine() + "\n");
                    bws[index].write(br1.readLine()+"\n");
                    
                }
                br1.close();
                for (int i = 0; i < bws.length; i++) {
                    bws[i].flush();
                    bws[i].close();
                }
                    
                       
        }
                
        catch (Exception ex) {
               System.out.println(seq+"\t1234");        
               ex.printStackTrace();
               
        }
        });      
     }
    public static void main(String[] args) throws IOException, FileNotFoundException, BiffException{
        new DemoSample();
    }
}




