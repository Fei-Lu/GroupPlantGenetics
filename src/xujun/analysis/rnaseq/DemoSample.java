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
import xuebo.analysis.annotation.IOUtils;
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
        
        nameSet.parallelStream().forEach(p -> {
            String infile1 = new File (inputDirS, p+"_1.clean.fq").getAbsolutePath();
            String infile2 = new File (inputDirS, p+"_2.clean.fq").getAbsolutePath();
            String output=new File(outputDirS,"wrong").getAbsolutePath();
            String seq = null;
            String seq2=null;
            try {
                BufferedReader br1 = utils.IOUtils.getTextReader(infile1);  
                BufferedReader br2 = utils.IOUtils.getTextReader(infile2); 
                BufferedWriter[] bws = new BufferedWriter[rowNumber];
                BufferedWriter bw=utils.IOUtils.getTextWriter(output);
                for (int i = 0; i < bws.length; i++) {
                    String outfileS = new File (outputDirS, nameList.get(i)).getAbsolutePath();
                    bws[i] = IOUtils.getTextWriter(outfileS);
                }
                String temp=null; 
                String temp2=null;
                int i=0;
                while ((temp = br1.readLine()) != null&i<100){
                    i+=1;
                    seq=br1.readLine();                    
                    Integer index = barcodeIndexMap.get(seq.substring(0, 8));
                    br1.readLine();br1.readLine();                    
                    if (index == null) {                        
                        bw.write(br2.readLine());bw.write(br2.readLine());bw.write(br2.readLine());bw.write(br2.readLine());
                        continue;
                    }else{
                        temp2=br2.readLine();
                        seq2=br2.readLine();
                        int count=0;
                        for(int a=0;a<seq2.length();a++){
                            if(seq2.charAt(a)=='A'){
                                count++;
                            }else{
                                count=0;
                            }
                            if(count>7){
                                seq2=seq2.substring(0, a-7);
                            }
                        }                    
                        bws[index].write(temp2 + "\n");
                        bws[index].write(seq2 + "\n");
                        bws[index].write(br2.readLine() + "\n");
                        bws[index].write(br2.readLine()+"\n");
                    }                    
                }
                br1.close();
                br2.close();
                    
                       
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




