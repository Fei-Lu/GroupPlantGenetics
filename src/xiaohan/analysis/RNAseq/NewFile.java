/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import pgl.infra.table.RowTable;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import pgl.infra.utils.IOUtils;

/**
 *
 * @author yxh
 */
public class NewFile {
    String dir = "/Users/yxh/Documents/RNA-seq/test";
    RowTable<String> t = new RowTable<>("/Users/yxh/programs/RNA-seq-root-20180315-1.txt");
    
    public NewFile(){
        this.mkfile();
    }
    public void creat(){
    }
//    public void creatNewFile(String dir,RowTable t){
//        for(int i=0;i<t.getRowNumber();i++){
//         File.(dir,t.getCell(i, 0).toString()+"_R1.fq");
//            }
//        for(int i=0;i<t.getRowNumber();i++){
//            
//            File file = new File(dir,t.getCell(i, 0).toString()+"_R2.fq");
//            
//            }
//        }
    public void mkfile() {
        try{
            for(int i=0;i<t.getRowNumber();i++){
                String name = t.getCell(i, 0).toString() + ".rate";
                String infileS = new File(dir,name).getAbsolutePath();
                BufferedWriter bw = IOUtils.getTextWriter(infileS);
                bw.flush(); bw.close();
            }
        }
        catch(Exception e){
            e.printStackTrace();
            System.exit(1);
        }
    }
    
    public static void main(String[] main){
        new NewFile();
    }
    }

