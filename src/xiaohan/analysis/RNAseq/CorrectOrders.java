/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import pgl.infra.table.RowTable;

/**
 *
 * @author yxh
 */
public class CorrectOrders {
    public CorrectOrders(){
        this.Correct();
    }
    public void Correct(){
        String infile1 = "";
        String infile2 = "";
        HashMap NameStrain = new HashMap() ;
        HashMap SampleStrain = new HashMap() ;
        List<String> Size = new ArrayList<String>() ;
        RowTable<String> t1 = new RowTable<>(infile1);
        for (int i = 0 ; i < t1.getRowNumber();i++){
            NameStrain.put(t1.getCell(i, 7).substring(1),t1.getCell(i, 0));
            SampleStrain.put(t1.getCell(i, 4).substring(1),t1.getCell(i,0));
            Size.add(t1.getCell(i, 0).substring(1));
        }
        RowTable<String> t2 = new RowTable<>(infile2);
        for (int i = 0 ; i < t2.getRowNumber();i++){
            List<String> Sample = new ArrayList<String>();
            Sample.add(t2.getCell(i, 1).substring(1));
            String Sample1 = Sample.get(i).toString();
            if (SampleStrain.containsValue(Sample1)){
                System.out.print(Sample1);
                System.out.print(NameStrain.get(i));
            }
        }
    }
            
    public static void main(String[] args){
        new CorrectOrders();
        
    }
}
