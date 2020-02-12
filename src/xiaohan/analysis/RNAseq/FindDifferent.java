/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;


import java.io.BufferedReader;
import static java.sql.DriverManager.println;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import pgl.infra.table.RowTable;

/**
 *
 * @author yxh
 */
public class FindDifferent {
    public FindDifferent(){
        this.FindDiff();
    }
    
    public void FindDiff(){
        String infile1 = "/Users/yxh/Desktop/Book2.txt";
        RowTable<String> t = new RowTable<>(infile1);
        List<String> FirstList = new ArrayList<String>() ;
        List<String> SecondList = new ArrayList<String>() ;
        List<String> Diff = new ArrayList<String>() ;
        for (int i = 0; i < t.getRowNumber(); i++) {
        FirstList.add(t.getCell(i, 0));
        SecondList.add(t.getCell(i, 1));
        String First = FirstList.get(i).toString();
        String Second = SecondList.get(i).toString();
            if (First.equals(Second)){}
                else {
                      Diff.add(First);  
            }
        String Different = Diff.toString();
        System.out.println(Different);
        }
        
    }
    
    public static void main(String[] args){
        new FindDifferent();
    }

}
