/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package aoyue.analysis.WheatGeneticLoad;

import pgl.infra.table.RowTable;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import pgl.infra.utils.IOFileFormat;

/**
 *
 * @author Aoyue
 */
public class MapMake {

    public MapMake() {
        //this.sortRs_ID();
        //this.countCaseInLuLab();
        //this.countCaseInJiaoLab();
        //this.countCaseInLingLab();
        this.countCatgory();
    }
    
    public void countCatgory(){
        String infileS = "/Users/Aoyue/Documents/sum.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(4);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
    public void sortRs_ID(){
        String infileS = "/Users/Aoyue/Documents/Data/project/wheatVMapII/003_dataAnalysis/test.txt";
        String outfileS = "/Users/Aoyue/Documents/ABD_ltest.txt";
        RowTable<String> t = new RowTable<>(infileS);
        t.sortAsNumber(0);
        t.writeTextTable(outfileS, IOFileFormat.Text);
    }
    
    public void countCaseInLuLab(){
        String infileS = "/Users/Aoyue/Documents/ABDLuLabmap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(3);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
    public void countCaseInJiaoLab(){
        String infileS = "/Users/Aoyue/Documents/ABDJiaoLabmap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(3);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
    public void countCaseInLingLab(){
        String infileS = "/Users/Aoyue/Documents/ABDLingLabmap.txt";
        RowTable<String> t = new RowTable<>(infileS);
        List<String> l = t.getColumn(3);
        System.out.println(l.size() + " list个数");
        Set<String> s = new HashSet<>(l);
        System.out.println(s.size() + " set个数");
        System.out.println(s);
        for(String a : s){
            System.out.println(a + "    " + Collections.frequency(l, a));
        }
    }
    
}
