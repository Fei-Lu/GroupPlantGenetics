/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xiaohan.analysis.RNAseq;

import java.io.File;

/**
 *
 * @author yxh
 */
public class test {
    public void test(){
        this.testforfun();
    }
    public void testforfun(){
    String inputDirS = "/Users/yxh/Documents/RNA-seq/test";
    File[] fs = new File(inputDirS).listFiles();
    System.out.println(fs.toString());
    }
    public static void main(String[] args){
        new test();
        
    }
}
    

