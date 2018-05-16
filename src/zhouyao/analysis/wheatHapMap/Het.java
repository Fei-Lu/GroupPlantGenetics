/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;

/**
 *
 * @author yaozhou
 */
public class Het {
    public Het(String inFile){
        this.getHet(inFile);
    }
    private void getHet(String inFile){
        BufferedReader br = YaoIOUtils.getTextReader(inFile);
        String outFile = inFile.replace(".vcf", ".het");
        BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
        
    }
}
