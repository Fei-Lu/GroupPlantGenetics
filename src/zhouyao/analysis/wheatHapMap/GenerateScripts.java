/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class GenerateScripts {
    public GenerateScripts(){
        this.getRScripts();
    }
    public void getRScripts(){
         BufferedWriter bw = YaoIOUtils.getTextWriter("scripts.sh");
            try {
                for (int i = 10; i < 15;i++){
                   bw.write("nohup Rscript FullModle_AD.R "+ i+ " > ./log/AD_BLUPres_"+ i +".log &\n");
                }
                bw.flush();
                bw.close();
           
            } catch (IOException ex) {
                Logger.getLogger(WheatHmpEntraince.class.getName()).log(Level.SEVERE, null, ex);
            };
    }
}
