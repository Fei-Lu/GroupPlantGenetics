/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package xuebo.analysis.annotation;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import xuebo.analysis.data4CandChIA_PET.XueboIOUtils;

/**
 *
 * @author xuebozhao
 */
public class generateScripts {
    public generateScripts(){
        this.getScripts();
    }
    private void getScripts(){
        
        BufferedWriter bw;
        for (int i = 41; i < 50; i++){
            StringBuilder outfile = new StringBuilder();
            outfile.append(Integer.toString(i));
            outfile.append(".script");
            bw = XueboIOUtils.getTextWriter(outfile.toString());
            try {
                bw.write("#BSUB -n 8\n");
                bw.write("#BSUB -q low\n");
                bw.write("#BSUB -o "+ i +".out\n");
                bw.write("#BSUB -J " + i + "_bs\n\n");
                bw.write("bismark --bowtie2 -n 2 -L 25 -p 8 -o outputSRR5210"
                        + i +"_dir --gzip /public-supool/home/xbzhao/methylation/bismark "
                        + " -1 SRR5210"+i+"_1.fastq.gz -2 SRR5210"+ i + "_2.fastq.gz -N 0 > bismarkSRR5210" 
                                + i +".log &\n");        
                bw.flush();
                bw.close();
            } catch (IOException ex) {
                Logger.getLogger(generateScripts.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
    }
}
