/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package zhouyao.analysis.wheatHapMap;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author yaozhou
 */
public class modifyPsl {
    public modifyPsl(String inFile){
        this.getModified(inFile);
    }
    private void getModified(String path){
        File all = new File (path);
        File[] fs = YaoIOUtils.listRecursiveFiles(all);
//        File[] fs = YaoIOUtils.listRecursiveFiles(new File(path));
        File[] subFiles = YaoIOUtils.listFilesEndsWith(fs, ".psl");
        for(File f : subFiles){
           
            BufferedReader br = YaoIOUtils.getTextReader(f.toString());
            String outFile = f.toString().replace(".psl",".modified.psl");
            BufferedWriter bw = YaoIOUtils.getTextWriter(outFile);
            String temp = null;
            String[] tname = outFile.split("_");
//            System.out.println(tname[tname.length-1]);
            String[] tN = tname[tname.length-1].split("\\.");
            String tNAME = tN[0];
             System.out.println(tNAME+"\n");
            try {
                while((temp = br.readLine())!=null){
                    if(temp.startsWith("#")){
                        bw.write(temp);
                        bw.flush();
                        bw.newLine();
                    }else{
                        String[] tem = temp.split("\t");
                        tem[13] = tNAME; //tName
                        tem[14] = Integer.toString(Integer.parseInt(tem[14]) - 4); // tSize
                        tem[15] = Integer.toString(Integer.parseInt(tem[15]) - 4); // tStart
                        tem[16] = Integer.toString(Integer.parseInt(tem[16]) - 4); //tEnd
                        String[] tStarts = tem[20].split(",");
                        StringBuilder st = new StringBuilder();
                        for(int i = 0; i< tStarts.length; i++){
                            tStarts[i]=Integer.toString(Integer.parseInt(tStarts[i])-4);
                            st.append(tStarts[i]);
                            st.append(",");
                        }
                        tem[20] = st.toString(); // tStarts
                        StringBuilder newTem = new StringBuilder();
                        for(int i = 0; i< tem.length-1; i++){
                           newTem.append(tem[i]);
                           newTem.append("\t");
                        }
                        newTem.append(tem[20]);
                        bw.write(newTem.toString());
                        bw.newLine();
                        bw.flush();
                    }
                }
                bw.flush();
                bw.close();
            } catch (IOException ex) {
                ex.printStackTrace();
            }
        }
    }
    
}
