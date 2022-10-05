package daxing.v2.ancestryHmm;

import daxing.common.utiles.IOTool;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class Split {

    public static void splitPanelToChr(String panelFile, String outDir){
        try (BufferedReader br = IOTool.getReader(panelFile)) {
            BufferedWriter bw = null;
            String line;
            List<String> temp;
            String currentChr = "NA";
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line.substring(0,30));
                if (currentChr.equals("NA")){
                    currentChr = temp.get(0);
                    bw = IOTool.getWriter(new File(outDir, "chr"+temp.get(0)+"_vmap2.1_onlyGenotype_haploid_pruned" +
                            ".panel"));
                }
                if (currentChr.equals(temp.get(0))){
                    bw.write(line);
                    bw.newLine();
                }else {
                    bw.flush();
                    bw.close();
                    currentChr = temp.get(0);
                    bw = IOTool.getWriter(new File(outDir, "chr"+temp.get(0)+"_vmap2.1_onlyGenotype_haploid_pruned" +
                            ".panel"));
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
