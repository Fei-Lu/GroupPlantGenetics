package daxing.vmapII_1000.variantsAnnotation;

import daxing.common.IOTool;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class PhyloP {

    public static void split_toChrID(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, "phyloP_CON.bedg.gz");
        Map<Integer, BufferedWriter> chrBWMap = new HashMap<>();
        BufferedWriter bw = null;
        String chr;
        int[] chrIDArray;
        File outFile;
        try {
            for (int i = 0; i < files.size(); i++){
                chr = files.get(i).getName().substring(0,2);
                chrIDArray= RefV1Utils.getChrIDs();
                for (int e : chrIDArray){
                    outFile = new File(outDir, files.get(i).getName().replaceAll(chr,
                            "chr"+PStringUtils.getNDigitNumber(3
                                    , e)));
                    bw = IOTool.getWriter(outFile);
                    bw.write("ChrID\tPos\tPhyloP");
                    bw.newLine();
                    chrBWMap.put(e, bw);
                }
            }
            BufferedReader br;
            String line, refChr;
            List<String> temp;
            int chrIDEnd, posEnd, refEndPos ;
            double gerpScore;
            List<Integer> chrIDList;
            StringBuilder sb = new StringBuilder();
            for (File file: files){
                br = IOTool.getReader(file);
                chrIDList = new ArrayList<>();
                while ((line = br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    refChr = temp.get(0).substring(3,5);
                    refEndPos = Integer.parseInt(temp.get(2));
                    gerpScore = Double.parseDouble(temp.get(3));
                    chrIDEnd = RefV1Utils.getChrID(refChr, refEndPos);
                    posEnd = RefV1Utils.getPosOnChrID(refChr, refEndPos);
                    chrIDList.add(chrIDEnd);
                    bw = chrBWMap.get(chrIDEnd);
                    sb.setLength(0);
                    sb.append(chrIDEnd).append("\t").append(posEnd).append("\t");
                    sb.append(gerpScore);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
                for (Integer integer: chrIDList){
                    chrBWMap.get(integer).flush();
                    chrBWMap.get(integer).close();
                }
            }
//            for (Map.Entry<Integer, BufferedWriter> entry : chrBWMap.entrySet()){
//                entry.getValue().flush();
//                entry.getValue().close();
//            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
