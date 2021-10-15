package daxing.vmapII_1000.variantsAnnotation;

import daxing.common.utiles.IOTool;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.util.List;

public class Ancestral {

    public static void getAncestral(String estSFSInputDir, double threshold, String outFileName){
        List<File> files = IOTool.getFileListInDirEndsWith(estSFSInputDir, ".bedGraph.gz");
        try {
            BufferedWriter bw = IOTool.getWriter(outFileName);
            bw.write("Chr\tPos\tAncestral");
            bw.newLine();
            BufferedReader br;
            int total=0;
            for (int i = 0; i < files.size(); i++) {
                br = IOTool.getReader(files.get(i));
                String line, chr, pos;
                double p;
                List<String> temp, tem;
                TDoubleArrayList acgtCountList;
                StringBuilder sb = new StringBuilder();
                String[] acgt ={"A","C","G","T"};
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    p =Double.parseDouble(temp.get(6));
                    if (p < threshold) continue;
                    total++;
                    chr = temp.get(0).substring(3,5);
                    pos = temp.get(2);
                    tem = PStringUtils.fastSplit(temp.get(3), ",");
                    acgtCountList = new TDoubleArrayList();
                    for (int j = 0; j < tem.size(); j++) {
                        acgtCountList.add(Integer.parseInt(tem.get(j)));
                    }
                    int majorAlleleIndex=acgtCountList.indexOf(acgtCountList.max());
                    sb.setLength(0);
                    sb.append(chr).append("\t").append(pos).append("\t").append(acgt[majorAlleleIndex]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                br.close();
            }
            bw.flush();
            bw.close();
            System.out.println("Total "+total+ " snp have ancestral state");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
