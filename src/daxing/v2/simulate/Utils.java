package daxing.v2.simulate;

import daxing.common.utiles.IOTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

public class Utils {

    public static void addAncestral(String simulatedGenotype, String outFile){
        try (BufferedReader br = IOTool.getReader(simulatedGenotype);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            StringBuilder sb = new StringBuilder();
            sb.append(line).append("\t").append("ancestral");
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                sb.setLength(0);
                sb.append(line).append("\t").append(0);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
