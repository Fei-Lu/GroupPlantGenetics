package daxing.v2.loter;

import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.text.NumberFormat;

public class LoterMisc {

    /**
     *
     * @param inputFile 有header
     * @param rate such as 0.01
     * @param subsetFile 有header
     */
    public static void getSubsetFromFile_loterAncestryProportionPerSite(String inputFile, double rate, String subsetFile){
        long start=System.nanoTime();
        BufferedReader br;
        BufferedWriter bw;
        try {
            if (inputFile.endsWith("gz")){
                br= IOUtils.getTextGzipReader(inputFile);
                bw=IOUtils.getTextGzipWriter(subsetFile);
            }else {
                br=IOUtils.getTextReader(inputFile);
                bw=IOUtils.getTextWriter(subsetFile);
            }
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            String line;
            double r;
            int count=0;
            int total=0;
            NumberFormat numberFormat=NumberFormat.getInstance();
            numberFormat.setGroupingUsed(true);
            int chrID = Integer.parseInt(new File(inputFile).getName().substring(3,6));
            String subgenome = RefV1Utils.getSubgenomeFromChrID(chrID);
            while ((line=br.readLine())!=null){
                total++;
                r=Math.random();
                if (r>rate){
                    if (subgenome.equals("A") || subgenome.equals("B")){
                        br.readLine();
                        br.readLine();
                        continue;
                    }else {
                        br.readLine();
                        continue;
                    }
                }
                count++;
                if (subgenome.equals("A") || subgenome.equals("B")){
                    bw.write(line);
                    bw.newLine();
                    line = br.readLine();
                    bw.write(line);
                    bw.newLine();
                    line=br.readLine();
                    bw.write(line);
                    bw.newLine();
                }else {
                    bw.write(line);
                    bw.newLine();
                    line=br.readLine();
                    bw.write(line);
                    bw.newLine();
                }
            }
            br.close();
            bw.flush();
            bw.close();
            System.out.println("samping "+numberFormat.format(count)+"("+numberFormat.format(total)+") row from "
                    +new File(inputFile).getName()+" into "+new File(subsetFile).getName()+" in "
                    + Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
