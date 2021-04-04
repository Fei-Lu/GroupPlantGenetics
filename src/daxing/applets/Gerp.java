package daxing.applets;

import com.ibm.icu.text.NumberFormat;
import daxing.common.IOTool;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class Gerp {

    public static void calculateConstraintRegionSize(String inputGERPDir, String outFile){
        List<File> files= IOTool.getVisibleDir(inputGERPDir);
        TDoubleArrayList gerpRes=new TDoubleArrayList();
        for (File file: files){
            gerpRes.add(calculateGERPGreaterThanZeroCount(file.getAbsolutePath()));
        }
        double[] gerpArray= gerpRes.toArray();
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setGroupingUsed(false);
        numberFormat.setMinimumFractionDigits(5);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < gerpArray.length; i++) {
                sb.setLength(0);
                sb.append("chr").append(PStringUtils.getNDigitNumber(3, i+1)).append("\t").append(gerpArray[i]);
                sb.append("\t").append(RefV1Utils.getChrIDLength(i+1)).append("\t");
                sb.append(numberFormat.format(gerpArray[i]/RefV1Utils.getChrIDLength(i+1)));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int calculateGERPGreaterThanZeroCount(String gerpFile){
        int count=0;
        try (BufferedReader br = IOTool.getReader(gerpFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                if (Double.parseDouble(temp.get(2)) > 0){
                    count++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return count;
    }
}
