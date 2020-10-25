package daxing.individualIntrogression;

import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * this class using fdByIndividual目录计算每个六倍体种质的PGI
 */
public class PGI {

    public static void calculatePGI(String fdByIndividualDir, String coordinateFile,
                                    String taxaInfoDBFile, String outFilePGI){
        List<File> fileList= IOUtils.getVisibleFileListInDir(fdByIndividualDir);
        Map<String, String> groupBySubspeciesMap= RowTableTool.getMap(taxaInfoDBFile, 23, 25);
        Map<String, String> groupBySubcontinentMap= RowTableTool.getMap(taxaInfoDBFile, 23, 24);
        try (BufferedWriter bw = IOTool.getWriter(outFilePGI)) {
            bw.write("IntrogressionID\tPGI\tP3\tGroupBySubspecies\tGroupBySubcontinent");
            bw.newLine();
            double[] pgi;
            String introgressionID;
            String[] p3Array={"WE","DE","FT","AT"};
            Arrays.sort(p3Array);
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < fileList.size(); i++) {
                introgressionID=fileList.get(i).getName().substring(0,5);
                pgi=calculatePGI(fileList.get(i), coordinateFile);
                for (int j = 0; j < pgi.length; j++) {
                    sb.setLength(0);
                    sb.append(introgressionID).append("\t").append(pgi[j]).append("\t");
                    sb.append(p3Array[j]).append("\t").append(groupBySubspeciesMap.get(introgressionID)).append("\t");
                    sb.append(groupBySubcontinentMap.get(introgressionID));
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double[] calculatePGI(File fdByIndividualFile, String coordinateFile){
        // WE DE FT AT
        String[] p3Array={"WE","DE","FT","AT"};
        double[] pgiArray=new double[4];
        double genomeSize=getGenomeSizeFromCoordinateFile(coordinateFile);
        Arrays.sort(p3Array);
        try (BufferedReader br = IOTool.getReader(fdByIndividualFile)) {
            br.readLine();
            String line, p3;
            List<String> temp;
            int start, end, p3Index;
            double fd;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                p3=temp.get(11).substring(0,2);
                fd=Double.parseDouble(temp.get(9));
                p3Index=Arrays.binarySearch(p3Array, p3);
                pgiArray[p3Index]+=fd * (end-start+1);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[] pgi=new double[pgiArray.length];
        for (int i = 0; i < pgiArray.length; i++) {
            pgi[i]=pgiArray[i]/genomeSize;
        }
        return pgi;
    }

    private static double getGenomeSizeFromCoordinateFile(String coordinateFile){
        double genomeSize=0;
        try (BufferedReader brCoordinate = IOTool.getReader(coordinateFile)) {
            brCoordinate.readLine();
            String line;
            List<String> temp;
            int start, end;
            while ((line= brCoordinate.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                genomeSize+=end-start+1;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return genomeSize;
    }
}
