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
                                    String taxaInfoDBFile, String outFilePGI, double threshForIndividualFd){
        List<File> fileList= IOUtils.getVisibleFileListInDir(fdByIndividualDir);
        Map<String, String> groupBySubspeciesMap= RowTableTool.getMap(taxaInfoDBFile, 23, 25);
        Map<String, String> groupBySubcontinentMap= RowTableTool.getMap(taxaInfoDBFile, 23, 24);
        try (BufferedWriter bw = IOTool.getWriter(outFilePGI)) {
            bw.write("IntrogressionID\tPGI\tP3\tSub\tGroupBySubspecies\tGroupBySubcontinent");
            bw.newLine();
            double[][] pgi;
            String introgressionID;
            String[] p3Array={"WE","DE","FT","AT"};
            Arrays.sort(p3Array);
            StringBuilder sb=new StringBuilder();
            String sub=null;
            for (int i = 0; i < fileList.size(); i++) {
                introgressionID=fileList.get(i).getName().substring(0,5);
                pgi=calculatePGI(fileList.get(i), coordinateFile, threshForIndividualFd);
                for (int j = 0; j < pgi.length; j++) {
                    for (int k = 0; k < pgi[j].length; k++) {
                        if (pgi[j][k] < 0) continue;
                        if (j==0 && k==0){
                            sub="D";
                        }else if (k==0){
                            sub="A";
                        }else if (k==1){
                            sub="B";
                        }
                        sb.setLength(0);
                        sb.append(introgressionID).append("\t").append(pgi[j][k]).append("\t");
                        sb.append(p3Array[j]).append("\t").append(sub).append("\t");
                        sb.append(groupBySubspeciesMap.get(introgressionID)).append("\t");
                        sb.append(groupBySubcontinentMap.get(introgressionID));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static double[][] calculatePGI(File fdByIndividualFile, String coordinateFile, double threshForIndividualFd){
        // WE DE FT AT
        String[] p3Array={"WE","DE","FT","AT"};
        double[][] pgiArray=new double[4][2];
        double genomeSize=getGenomeSizeFromCoordinateFile(coordinateFile);
        Arrays.sort(p3Array);
        try (BufferedReader br = IOTool.getReader(fdByIndividualFile)) {
            br.readLine();
            String line, p3, chrSub;
            List<String> temp;
            int start, end, p3Index, p3SubIndex=-1;
            double fd;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                chrSub=temp.get(0).substring(1,2);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                p3=temp.get(11).substring(0,2);
                p3Index=Arrays.binarySearch(p3Array, p3);
                fd=Double.parseDouble(temp.get(9));
                if (fd < threshForIndividualFd) continue;
                switch (chrSub){
                    case "A":
                        p3SubIndex=0;
                        break;
                    case "B":
                        p3SubIndex=1;
                        break;
                    case "D":
                        p3SubIndex=0;
                        break;
                }
                pgiArray[p3Index][p3SubIndex]+=fd * (end-start+1);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        double[][] pgi=new double[pgiArray.length][2];
        for (int i = 0; i < pgi.length; i++) {
            Arrays.fill(pgi[i],-1);
        }
        for (int i = 0; i < pgiArray.length; i++) {
            for (int j = 0; j < pgiArray[i].length; j++) {
                if (i==0 && j ==1) continue;
                pgi[i][j]=pgiArray[i][j]/genomeSize;
            }
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

    public static void start() {
        String fdByIndividualDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/006_fdResByIndividual/003_fdByIndividual.newMethod";
        String coordinateFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/006_fdResByIndividual/fdLoadBySubspecies100SNPwindow_50Step.txt";
        String taxaInfoDBFile="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String outFilePGI="/Users/xudaxing/RStudio/FdLoad/005_fdByIndividual.newMethod.PGI/PGI100SNPwindow_50Step" +
                ".newMethodBySub.thresh0.txt";
        double threshForIndividualFd=0.5;
        PGI.calculatePGI(fdByIndividualDir, coordinateFile, taxaInfoDBFile, outFilePGI, threshForIndividualFd);
    }
}
