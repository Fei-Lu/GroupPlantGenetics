package daxing.load.ancestralSite;

import daxing.common.ArrayTool;
import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class TriadLoad {

    String[] triadID;
    double[][] normalizedDerivedNum;  // synA synB synD nonsynA nonsynB nonsynD delA delB delD
    String[][] region;                //  synRegion nonsynRegion delRegion

    public TriadLoad(String triadListDir){
        List<File> files= IOUtils.getVisibleFileListInDir(triadListDir);
        Set<String> triadIDSet=new HashSet<>();
        for (int i = 0; i < files.size(); i++) {
            triadIDSet.addAll(RowTableTool.getColumnSet(files.get(i).getAbsolutePath(), 0));
        }
        this.triadID=triadIDSet.stream().sorted().toArray(String[]::new);
        normalizedDerivedNum=new double[triadID.length][];
        for (int i = 0; i < normalizedDerivedNum.length; i++) {
            normalizedDerivedNum[i]=new double[9];
        }
        region=new String[triadID.length][];
        for (int i = 0; i < region.length; i++) {
            region[i]=new String[3];
        }
    }

    private int binarySearch(String triadID){
        return Arrays.binarySearch(this.triadID, triadID);
    }

    private void add(int triadIDIndex, double[] normalizedDerivedNum){
        double[] sum=ArrayTool.add(this.normalizedDerivedNum[triadIDIndex], normalizedDerivedNum);
        this.normalizedDerivedNum[triadIDIndex]=sum;
    }

    public void mergeMutipleTriad(List<File> normalizedTriadFiles, String outFile){
        BufferedReader br;
        try {
            for (int i = 0; i < normalizedTriadFiles.size(); i++) {
                br = IOTool.getReader(normalizedTriadFiles.get(i));
                br.readLine();
                String line;
                List<String> temp;
                int index;
                double[] normalizedDerivedNum;
                while ((line=br.readLine())!=null){
                    temp= PStringUtils.fastSplit(line);
                    index=binarySearch(temp.get(0));
                    normalizedDerivedNum=new double[9];
                    normalizedDerivedNum[0]=Double.parseDouble(temp.get(1));
                    normalizedDerivedNum[1]=Double.parseDouble(temp.get(2));
                    normalizedDerivedNum[2]=Double.parseDouble(temp.get(3));
                    normalizedDerivedNum[3]=Double.parseDouble(temp.get(5));
                    normalizedDerivedNum[4]=Double.parseDouble(temp.get(6));
                    normalizedDerivedNum[5]=Double.parseDouble(temp.get(7));
                    normalizedDerivedNum[6]=Double.parseDouble(temp.get(9));
                    normalizedDerivedNum[7]=Double.parseDouble(temp.get(10));
                    normalizedDerivedNum[8]=Double.parseDouble(temp.get(11));
                    this.add(index, normalizedDerivedNum);
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        average(normalizedTriadFiles.size());
        determineRegion();
        write(outFile);
    }

    private void average(int fileNum){
        double[][] normalizedDerivedNum=this.normalizedDerivedNum;
        double[][] averaged=new double[normalizedDerivedNum.length][];
        for (int i = 0; i < averaged.length; i++) {
            averaged[i]=new double[9];
        }
        for (int i = 0; i < normalizedDerivedNum.length; i++) {
            for (int j = 0; j < normalizedDerivedNum[i].length; j++) {
                averaged[i][j]= NumberTool.format(normalizedDerivedNum[i][j]/fileNum, 5);
            }
        }
        this.normalizedDerivedNum=averaged;
    }

    private void determineRegion(){
        double[][] normalizedDerivedNum=this.normalizedDerivedNum;
        double[] synDerivedNum, nonsynDerivedNum, delDerivedNum;
        String[] region;
        for (int i = 0; i < normalizedDerivedNum.length; i++) {
            synDerivedNum=Arrays.copyOfRange(normalizedDerivedNum[i], 0,3);
            nonsynDerivedNum=Arrays.copyOfRange(normalizedDerivedNum[i], 3,6);
            delDerivedNum=Arrays.copyOfRange(normalizedDerivedNum[i], 6,9);
            region=new String[3];
            region[0]=Standardization.getNearestPointIndex(synDerivedNum).getRegion();
            region[1]=Standardization.getNearestPointIndex(nonsynDerivedNum).getRegion();
            region[2]=Standardization.getNearestPointIndex(delDerivedNum).getRegion();
            this.region[i]=region;
        }
    }

    private void write(String outFile){
        try (BufferedWriter bw = IOTool.getTextGzipWriter(outFile)) {
            bw.write("TriadID\tnormalizedNumDerivedInSynA\tnormalizedNumDerivedInSynB\tnormalizedNumDerivedInSynD" +
                    "\tsynRegion\tnormalizedNumDerivedInNonsynA\tnormalizedNumDerivedInNonsynB" +
                    "\tnormalizedNumDerivedInNonsynD\tnonsynRegion\tnormalizedNumDerivedInHGDeleteriousA" +
                    "\tnormalizedNumDerivedInHGDeleteriousB\tnormalizedNumDerivedInHGDeleteriousD\tdelRegion");
            bw.newLine();
            StringBuilder sb;
            double[] normalizedDerivedNum;
            String[] region;
            for (int i = 0; i < this.triadID.length; i++) {
                normalizedDerivedNum=this.normalizedDerivedNum[i];
                region=this.region[i];
                sb=new StringBuilder();
                sb.append(triadID[i]).append("\t");
                for (int j = 0; j < normalizedDerivedNum.length; j++) {
                    sb.append(normalizedDerivedNum[j]).append("\t");
                    if (j==2){
                        sb.append(region[0]).append("\t");
                    }
                    if (j==5){
                        sb.append(region[1]).append("\t");
                    }
                    if (j==8){
                        sb.append(region[2]).append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


}
