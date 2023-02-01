package daxing.v2.localAncestryInfer.simulation;

import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SimulationMetadata {

    String[] demesID;
    String[] demesPath;
    String[] admixedPop;
    int[] admixedPopSampleSize;
    String[] nativePop;
    int[] nativePopSampleSize;
    List<String>[] introgressedPop;
    IntList[] introgressedPopSampleSize;
    int[] sequenceLen;
    double[] recombinationRate;

    public SimulationMetadata(String simulationMetadataFile){
        List<String> demesIDList = new ArrayList<>();
        List<String> demesPathList = new ArrayList<>();
        List<String> admixedPopList = new ArrayList<>();
        IntList admixedPopSampleSizeList = new IntArrayList();
        List<String> nativePopList = new ArrayList<>();
        IntList nativePopSampleSizeList = new IntArrayList();
        List<List<String>> introgressedPopList = new ArrayList<>();
        List<IntList> introgressedPopSampleSizeList = new ArrayList<>();
        IntList seqLenList = new IntArrayList();
        DoubleList recombinationRateList = new DoubleArrayList();
        try (BufferedReader br = IOTool.getReader(simulationMetadataFile)) {
            br.readLine();
            String line;
            List<String> temp, tem, te;
            IntList intList;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                demesIDList.add(temp.get(0));
                demesPathList.add(temp.get(1));
                tem = PStringUtils.fastSplit(temp.get(2), ":");
                admixedPopList.add(tem.get(0));
                admixedPopSampleSizeList.add(Integer.parseInt(tem.get(1)));
                tem = PStringUtils.fastSplit(temp.get(3), ":");
                nativePopList.add(tem.get(0));
                nativePopSampleSizeList.add(Integer.parseInt(tem.get(1)));
                tem = PStringUtils.fastSplit(temp.get(4), ":");
                te = PStringUtils.fastSplit(tem.get(0), ",");
                introgressedPopList.add(te);
                te = PStringUtils.fastSplit(tem.get(1),",");
                intList =new IntArrayList();
                for (int i = 0; i < te.size(); i++) {
                    intList.add(Integer.parseInt(te.get(i)));
                }
                introgressedPopSampleSizeList.add(intList);
                seqLenList.add(Integer.parseInt(temp.get(5)));
                recombinationRateList.add(Double.parseDouble(temp.get(6)));
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        this.demesID=demesIDList.toArray(new String[0]);
        this.demesPath = demesPathList.toArray(new String[0]);
        this.admixedPop=admixedPopList.toArray(new String[0]);
        this.admixedPopSampleSize = admixedPopSampleSizeList.toIntArray();
        this.nativePop=nativePopList.toArray(new String[0]);
        this.nativePopSampleSize=nativePopSampleSizeList.toIntArray();
        this.introgressedPop = introgressedPopList.toArray(new List[0]);
        this.introgressedPopSampleSize = introgressedPopSampleSizeList.toArray(new IntList[0]);
        this.sequenceLen=seqLenList.toIntArray();
        this.recombinationRate=recombinationRateList.toDoubleArray();
    }

    public String[] getNativePop() {
        return nativePop;
    }

    public String[] getAdmixedPop() {
        return admixedPop;
    }

    public double[] getRecombinationRate() {
        return recombinationRate;
    }

    public int[] getAdmixedPopSampleSize() {
        return admixedPopSampleSize;
    }

    public int[] getNativePopSampleSize() {
        return nativePopSampleSize;
    }

    public int[] getSequenceLen() {
        return sequenceLen;
    }

    public IntList[] getIntrogressedPopSampleSize() {
        return introgressedPopSampleSize;
    }

    public List<String>[] getIntrogressedPop() {
        return introgressedPop;
    }

    public String[] getDemesID() {
        return demesID;
    }

    public String[] getDemesPath() {
        return demesPath;
    }

    public int[] get_nWayAdmixture(){
        String[] demePathArray = this.getDemesPath();
        int[] nWayAdmixture = new int[demePathArray.length];
        Arrays.fill(nWayAdmixture, -1);
        for (int i = 0; i < demePathArray.length; i++) {
            nWayAdmixture[i] = Integer.parseInt(new File(demePathArray[i]).getName().substring(4,5));
        }
        return nWayAdmixture;
    }

    public int[] getTimeSinceAdmixture(){
        String[] demePathArray = this.getDemesPath();
        int[] timeSinceAdmixture = new int[demePathArray.length];
        Arrays.fill(timeSinceAdmixture, -1);
        List<String> temp;
        int index;
        for (int i = 0; i < demePathArray.length; i++) {
            temp = PStringUtils.fastSplit(new File(demePathArray[i]).getName().substring(4,5), "_");
            index = temp.indexOf("at");
            timeSinceAdmixture[i] = Integer.parseInt(temp.get(index+1));
        }
        return timeSinceAdmixture;
    }

    public int[] getChrID(){
        int[] chrID = new int[this.demesPath.length];
        Arrays.fill(chrID, 1);
        return chrID;
    }


}
