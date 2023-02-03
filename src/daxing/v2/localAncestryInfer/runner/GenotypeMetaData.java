package daxing.v2.localAncestryInfer.runner;

import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.simulation.SimulationMetadata;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class GenotypeMetaData {

    String[] genotypeID;

    String[] genotypePath;

    String[] taxaInfoPath;

    int[] nWayAdmixture;

    String[] admixedPop;

    String[] nativePop;

    List<String>[] introgressedPop;

    List<String>[] referencePopList;

    int[] timeSinceAdmixture;   // -1 means unknown

    int[] chrID;

    String[] recombinationMap;

    /**
     * time since admixture is set to -1 by default
     * @param simulationMetadataPath
     * @param simulationOutDir
     */
    public GenotypeMetaData(String simulationMetadataPath, String simulationOutDir){
        SimulationMetadata simulationMetadata = new SimulationMetadata(simulationMetadataPath);
        this.genotypeID=simulationMetadata.getDemesID();
        String[] demeID = simulationMetadata.getDemesID();
        String[] genotypePath = new String[demeID.length];
        String[] taxaInfoPath = new String[demeID.length];
        String[] recombinationMap= new String[demeID.length];
        for (int i = 0; i < demeID.length; i++) {
            genotypePath[i] = new File(simulationOutDir, demeID[i]+".vcf").getAbsolutePath();
            taxaInfoPath[i] = new File(simulationOutDir, demeID[i]+".taxaInfo").getAbsolutePath();
            recombinationMap[i] = new File(simulationOutDir, demeID[i]+".recombinationMap").getAbsolutePath();
        }
        this.genotypePath=genotypePath;
        this.taxaInfoPath=taxaInfoPath;
        this.recombinationMap=recombinationMap;
        this.nWayAdmixture=simulationMetadata.get_nWayAdmixture();
        this.admixedPop= simulationMetadata.getAdmixedPop();
        this.nativePop= simulationMetadata.getNativePop();
        this.introgressedPop=simulationMetadata.getIntrogressedPop();
        int[] timeSinceAdmixture=new int[this.genotypeID.length];
        Arrays.fill(timeSinceAdmixture, -1);
        this.timeSinceAdmixture=timeSinceAdmixture;
        this.chrID=simulationMetadata.getChrID();
        List<String>[] referencePopList = new List[genotypeID.length];
        for (int i = 0; i < referencePopList.length; i++) {
            referencePopList[i] = new ArrayList<>();
            referencePopList[i].addAll(simulationMetadata.getIntrogressedPop()[i]);
            referencePopList[i].add(simulationMetadata.getNativePop()[i]);
        }
        this.referencePopList=referencePopList;
    }

    /**
     * time since admixture was supplied by genotypePathFile
     * @param genotypePathFile
     */
    public GenotypeMetaData(String genotypePathFile){
        List<String> genotypeIDList = new ArrayList<>();
        List<String> genotypePathList=new ArrayList<>();
        List<String> taxaInfoPathList = new ArrayList<>();
        IntList nWayAdmixtureList = new IntArrayList();
        List<String> admixedPopList = new ArrayList<>();
        List<String> nativePopList = new ArrayList<>();
        List<List<String>> introgressedPopList = new ArrayList<>();
        IntList timeSinceAdmixture=new IntArrayList();
        IntList chrIDList = new IntArrayList();
        List<String> recombinationMap = new ArrayList<>();
        String line;
        List<String> temp, tem;
        try (BufferedReader br = IOTool.getReader(genotypePathFile)) {
            br.readLine();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                genotypeIDList.add(temp.get(0));
                genotypePathList.add(temp.get(1));
                taxaInfoPathList.add(temp.get(2));
                nWayAdmixtureList.add(Integer.parseInt(temp.get(3)));
                admixedPopList.add(temp.get(4));
                nativePopList.add(temp.get(5));
                tem = PStringUtils.fastSplit(temp.get(6), ",");
                introgressedPopList.add(tem);
                int time = Integer.parseInt(temp.get(7));
                timeSinceAdmixture.add(time < 0 ? -1 : time);
                chrIDList.add(Integer.parseInt(temp.get(8)));
                recombinationMap.add(temp.get(9));
            }
            br.close();
            this.genotypeID = genotypeIDList.toArray(new String[0]);
            this.genotypePath = genotypePathList.toArray(new String[0]);
            this.nWayAdmixture=nWayAdmixtureList.toIntArray();
            this.admixedPop = admixedPopList.toArray(new String[0]);
            this.timeSinceAdmixture = timeSinceAdmixture.toIntArray();
            this.chrID=chrIDList.toIntArray();
            this.recombinationMap=recombinationMap.toArray(new String[0]);
            List<String>[] referencePopList = new List[genotypeID.length];
            for (int i = 0; i < referencePopList.length; i++) {
                referencePopList[i] = new ArrayList<>();
                referencePopList[i].addAll(introgressedPopList.get(i));
                referencePopList[i].add(nativePopList.get(i));
            }
            this.referencePopList=referencePopList;
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public int[] getnWayAdmixture() {
        return nWayAdmixture;
    }

    public String[] getGenotypePath() {
        return genotypePath;
    }

    public String[] getTaxaInfoPath() {
        return taxaInfoPath;
    }

    public int[] getTimeSinceAdmixture() {
        return timeSinceAdmixture;
    }

    public String[] getGenotypeID() {
        return genotypeID;
    }

    public List<String>[] getReferencePopList() {
        return referencePopList;
    }

    public String[] getAdmixedPop() {
        return admixedPop;
    }

    public String[] getNativePop() {
        return nativePop;
    }

    public List<String>[] getIntrogressedPop() {
        return introgressedPop;
    }

    public int[] getChrID() {
        return chrID;
    }

    public String[] getRecombinationMap() {
        return recombinationMap;
    }

    public String getGenotypeID(int indexOfRun){
        return this.getGenotypeID()[indexOfRun];
    }

    public String getGenotypePath(int indexOfRun){
        return this.getGenotypePath()[indexOfRun];
    }

    public TaxaInfo getTaxaInfo(int indexOfRun){
        return new TaxaInfo(this.getTaxaInfoPath()[indexOfRun]);
    }

    public int getNwayAdmixture(int indexOfRun){
        return this.getnWayAdmixture()[indexOfRun];
    }

    public String getAdmixedPop(int indexOfRun){
        return this.getAdmixedPop()[indexOfRun];
    }

    public String getNativePop(int indexOfRun){
        return this.getNativePop()[indexOfRun];
    }

    public List<String> getIntrogressedPop(int indexOfRun){
        return this.getIntrogressedPop()[indexOfRun];
    }

    public List<String> getReferencePopList(int indexOfRun){
        return this.getReferencePopList()[indexOfRun];
    }

    public int getTimeSinceAdmixture(int indexOfRun){
        return this.getTimeSinceAdmixture()[indexOfRun];
    }

    public int getChrID(int indexOfRun){
        return this.getChrID()[indexOfRun];
    }

    public String getRecombinationMap(int indexOfRun){
        return this.getRecombinationMap()[indexOfRun];
    }
}
