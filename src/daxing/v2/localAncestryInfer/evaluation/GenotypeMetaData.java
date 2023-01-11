package daxing.v2.localAncestryInfer.evaluation;

import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class GenotypeMetaData {

    String[] genotypeID;

    String[] genotypePath;

    int[] nWayAdmixture;

    String[] admixedPop;

    List<String>[] referencePopList;

    int[] timeSinceAdmixture;   // -1 means unknown

    int[] chrID;

    String[] recombinationMap;

    public GenotypeMetaData(String genotypePathFile){
        List<String> genotypeIDList = new ArrayList<>();
        List<String> genotypePathList=new ArrayList<>();
        IntList nWayAdmixtureList = new IntArrayList();
        List<String> admixedPopList = new ArrayList<>();
        List<List<String>> refPopList = new ArrayList<>();
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
                nWayAdmixtureList.add(Integer.parseInt(temp.get(2)));
                admixedPopList.add(temp.get(3));
                tem = PStringUtils.fastSplit(temp.get(4), ",");
                refPopList.add(tem);
                int time = Integer.parseInt(temp.get(5));
                timeSinceAdmixture.add(time < 0 ? -1 : time);
                chrIDList.add(Integer.parseInt(temp.get(6)));
                recombinationMap.add(temp.get(7));
            }
            br.close();
            this.genotypeID = genotypeIDList.toArray(new String[0]);
            this.genotypePath = genotypePathList.toArray(new String[0]);
            this.nWayAdmixture=nWayAdmixtureList.toIntArray();
            this.admixedPop = admixedPopList.toArray(new String[0]);
            this.referencePopList = refPopList.toArray(new List[0]);
            this.timeSinceAdmixture = timeSinceAdmixture.toIntArray();
            this.chrID=chrIDList.toIntArray();
            this.recombinationMap=recombinationMap.toArray(new String[0]);
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

    public int[] getChrID() {
        return chrID;
    }

    public String[] getRecombinationMap() {
        return recombinationMap;
    }
}
