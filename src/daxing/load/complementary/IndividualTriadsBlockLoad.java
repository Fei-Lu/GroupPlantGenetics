package daxing.load.complementary;

import daxing.common.LoadType;
import daxing.common.NumberTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IndividualTriadsBlockLoad {

    TriadsBlock triadsBlocks;
    List<String>[] genotypedGeneList;
    //    int[] blockCDSTotalLen;
    int[][] subLoadInfo_count; //dim 1 2: subgenome, loadInfo
    double[][] abdLoadType_load; //dim 1 2: subgenome, loadType
    /**
     * -1: syn derived count is 0
     */
    double[] stronglyLoad; // dim 1: subgenome
    double[] slightlyLoad; // dim 1: subgenome

    public IndividualTriadsBlockLoad(TriadsBlock triadsBlock, List<String> geneNameList, List<int[]> loadInfoList, TIntArrayList[] abdIndexArray){
        this.triadsBlocks=triadsBlock;
        this.genotypedGeneList=new List[WheatLineage.values().length];
        this.subLoadInfo_count=new int[WheatLineage.values().length][];
        this.abdLoadType_load=new double[WheatLineage.values().length][];
        for (int i = 0; i < WheatLineage.values().length; i++) {
            this.genotypedGeneList[i]=new ArrayList<>();
            this.subLoadInfo_count[i]=new int[9];
            this.abdLoadType_load[i]=new double[LoadType.values().length];
            Arrays.fill(this.subLoadInfo_count[i], -1);
            Arrays.fill(this.abdLoadType_load[i], -1);
        }
        this.stronglyLoad=new double[WheatLineage.values().length];
        this.slightlyLoad=new double[WheatLineage.values().length];
        Arrays.fill(this.stronglyLoad, -1);
        Arrays.fill(this.slightlyLoad, -1);
        this.mergeLoadInfo(geneNameList, loadInfoList, abdIndexArray);
    }

    public TriadsBlock getTriadsBlocks() {
        return triadsBlocks;
    }

    public List<String>[] getGenotypedGeneNumArray() {
        return genotypedGeneList;
    }

    public int[] getGenotypedGeneList(){
        List<String>[] genotypedGeneNumArray=this.getGenotypedGeneNumArray();
        int[] genotypedGeneNum=new int[WheatLineage.values().length];
        for (int i = 0; i < genotypedGeneNumArray.length; i++) {
            genotypedGeneNum[i]=genotypedGeneNumArray[i].size();
        }
        return genotypedGeneNum;
    }

    public int[][] getSubLoadInfo_count() {
        return subLoadInfo_count;
    }

    public double[][] getAbdLoadType_load() {
        return abdLoadType_load;
    }

    public double[] getSlightlyLoad() {
        return slightlyLoad;
    }

    public double[] getStronglyLoad() {
        return stronglyLoad;
    }

    private void mergeLoadInfo(List<String> geneNameList, List<int[]> loadInfoList, TIntArrayList[] abdIndexArray){
        for (int i = 0; i < WheatLineage.values().length; i++) {
            for (int j = 0; j < abdIndexArray[i].size(); j++) {
                this.genotypedGeneList[i].add(geneNameList.get(abdIndexArray[i].get(j)));
            }
        }
        int[] loadInfo;
        int[] loadInfoCum;
        for (int i = 0; i < abdIndexArray.length; i++) {
            loadInfoCum=new int[9];
            for (int j = 0; j < abdIndexArray[i].size(); j++) {
                loadInfo=loadInfoList.get(abdIndexArray[i].get(j));
                for (int k = 0; k < loadInfo.length; k++) {
                    loadInfoCum[k]+=loadInfo[k];
                }
            }
            this.subLoadInfo_count[i]=loadInfoCum;
        }
        this.calculateLoad();
    }

    private void calculateLoad(){
        int[][] subLoadInfo_count=this.getSubLoadInfo_count();
        double[] loadType_load;
        int[] loadInfo;
        for (int i = 0; i < subLoadInfo_count.length; i++) {
            loadInfo=subLoadInfo_count[i];
            loadType_load=new double[3];
            Arrays.fill(loadType_load, -1);
            for (int j = 0; j < LoadType.values().length; j++) {
                if (loadInfo[j*3]==0) continue;
                loadType_load[j]=(loadInfo[j*3+1]+0.5*loadInfo[j*3+2])/loadInfo[j*3];
            }
            this.abdLoadType_load[i]=loadType_load;
        }
        for (int i = 0; i < abdLoadType_load.length; i++) {
            if (abdLoadType_load[i][0] <= 0) continue;
            if (abdLoadType_load[i][1] < 0) continue;
            if (abdLoadType_load[i][2] < 0) continue;
            slightlyLoad[i]=abdLoadType_load[i][1]/abdLoadType_load[i][0];
            stronglyLoad[i]=abdLoadType_load[i][2]/abdLoadType_load[i][0];
        }
    }

    /**
     *
     * @return slightly and strongly
     */
    public double[] calculateAdditiveLoad(){
        double[] slightlyLoad=this.getSlightlyLoad();
        double[] stronglyLoad=this.getStronglyLoad();
        double[] slightlyStrongly=new double[2];
        slightlyStrongly[0]= Arrays.stream(slightlyLoad).sum();
        slightlyStrongly[1]=Arrays.stream(stronglyLoad).sum();
        return slightlyStrongly;
    }

    /**
     *
     * @return slightly and strongly
     */
    public double[] calculateDominanceLoad(){
        double[] slightlyLoad=this.getSlightlyLoad();
        double[] stronglyLoad=this.getStronglyLoad();
        double[] slightlyStrongly=new double[2];
        slightlyStrongly[0]=Arrays.stream(slightlyLoad).min().getAsDouble();
        slightlyStrongly[1]=Arrays.stream(stronglyLoad).min().getAsDouble();
        return slightlyStrongly;
    }

    public String toString(){
        StringBuilder sb=new StringBuilder();
        List<String>[] blockGeneName=this.getTriadsBlocks().getBlockGeneName();
        int[] abdGenotypedGeneNum=this.getGenotypedGeneList();
        double[] slightlyLoad=this.getSlightlyLoad();
        double[] stronglyLoad=this.getStronglyLoad();
        sb.append(this.getTriadsBlocks().getTriadsID()).append("\t");
        for (int i = 0; i < blockGeneName.length; i++) {
            sb.append(blockGeneName[i].size()).append(",");
        }
        sb.deleteCharAt(sb.length()-1).append("\t");
        for (int i = 0; i < abdGenotypedGeneNum.length; i++) {
            sb.append(abdGenotypedGeneNum[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1).append("\t");
        for (int i = 0; i < slightlyLoad.length; i++) {
            if (slightlyLoad[i] < 0){
                sb.append("NA").append(",");
            }else {
                sb.append(NumberTool.format(slightlyLoad[i], 5)).append(",");
            }
        }
        sb.deleteCharAt(sb.length()-1).append("|");
        for (int i = 0; i < stronglyLoad.length; i++) {
            if (stronglyLoad[i] < 0){
                sb.append("NA").append(",");
            }else {
                sb.append(NumberTool.format(stronglyLoad[i], 5)).append(",");
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }
}
