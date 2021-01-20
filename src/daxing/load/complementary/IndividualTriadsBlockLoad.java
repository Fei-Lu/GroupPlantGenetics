package daxing.load.complementary;

import daxing.common.*;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class IndividualTriadsBlockLoad {

    TriadsBlock triadsBlocks;
    List<String>[] genotypedGeneList;
    //    int[] blockCDSTotalLen;
    int[][] subLoadInfo_count; //dim 1 2: subgenome, loadInfo
    double[][] abdLoadType_load; //dim 1 2: subgenome, loadType

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

    /**
     * NA: missing
     * NaN: 0/0
     * @return
     */
    public String[] getSlightlyLoad() {
        String[] slightlyLoad=new String[WheatLineage.values().length];
        Arrays.fill(slightlyLoad, "NA");
        NumberFormat numberFormat = NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(5);
        for (int i = 0; i < abdLoadType_load.length; i++) {
            if (abdLoadType_load[i][0] < 0) continue;
            if (abdLoadType_load[i][1] < 0) continue;
            if (abdLoadType_load[i][0] == 0){
                slightlyLoad[i]=Double.toString(abdLoadType_load[i][1]/abdLoadType_load[i][0]);
            }else {
                slightlyLoad[i]=numberFormat.format(abdLoadType_load[i][1]/abdLoadType_load[i][0]);
            }
        }
        return slightlyLoad;
    }

    /**
     * NA: missing
     * NaN: 0/0
     * @return
     */
    public String[] getStronglyLoad() {
        String[] stronglyLoad=new String[WheatLineage.values().length];
        Arrays.fill(stronglyLoad, "NA");
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(5);
        for (int i = 0; i < abdLoadType_load.length; i++) {
            if (abdLoadType_load[i][0] < 0) continue;
            if (abdLoadType_load[i][2] < 0) continue;
            if (abdLoadType_load[i][0] == 0){
                stronglyLoad[i]=Double.toString(abdLoadType_load[i][2]/abdLoadType_load[i][0]);
            }else {
                stronglyLoad[i]=numberFormat.format(abdLoadType_load[i][2]/abdLoadType_load[i][0]);
            }
        }
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
//                loadType_load[j]=(loadInfo[j*3+1]+0.5*loadInfo[j*3+2])/loadInfo[j*3];
                loadType_load[j]=loadInfo[j*3+1]+0.5*loadInfo[j*3+2];
            }
            this.abdLoadType_load[i]=loadType_load;
        }
    }

    public String toString(){
        StringBuilder sb=new StringBuilder();
        List<String>[] blockGeneName=this.getTriadsBlocks().getBlockGeneName();
        int[] abdGenotypedGeneNum=this.getGenotypedGeneList();
        String[] slightlyLoad=this.getSlightlyLoad();
        String[] stronglyLoad=this.getStronglyLoad();
        sb.append(this.getTriadsBlocks().getTriadsID()).append("\t");
        for (int i = 0; i < blockGeneName.length; i++) {
            sb.append(blockGeneName[i].size()).append(",");
        }
        sb.deleteCharAt(sb.length()-1).append("\t");
        for (int i = 0; i < abdGenotypedGeneNum.length; i++) {
            sb.append(abdGenotypedGeneNum[i]).append(",");
        }
        sb.deleteCharAt(sb.length()-1).append("\t");
        sb.append(String.join(",", slightlyLoad)).append("|");
        sb.append(String.join(",", stronglyLoad));
        return sb.toString();
    }

    /**
     *
     * @param individualTriadsBlockFile
     * @return dim 1 2: additive or dominance, slightly or strongly
     */
    public static List<String>[][] getAdditiveDominanceSlightlyStronglyLoad(File individualTriadsBlockFile){
        List<String>[][] res=new List[2][];
        for (int i = 0; i < res.length; i++) {
            res[i]=new List[2];
            for (int j = 0; j < res[i].length; j++) {
                res[i][j]=new ArrayList<>(19000);
            }
        }
        try (BufferedReader br = IOTool.getReader(individualTriadsBlockFile)) {
            br.readLine();
            String line, slightlyLoadAdditive = null, slightlyLoadDominance=null, stronglyLoadAdditive=null,
                    stronglyLoadDominance=null;
            List<String> temp, tem, teSlightly, teStrongly;
            double[] slightlyLoadArray;
            double[] stronglyLoadArray;
            NumberFormat numberFormat=NumberFormat.getInstance();
            numberFormat.setMaximumFractionDigits(5);
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                tem=PStringUtils.fastSplit(temp.get(3), "|");
                teSlightly=PStringUtils.fastSplit(tem.get(0), ",");
                teStrongly=PStringUtils.fastSplit(tem.get(1), ",");
                slightlyLoadArray=new double[3];
                stronglyLoadArray=new double[3];
                Arrays.fill(slightlyLoadArray, -1);
                Arrays.fill(stronglyLoadArray, -1);
                for (int i = 0; i < teSlightly.size(); i++) {
                    if (!StringTool.isNumeric(teSlightly.get(i))) continue;
                    slightlyLoadArray[i]=Double.parseDouble(teSlightly.get(i));
                }
                for (int i = 0; i < teStrongly.size(); i++) {
                    if (!StringTool.isNumeric(teStrongly.get(i))) continue;
                    stronglyLoadArray[i]=Double.parseDouble(teStrongly.get(i));
                }
                for (int i = 0; i < slightlyLoadArray.length; i++) {
                    if (slightlyLoadArray[i] < 0){
                        slightlyLoadAdditive="NA";
                        slightlyLoadDominance="NA";
                        break;
                    }else {
                        slightlyLoadAdditive=String.valueOf(numberFormat.format(Arrays.stream(slightlyLoadArray).sum()));
                        slightlyLoadDominance= String.valueOf(numberFormat.format(Arrays.stream(slightlyLoadArray).min().getAsDouble()));
                    }
                }
                for (int i = 0; i < stronglyLoadArray.length; i++) {
                    if (stronglyLoadArray[i] < 0){
                        stronglyLoadAdditive="NA";
                        stronglyLoadDominance="NA";
                        break;
                    }else {
                        stronglyLoadAdditive=String.valueOf(numberFormat.format(Arrays.stream(stronglyLoadArray).sum()));
                        stronglyLoadDominance= String.valueOf(numberFormat.format(Arrays.stream(stronglyLoadArray).min().getAsDouble()));
                    }
                }
                res[0][0].add(slightlyLoadAdditive);
                res[0][1].add(stronglyLoadAdditive);
                res[1][0].add(slightlyLoadDominance);
                res[1][1].add(stronglyLoadDominance);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    /**
     *
     * @param individualTriadsBlockFile
     * @return dim 1 2: slightly or strongly
     */
    public static List<String>[] getSlightlyStronglyLoad(File individualTriadsBlockFile){
        List<String>[] res=new List[2];
        for (int i = 0; i < res.length; i++) {
            res[i]=new ArrayList<>(19000);
        }
        res[0]=RowTableTool.getColumnList(individualTriadsBlockFile.getAbsolutePath(), 4);
        res[1]=RowTableTool.getColumnList(individualTriadsBlockFile.getAbsolutePath(), 5);
        return res;
    }
}
