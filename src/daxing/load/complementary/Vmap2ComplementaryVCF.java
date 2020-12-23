package daxing.load.complementary;

import daxing.common.ChrRange;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.inference.TestUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.DoublePredicate;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class Vmap2ComplementaryVCF {

    static List<TriadsBlockRecord> triadsBlockRecordList;
    static List<String> taxonList;
    static String[] groupBySubcontinent={"LR_America","LR_Africa","LR_EU","LR_WA","LR_CSA","LR_EA","Cultivar"};

    public Vmap2ComplementaryVCF(String vmap2ComplementaryVCF){
        this.taxonList=new ArrayList<>();
        this.initialize(vmap2ComplementaryVCF);
    }

    private void initialize(String vmap2ComplementaryVCF){
        this.triadsBlockRecordList =new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(vmap2ComplementaryVCF)) {
            String line, triadsBlockID;
            List<String> temp, tem, te, t;
            byte[] blockGeneNum;
            byte[] genotypedGeneNum;
            ChrRange[] chrRange;
            int[] cdsLen;
            String chr;
            int start, end;
            TriadsBlockRecord triadsBlockRecord;
            temp=PStringUtils.fastSplit(br.readLine());
            for (int i = 5; i < temp.size(); i++) {
                this.taxonList.add(temp.get(i));
            }
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                blockGeneNum=new byte[WheatLineage.values().length];
                genotypedGeneNum=new byte[WheatLineage.values().length];
                chrRange=new ChrRange[WheatLineage.values().length];
                cdsLen=new int[WheatLineage.values().length];
                triadsBlockID=temp.get(0);
                tem=PStringUtils.fastSplit(temp.get(1), ",");
                for (int i = 0; i < tem.size(); i++) {
                    blockGeneNum[i]=Byte.parseByte(tem.get(i));
                }
                tem=PStringUtils.fastSplit(temp.get(2), ",");
                for (int i = 0; i < tem.size(); i++) {
                    genotypedGeneNum[i]=Byte.parseByte(tem.get(i));
                }
                tem=PStringUtils.fastSplit(temp.get(3), ";");
                for (int i = 0; i < tem.size(); i++) {
                    te=PStringUtils.fastSplit(tem.get(i), ":");
                    chr=te.get(0);
                    t=PStringUtils.fastSplit(te.get(1), ",");
                    start=Integer.parseInt(t.get(0));
                    end=Integer.parseInt(t.get(1))+1;
                    chrRange[i]=new ChrRange(chr, start, end);
                }
                tem=PStringUtils.fastSplit(temp.get(4), ",");
                for (int i = 0; i < tem.size(); i++) {
                    cdsLen[i]=Integer.parseInt(tem.get(i));
                }
                triadsBlockRecord= new TriadsBlockRecord(triadsBlockID, blockGeneNum, genotypedGeneNum, chrRange, cdsLen);
                for (int i = 5; i < temp.size(); i++) {
                    triadsBlockRecord.addLoadINFO(temp.get(i));
                }
                this.triadsBlockRecordList.add(triadsBlockRecord);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static int getTaxonNum(){
        return Vmap2ComplementaryVCF.taxonList.size();
    }

    public static List<String> getHexaploidTaxonList(String pseudoInfoFile){
        List<String> pseudoTaxonNameList=RowTableTool.getColumnList(pseudoInfoFile,0);
        List<String> taxonList=Vmap2ComplementaryVCF.taxonList;
        Predicate<String> pseudoPredict=taxon ->pseudoTaxonNameList.contains(taxon);
        return taxonList.stream().filter(pseudoPredict.negate()).collect(Collectors.toList());
    }

    public static List<String> getPseudoTaxonList(String pseudoInfoFile){
        return RowTableTool.getColumnList(pseudoInfoFile,0);
    }

    public static TIntArrayList getTaxonIndex(List<String> taxonNameList){
        TIntArrayList indexList = new TIntArrayList();
        for (int i = 0; i < taxonList.size(); i++) {
            if (!taxonNameList.contains(taxonList.get(i))) continue;
            indexList.add(i);
        }
        return indexList;
    }

    public static DoublePredicate getNonNAInfNaNPredict(){
        DoublePredicate isNA=value -> value < 0;
        DoublePredicate isInf=Double::isInfinite;
        DoublePredicate isNaN=Double::isNaN;
        return isNA.negate().and(isInf.negate()).and(isNaN.negate());
    }

    /**
     * 获取伪六倍体Index
     * @param pseudoHexaploidInfo
     * @return
     */
    private static TIntArrayList getPseudoIndexList(String pseudoHexaploidInfo){
        List<String> pseudoNameList=Vmap2ComplementaryVCF.getPseudoTaxonList(pseudoHexaploidInfo);
        return Vmap2ComplementaryVCF.getTaxonIndex(pseudoNameList);
    }

    /**
     * 获取六倍体index
     * @param pseudoHexaploidInfo
     * @return
     */
    private static TIntArrayList getHexaploidIndexList(String pseudoHexaploidInfo){
        List<String> hexaploidNameList=Vmap2ComplementaryVCF.getHexaploidTaxonList(pseudoHexaploidInfo);
        return Vmap2ComplementaryVCF.getTaxonIndex(hexaploidNameList);
    }

    /**
     * 获取六倍体index
     * @param taxaInfoDB
     * @return
     */
    private static TIntArrayList[] getGroupBySubcontinentIndexList(String taxaInfoDB){
        String[] groupBySubcontinent=Vmap2ComplementaryVCF.groupBySubcontinent;
        Arrays.sort(groupBySubcontinent);
        Map<String, String> taxonGroupBySubcontinentMap=RowTableTool.getMap(taxaInfoDB, 0, 24);
        List<String> taxonList=Vmap2ComplementaryVCF.taxonList;
        TIntArrayList[] groupBySubcontinentIndexList=new TIntArrayList[groupBySubcontinent.length];
        for (int i = 0; i < groupBySubcontinentIndexList.length; i++) {
            groupBySubcontinentIndexList[i]=new TIntArrayList();
        }
        int index;
        String value, group;
        for (int i = 0; i < groupBySubcontinentIndexList.length; i++) {
            group=groupBySubcontinent[i];
            for (int j = 0; j < taxonList.size(); j++) {
                if (!taxonGroupBySubcontinentMap.containsKey(taxonList.get(j))) continue;
                value=taxonGroupBySubcontinentMap.get(taxonList.get(j));
                index=Arrays.binarySearch(groupBySubcontinent, value);
                if (index < 0) continue;
                if (!groupBySubcontinent[index].equals(group)) continue;
                groupBySubcontinentIndexList[i].add(j);
            }
        }
        return groupBySubcontinentIndexList;
    }

    public void calculateTwoSampleTTestStatics(String tTestStaticsOutFile, String taxaInfoDB,
                                               String pseudoHexaploidInfo, SubgenomeCombination subgenomeCombination,
                                               WheatLineage positionBySub){
        String[] groupBySubcontinent=Vmap2ComplementaryVCF.groupBySubcontinent;
        Arrays.sort(groupBySubcontinent);
        try (BufferedWriter bw = IOTool.getWriter(tTestStaticsOutFile)) {
            TIntArrayList[] groupBySubcontinentIndexList= Vmap2ComplementaryVCF.getGroupBySubcontinentIndexList(taxaInfoDB);
            TIntArrayList pseudhoIndexList=Vmap2ComplementaryVCF.getPseudoIndexList(pseudoHexaploidInfo);
            bw.write("TriadsBlockID\tGroupBySubcontinent\tGenotypedHexaploidTaxaNum" +
                    "\tGenotypedPseudoTaxaNum\tSlightlyOrStrongly\tAdditiveOrDominance\tT_notEqualVariances\tp_notEqualVariances");
            bw.newLine();
            List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
            double[][][] slightlyStronglyAdditiveDominanceTaxonLoad, slightlyStronglyAdditiveDominancePseudoLoad;
            double[] taxonLoadGreaterOrEqualThan0, pseudoLoadGreaterOrEqualThan0;
            double t, p;
            StringBuilder sb=new StringBuilder();
            SlightlyOrStrongly slightlyOrStrongly;
            AdditiveOrDominance additiveOrDominance;
            DoublePredicate lessThan0=num -> num < 0;
            String[] na=new String[6];
            Arrays.fill(na, "NA");
            ChrRange chrRange;
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                for (int j = 0; j < groupBySubcontinentIndexList.length; j++) {
                    slightlyStronglyAdditiveDominanceTaxonLoad=
                            triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxonListLoad(groupBySubcontinentIndexList[j], subgenomeCombination);
                    slightlyStronglyAdditiveDominancePseudoLoad=
                            triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxonListLoad(pseudhoIndexList, subgenomeCombination);
                    for (int l = 0; l < SlightlyOrStrongly.values().length; l++) {
                        slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(l);
                        for (int m = 0; m < AdditiveOrDominance.values().length; m++) {
                            additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(m);
                            chrRange=triadsBlockRecordList.get(i).getChrRange()[positionBySub.getIndex()];
                            sb.setLength(0);
                            sb.append(triadsBlockRecordList.get(i).getTriadsBlockID()).append("\t");
                            sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
                            sb.append(chrRange.getEnd()-1).append("\t");
                            sb.append(groupBySubcontinent[j]).append("\t");
                            if (Arrays.stream(slightlyStronglyAdditiveDominanceTaxonLoad[l][m]).allMatch(lessThan0)){
                                sb.append(String.join("\t", na));
                                bw.write(sb.toString());
                                bw.newLine();
                                continue;
                            }
                            if (Arrays.stream(slightlyStronglyAdditiveDominancePseudoLoad[l][m]).allMatch(lessThan0)){
                                sb.append(String.join("\t", na));
                                bw.write(sb.toString());
                                bw.newLine();
                                continue;
                            }
                            taxonLoadGreaterOrEqualThan0=
                                    Arrays.stream(slightlyStronglyAdditiveDominanceTaxonLoad[l][m]).filter(Vmap2ComplementaryVCF.getNonNAInfNaNPredict()).toArray();
                            pseudoLoadGreaterOrEqualThan0=
                                    Arrays.stream(slightlyStronglyAdditiveDominancePseudoLoad[l][m]).filter(Vmap2ComplementaryVCF.getNonNAInfNaNPredict()).toArray();
                            if (taxonLoadGreaterOrEqualThan0.length < 2 || pseudoLoadGreaterOrEqualThan0.length < 2){
                                sb.append(taxonLoadGreaterOrEqualThan0.length).append("\t");
                                sb.append(pseudoLoadGreaterOrEqualThan0.length).append("\t");
                                sb.append("NA").append("\t").append("NA").append("\t");
                                sb.append("NA").append("\t").append("NA");
                                bw.write(sb.toString());
                                bw.newLine();
                                continue;
                            }
                            sb.append(taxonLoadGreaterOrEqualThan0.length).append("\t");
                            sb.append(pseudoLoadGreaterOrEqualThan0.length).append("\t");
                            sb.append(slightlyOrStrongly.getValue()).append("\t");
                            sb.append(additiveOrDominance.getValue()).append("\t");
                            t=TestUtils.t(taxonLoadGreaterOrEqualThan0, pseudoLoadGreaterOrEqualThan0);
                            p=TestUtils.tTest(taxonLoadGreaterOrEqualThan0, pseudoLoadGreaterOrEqualThan0);
                            sb.append(t).append("\t").append(p);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void calculateOneSampleTOrZScoreMatrix(String pseudohexaploidInfo,
                                              String tOrZScoreOutFile,
                                               SubgenomeCombination subgenomeCombination, WheatLineage positionBySub,
                                               boolean ifZScore){
        List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(3);
        numberFormat.setGroupingUsed(false);
        try (BufferedWriter bw = IOTool.getWriter(tOrZScoreOutFile)) {
            StringBuilder sb=new StringBuilder();
            sb.setLength(0);
            sb.append("TriadsBlockID\tChr\tPosStart\tPosEnd\tSlightlyOrStrongly\tAdditiveOrDominance\t");
            sb.append(String.join("\t", this.getHexaploidTaxonList(pseudohexaploidInfo)));
            bw.write(sb.toString());
            bw.newLine();
            ChrRange chrRange;
            double[][][] slightlyStronglyAdditiveDominance_OneSampleTOrZScore;
            double tOrZScore;
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                slightlyStronglyAdditiveDominance_OneSampleTOrZScore=
                        ifZScore ?
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo, subgenomeCombination, true) :
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo, subgenomeCombination, false);
                for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                    for (int k = 0; k < AdditiveOrDominance.values().length; k++) {
                        chrRange=triadsBlockRecordList.get(i).getChrRange()[positionBySub.getIndex()];
                        sb.setLength(0);
                        sb.append(triadsBlockRecordList.get(i).getTriadsBlockID()).append("\t");
                        sb.append(chrRange.getChr()).append("\t");
                        sb.append(chrRange.getStart()).append("\t");
                        sb.append(chrRange.getEnd()-1).append("\t");
                        sb.append(SlightlyOrStrongly.newInstanceFromIndex(j).getValue()).append("\t");
                        sb.append(AdditiveOrDominance.newInstanceFromIndex(k).getValue()).append("\t");
                        for (int l = 0; l < slightlyStronglyAdditiveDominance_OneSampleTOrZScore[j][k].length; l++) {
                            tOrZScore=slightlyStronglyAdditiveDominance_OneSampleTOrZScore[j][k][l];
                            if (Double.MIN_VALUE==tOrZScore){
                                sb.append("NA").append("\t");
                                continue;
                            }
                            if (Double.isInfinite(tOrZScore)|| Double.isNaN(tOrZScore)){
                                sb.append(tOrZScore).append("\t");
                                continue;
                            }
                            sb.append(numberFormat.format(tOrZScore)).append("\t");
                        }
                        sb.deleteCharAt(sb.length()-1);
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

    public void calculateOneSampleTorZScoreByTaxon(String pseudohexaploidInfo, String taxaInfoFile,
                                        String oneSampleTByTaxonOutFile,
                                 SubgenomeCombination subgenomeCombination, boolean ifZScore){
        List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(3);
        numberFormat.setGroupingUsed(false);
        try (BufferedWriter bw = IOTool.getWriter(oneSampleTByTaxonOutFile)) {
            StringBuilder sb=new StringBuilder();
            sb.setLength(0);
            sb.append("Taxon\tGroupBySubcontinent\tSlightlyOrStrongly\tAdditiveOrDominance\tOneSampleT_sum" +
                    "\tN");
            bw.write(sb.toString());
            bw.newLine();
            List<String> hexaploidNameList=Vmap2ComplementaryVCF.getHexaploidTaxonList(pseudohexaploidInfo);
            double[][][] slightlyStronglyAdditiveDominance_OneSampleTOrZScore;
            double tOrZScore;
//            int[][][] slightlyStronglyAdditiveDominanceIndiv_countPositiveInf= new int[SlightlyOrStrongly.values().length][][];
//            int[][][] slightlyStronglyAdditiveDominanceIndiv_countNegativeInf= new int[SlightlyOrStrongly.values().length][][];
            int[][][] slightlyStronglyAdditiveDominanceIndiv_count= new int[SlightlyOrStrongly.values().length][][];
            double[][][] slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT=new double[SlightlyOrStrongly.values().length][][];
            for (int i = 0; i < slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT.length; i++) {
                slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[i]=new double[AdditiveOrDominance.values().length][];
                slightlyStronglyAdditiveDominanceIndiv_count[i]=new int[AdditiveOrDominance.values().length][];
//                slightlyStronglyAdditiveDominanceIndiv_countPositiveInf[i]= new int[AdditiveOrDominance.values().length][];
//                slightlyStronglyAdditiveDominanceIndiv_countNegativeInf[i]= new int[AdditiveOrDominance.values().length][];
                for (int j = 0; j < slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[i].length; j++) {
                    slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[i][j]=new double[hexaploidNameList.size()];
                    slightlyStronglyAdditiveDominanceIndiv_count[i][j]=new int[hexaploidNameList.size()];
//                    slightlyStronglyAdditiveDominanceIndiv_countPositiveInf[i][j]=new int[hexaploidNameList.size()];
//                    slightlyStronglyAdditiveDominanceIndiv_countNegativeInf[i][j]=new int[hexaploidNameList.size()];
                }
            }
//            double[] hexaploidMini=new double[hexaploidNameList.size()];
//            double[] hexaploidMax=new double[hexaploidNameList.size()];
//            double max=0, mini=0;
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                slightlyStronglyAdditiveDominance_OneSampleTOrZScore=
                        ifZScore ?
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo, subgenomeCombination, true) :
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo, subgenomeCombination, false);
                for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                    for (int k = 0; k < AdditiveOrDominance.values().length; k++) {
                        for (int l = 0; l < slightlyStronglyAdditiveDominance_OneSampleTOrZScore[j][k].length; l++) {

                            tOrZScore=slightlyStronglyAdditiveDominance_OneSampleTOrZScore[j][k][l];
                            if (Double.isNaN(tOrZScore) || tOrZScore==Double.MIN_VALUE || Double.isInfinite(tOrZScore)) continue;
//                            max=tOrZScore > max ? tOrZScore : max;
//                            mini= tOrZScore < mini ? tOrZScore : mini;
//                            hexaploidMax[l]=max;
//                            hexaploidMini[l]=mini;
                            slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[j][k][l]+=tOrZScore;
                            slightlyStronglyAdditiveDominanceIndiv_count[j][k][l]++;
                        }
                    }
                }
            }
            Map<String, String> taxaGroupBySubcontinentMap=RowTableTool.getMap(taxaInfoFile, 0, 24);
            SlightlyOrStrongly slightlyOrStrongly;
            AdditiveOrDominance additiveOrDominance;
            String taxonName, groupBySubcontinent;
            int count;
            for (int i = 0; i < slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[0][0].length; i++) {
                taxonName=hexaploidNameList.get(i);
                groupBySubcontinent=taxaGroupBySubcontinentMap.get(taxonName);
                if (groupBySubcontinent.equals("OtherHexaploid") | groupBySubcontinent.equals("NA") | groupBySubcontinent.equals("IndianDwarfWheat")) continue;
                for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                    slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(j);
                    for (int k = 0; k < AdditiveOrDominance.values().length; k++) {
                        additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(k);
                        count=slightlyStronglyAdditiveDominanceIndiv_count[j][k][i];
                        tOrZScore=slightlyStronglyAdditiveDominanceIndiv_cumOneSampleT[j][k][i];
                        sb.setLength(0);
                        sb.append(taxonName).append("\t").append(groupBySubcontinent).append("\t");
                        sb.append(slightlyOrStrongly.getValue()).append("\t");
                        sb.append(additiveOrDominance.getValue()).append("\t");
                        sb.append(numberFormat.format(tOrZScore)).append("\t");
                        sb.append(count);
//                        sb.append(slightlyStronglyAdditiveDominanceIndiv_countPositiveInf[j][k][i]).append("\t");
//                        sb.append(slightlyStronglyAdditiveDominanceIndiv_countNegativeInf[j][k][i]);
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

    /**
     *
     * @param pseudohexaploidInfo
     * @param taxaInfoFile
     * @param oneSampleTOrZscoreOutFile
     * @param subgenomeCombination
     * @param ifZScore
     */
    public void calculateOneSampleTorZScoreByTriads(String pseudohexaploidInfo, String taxaInfoFile,
                                   String oneSampleTOrZscoreOutFile, SubgenomeCombination subgenomeCombination,
                                                    WheatLineage positionBySub, boolean ifZScore){
        List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
        double[][][] slightlyStronglyAdditiveDominanceTaxon_TorZScore;
        DoublePredicate na=value -> value==Double.MIN_VALUE;
        DoublePredicate nan=Double::isNaN;
        DoublePredicate inf=Double::isInfinite;
        DoublePredicate nonNANaNInf=na.negate().and(nan.negate()).and(inf.negate());
        TIntArrayList[] groupBySubcontinentIndex=Vmap2ComplementaryVCF.getGroupBySubcontinentIndexList(taxaInfoFile);
        TIntArrayList indexList;
        SlightlyOrStrongly slightlyOrStrongly;
        AdditiveOrDominance additiveOrDominance;
        TDoubleArrayList tOrZScore_list;
        double[] tOrZScoreArray, tOrZScoreArrayNonNANaNInf;
        double tOrZScore;
        StringBuilder sb=new StringBuilder();
        sb.setLength(0);
        String groupBySubcontinent;
        String[] groupBySubcontinentArray=Vmap2ComplementaryVCF.groupBySubcontinent;
        Arrays.sort(groupBySubcontinentArray);
        ChrRange chrRange;
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(3);
        numberFormat.setGroupingUsed(false);
        try (BufferedWriter bw = IOTool.getWriter(oneSampleTOrZscoreOutFile)) {
            sb.append("TriadsBlock\tChr\tPosStart\tPosEnd\tGroupBySubcontinent\tSlightlyOrStrongly" +
                    "\tAdditiveOrDominance\tOneSampleT_sum" +
                    "\tN");
            bw.write(sb.toString());
            bw.newLine();
            for (TriadsBlockRecord triadsBlockRecord: triadsBlockRecordList){
                slightlyStronglyAdditiveDominanceTaxon_TorZScore=ifZScore ?
                        triadsBlockRecord.getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo,
                                subgenomeCombination, true) :
                        triadsBlockRecord.getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(pseudohexaploidInfo, subgenomeCombination, false);
                for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                    slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(i);
                    for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                        additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(j);
                        for (int k = 0; k < groupBySubcontinentIndex.length; k++) {
                            indexList=groupBySubcontinentIndex[k];
                            groupBySubcontinent=groupBySubcontinentArray[k];
                            tOrZScore_list=new TDoubleArrayList();
                            for (int l = 0; l < slightlyStronglyAdditiveDominanceTaxon_TorZScore[i][j].length; l++) {
                                if (!indexList.contains(l)) continue;
                                tOrZScore_list.add(slightlyStronglyAdditiveDominanceTaxon_TorZScore[i][j][l]);
                            }
                            tOrZScoreArray=tOrZScore_list.toArray();
                            tOrZScoreArrayNonNANaNInf= Arrays.stream(tOrZScoreArray).filter(nonNANaNInf).toArray();
                            if (tOrZScoreArrayNonNANaNInf.length < 2) continue;
                            tOrZScore= Arrays.stream(tOrZScoreArrayNonNANaNInf).sum();
                            chrRange=triadsBlockRecord.getChrRange()[positionBySub.getIndex()];
                            sb.setLength(0);
                            sb.append(triadsBlockRecord.getTriadsBlockID()).append("\t");
                            sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
                            sb.append(chrRange.getEnd()-1).append("\t");
                            sb.append(groupBySubcontinent).append("\t");
                            sb.append(slightlyOrStrongly.getValue()).append("\t");
                            sb.append(additiveOrDominance.getValue()).append("\t");
                            sb.append(numberFormat.format(tOrZScore)).append("\t").append(tOrZScoreArrayNonNANaNInf.length);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }


    public enum SlightlyOrStrongly{
        SLIGHTLY(0, "Slightly"), STRONGLY(1, "Strongly");

        int index;
        String value;

        SlightlyOrStrongly(int index, String value) {
            this.index=index;
            this.value=value;
        }

        public int getIndex() {
            return index;
        }

        public String getValue() {
            return value;
        }

        public static SlightlyOrStrongly newInstanceFromIndex(int index){
            switch (index){
                case 0:
                    return SlightlyOrStrongly.SLIGHTLY;
                case 1:
                    return SlightlyOrStrongly.STRONGLY;
                default:
                    System.out.println("please check your parameter index of SlightlyOrStrongly");
                    System.exit(1);
            }
            return null;
        }
    }

    public enum AdditiveOrDominance{
        ADDITIVE(0, "Additive"), DOMINANCE(1, "Dominance");

        int index;
        String value;
        AdditiveOrDominance(int index, String value) {
            this.index=index;
            this.value=value;
        }

        public String getValue() {
            return value;
        }

        public int getIndex() {
            return index;
        }

        public static AdditiveOrDominance newInstanceFromIndex(int index){
            switch (index){
                case 0:
                    return AdditiveOrDominance.ADDITIVE;
                case 1:
                    return AdditiveOrDominance.DOMINANCE;
                default:
                    System.out.println("please check your parameter index of AdditiveOrDominance");
                    System.exit(1);
            }
            return null;
        }
    }

    public enum SubgenomeCombination{
        AB("01"), AD("02"), BD("12"), ABD("012");

        String index;

        SubgenomeCombination(String index) {
            this.index=index;
        }

        public String getIndex() {
            return index;
        }

        public TIntArrayList getIndexList(){
            TIntArrayList res=new TIntArrayList();
            String indexStr=this.getIndex();
            for (int i = 0; i < indexStr.length(); i++) {
                res.add(Integer.parseInt(indexStr.substring(i,i+1)));
            }
            return res;
        }
    }

    private static class TriadsBlockRecord {

        String triadsBlockID;
        byte[] blockGeneNum;
        byte[] genotypedGeneNum;
        /**
         * 21 genes range
         */
        ChrRange[] chrRange;
        int[] cdsLen;
        List<LoadINFO> loadINFO;

        TriadsBlockRecord(String triadsBlockID, byte[] blockGeneNum, byte[] genotypedGeneNum, ChrRange[] chrRange, int[] cdsLen){
            this.triadsBlockID=triadsBlockID;
            this.blockGeneNum=blockGeneNum;
            this.genotypedGeneNum=genotypedGeneNum;
            this.chrRange=chrRange;
            this.cdsLen=cdsLen;
            this.loadINFO=new ArrayList<>();
        }

        public String getTriadsBlockID() {
            return triadsBlockID;
        }

        public ChrRange[] getChrRange() {
            return chrRange;
        }

        boolean addLoadINFO(String loadINFOStr){
            LoadINFO loadINFO= new LoadINFO(loadINFOStr);
            return this.loadINFO.add(loadINFO);
        }

        public List<LoadINFO> getLoadINFO() {
            return loadINFO;
        }

        public List<LoadINFO> getLoadINFO(TIntArrayList taxonIndexList){
            List<LoadINFO> loadINFOS=this.getLoadINFO();
            List<LoadINFO> res=new ArrayList<>();
            for (int i = 0; i < taxonIndexList.size(); i++) {
                res.add(loadINFOS.get(i));
            }
            return  res;
        }

        /**
         *
         * @return [SlightStrongly][AdditiveDominance] load of specified taxonList
         */
        private List<double[][]> getSlightStronglyAdditiveDominanceLoad(TIntArrayList taxonIndexList,
                                                                        SubgenomeCombination subgenomeCombination){
            List<LoadINFO> loadINFOS=this.getLoadINFO(taxonIndexList);
            List<double[][]> res=new ArrayList<>();
            for (int i = 0; i < loadINFOS.size(); i++) {
                res.add(loadINFOS.get(i).getSlightlyStronglyAdditiveDominance_Load(subgenomeCombination));
            }
            return res;
        }

        /**
         * -1: NA
         * @param taxonIndexList
         * @return
         */
        public double[][][] getSlightStronglyAdditiveDominanceTaxonListLoad(TIntArrayList taxonIndexList,
                                                                            SubgenomeCombination subgenomeCombination){
            List<double[][]> slightStronglyAdditiveDominanceTaxonLoad=
                    this.getSlightStronglyAdditiveDominanceLoad(taxonIndexList, subgenomeCombination);
            double[][][] slightlyStronglyAdditiveDominanceTaxonList_load=
                    new double[SlightlyOrStrongly.values().length][][];
            for (int i = 0; i < slightlyStronglyAdditiveDominanceTaxonList_load.length; i++) {
                slightlyStronglyAdditiveDominanceTaxonList_load[i]=new double[AdditiveOrDominance.values().length][];
                for (int j = 0; j < slightlyStronglyAdditiveDominanceTaxonList_load[i].length; j++) {
                    slightlyStronglyAdditiveDominanceTaxonList_load[i][j]=
                            new double[slightStronglyAdditiveDominanceTaxonLoad.size()];
                    Arrays.fill(slightlyStronglyAdditiveDominanceTaxonList_load[i][j], -1);
                }
            }
            for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                    for (int k = 0; k < slightStronglyAdditiveDominanceTaxonLoad.size(); k++) {
                        if (slightStronglyAdditiveDominanceTaxonLoad.get(k)[i][j] < 0) continue;
                        slightlyStronglyAdditiveDominanceTaxonList_load[i][j][k]=
                                slightStronglyAdditiveDominanceTaxonLoad.get(k)[i][j];
                    }
                }
            }
            return slightlyStronglyAdditiveDominanceTaxonList_load;
        }

        /**
         *
         * @param pseudoInfoFile
         * @param taxonIndexList 指定的六倍体
         * @param subgenomeCombination
         * @return
         */
        public double[][] getSlightStronglyAdditiveDominanceTaxon_TwoSampleT(String pseudoInfoFile, TIntArrayList taxonIndexList,
                                                                             SubgenomeCombination subgenomeCombination){
            TIntArrayList pseudoTaxonIndexList=Vmap2ComplementaryVCF.getPseudoIndexList(pseudoInfoFile);
            double[][] res=new double[SlightlyOrStrongly.values().length][];
            for (int i = 0; i < res.length; i++) {
                res[i]=new double[AdditiveOrDominance.values().length];
                Arrays.fill(res[i], Double.MIN_VALUE);
            }
            double[][][] slightlyStronglyAdditiveDominancePseudoLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(pseudoTaxonIndexList,subgenomeCombination);
            double[][][] slightlyStronglyAdditiveDominanceHexaploidLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(taxonIndexList,subgenomeCombination);
            double[] hexaploid, hexaploidRemovedNAInfNaN;;
            double[] pseudoLoad, pseudoLoadRemovedNAInfNaN;
            for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                    pseudoLoad=slightlyStronglyAdditiveDominancePseudoLoad[i][j];
                    pseudoLoadRemovedNAInfNaN= Arrays.stream(pseudoLoad).filter(Vmap2ComplementaryVCF.getNonNAInfNaNPredict()).toArray();
                    hexaploid=slightlyStronglyAdditiveDominanceHexaploidLoad[i][j];
                    hexaploidRemovedNAInfNaN= Arrays.stream(hexaploid).filter(Vmap2ComplementaryVCF.getNonNAInfNaNPredict()).toArray();
                    if (pseudoLoadRemovedNAInfNaN.length < 2 || hexaploidRemovedNAInfNaN.length < 2) continue;
                    res[i][j]=TestUtils.t(hexaploidRemovedNAInfNaN, pseudoLoadRemovedNAInfNaN);
                }
            }
            return res;
        }

        /**
         * Double.MIN_VALUE: NA
         * @param pseudoInfoFile
         * @param subgenomeCombination
         * @return
         */
        public double[][][] getSlightStronglyAdditiveDominanceTaxon_OneSampleTOrZScore(String pseudoInfoFile,
                                                                               SubgenomeCombination subgenomeCombination, boolean ifZScore){
            List<String> pseudoTaxonList= Vmap2ComplementaryVCF.getPseudoTaxonList(pseudoInfoFile);
            List<String> hexaploidTaxonList=Vmap2ComplementaryVCF.getHexaploidTaxonList(pseudoInfoFile);
            TIntArrayList pseudoTaxonIndexList=Vmap2ComplementaryVCF.getTaxonIndex(pseudoTaxonList);
            TIntArrayList hexaploidTaxonIndexList=Vmap2ComplementaryVCF.getTaxonIndex(hexaploidTaxonList);
            double[][][] res=new double[SlightlyOrStrongly.values().length][][];
            for (int i = 0; i < res.length; i++) {
                res[i]=new double[AdditiveOrDominance.values().length][];
                for (int j = 0; j < res[i].length; j++) {
                    res[i][j]=new double[hexaploidTaxonList.size()];
                    Arrays.fill(res[i][j], Double.MIN_VALUE);
                }
            }
            double mean, sd;
            double[][][] slightlyStronglyAdditiveDominancePseudoLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(pseudoTaxonIndexList,subgenomeCombination);
            double[][][] slightlyStronglyAdditiveDominanceHexaploidLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(hexaploidTaxonIndexList,subgenomeCombination);
            double[] pseudoLoad, pseudoLoadRemovedNAInfNaN;
            for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                    pseudoLoad=slightlyStronglyAdditiveDominancePseudoLoad[i][j];
                    pseudoLoadRemovedNAInfNaN= Arrays.stream(pseudoLoad).filter(Vmap2ComplementaryVCF.getNonNAInfNaNPredict()).toArray();
                    if (pseudoLoadRemovedNAInfNaN.length < 2) continue;
                    for (int k = 0; k < hexaploidTaxonList.size(); k++) {
                        if (slightlyStronglyAdditiveDominanceHexaploidLoad[i][j][k] < 0) continue;
                        if (ifZScore){
                            mean= StatUtils.mean(pseudoLoadRemovedNAInfNaN);
                            sd=StatUtils.populationVariance(pseudoLoadRemovedNAInfNaN, mean);
                            res[i][j][k]=(slightlyStronglyAdditiveDominanceHexaploidLoad[i][j][k]-mean)/sd;
                        }else {
                            res[i][j][k]=TestUtils.t(slightlyStronglyAdditiveDominanceHexaploidLoad[i][j][k], pseudoLoadRemovedNAInfNaN);
                        }
                    }
                }
            }
            return res;
        }
    }

    private static class LoadINFO{

        /**
         * -1: NA
         * dim1: slightly or strongly;
         * dim2: subgenome
         */
        double[][] slightlyStronglySubgenome_Load;

        LoadINFO(String loadInfo){
            List<String> temp=PStringUtils.fastSplit(loadInfo, "|");
            double[][] slightlyStronglyLoad=new double[SlightlyOrStrongly.values().length][];
            for (int i = 0; i < slightlyStronglyLoad.length; i++) {
                slightlyStronglyLoad[i]=new double[WheatLineage.values().length];
                Arrays.fill(slightlyStronglyLoad[i], -1);
            }
            List<String> tempSlightlyLod=PStringUtils.fastSplit(temp.get(0), ",");
            List<String> tempStronglyLoad=PStringUtils.fastSplit(temp.get(1), ",");
            String load;
            for (int i = 0; i < WheatLineage.values().length; i++) {
                load=tempSlightlyLod.get(i);
                slightlyStronglyLoad[0][i]=load.equalsIgnoreCase("NA") ? -1 : Double.parseDouble(load);
                load=tempStronglyLoad.get(i);
                slightlyStronglyLoad[1][i]=load.equalsIgnoreCase("NA") ? -1 : Double.parseDouble(load);
            }
            this.slightlyStronglySubgenome_Load =slightlyStronglyLoad;
        }

        public double[][] getSlightlyStronglySubgenome_Load() {
            return slightlyStronglySubgenome_Load;
        }

        public double[][] getSlightlyStronglyAdditiveDominance_Load(SubgenomeCombination subgenomeCombination){
            double[][] slightlyStronglySubgenome_Load=this.getSlightlyStronglySubgenome_Load();
            double[][] slightStronglyAdditiveDominance_SubgenomeCombiantionLoad=new double[SlightlyOrStrongly.values().length][];
            for (int i = 0; i < slightStronglyAdditiveDominance_SubgenomeCombiantionLoad.length; i++) {
                slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i]=new double[AdditiveOrDominance.values().length];
                Arrays.fill(slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i], -1);
            }
            TIntArrayList indexList=subgenomeCombination.getIndexList();
            double[] subgenome_Load, subgenomeCombiantion_Load;
            double mini, sum;
            DoublePredicate lessThan0 = value -> value < 0;
            for (int i = 0; i < slightlyStronglySubgenome_Load.length; i++) {
                subgenome_Load=slightlyStronglySubgenome_Load[i];
                subgenomeCombiantion_Load=new double[indexList.size()];
                Arrays.fill(subgenomeCombiantion_Load, -1);
                for (int j = 0; j < subgenomeCombiantion_Load.length; j++) {
                    subgenomeCombiantion_Load[j]=subgenome_Load[indexList.get(j)];
                }
                if (Arrays.stream(subgenomeCombiantion_Load).allMatch(lessThan0)) continue;
                if (!Arrays.stream(subgenomeCombiantion_Load).filter(lessThan0.negate()).min().isPresent()) continue;
                mini=Arrays.stream(subgenomeCombiantion_Load).filter(lessThan0.negate()).min().getAsDouble();
                if (Arrays.stream(subgenomeCombiantion_Load).min().getAsDouble() < 0){
                    slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][0]=-1;
                }else {
                    sum=Arrays.stream(subgenomeCombiantion_Load).filter(lessThan0.negate()).sum();
                    slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][0]=sum;
                }
                slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][1]=mini;
            }
            return slightStronglyAdditiveDominance_SubgenomeCombiantionLoad;
        }
    }
}
