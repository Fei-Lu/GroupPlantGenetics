package daxing.load.complementary;

import daxing.common.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.random.EmpiricalDistribution;
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

public class Vmap2ComplementaryVCF {

    public static List<TriadsBlockRecord> triadsBlockRecordList;
    public static List<String> taxonList;
    public static String[] groupBySubcontinent={"LR_America","LR_Africa","LR_EU","LR_WA","LR_CSA","LR_EA","Cultivar"};

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
            int count=0;
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
                count++;
                if (count%5000==0){
                    System.out.println("Reading in "+count+" triads into memory");
                }
            }
            System.out.println("Total "+count+" triads had been reading into memory");
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
        List<String> res=new ArrayList<>();
        for (int i = 0; i < taxonList.size(); i++) {
            if (pseudoTaxonNameList.contains(taxonList.get(i))) continue;
            res.add(taxonList.get(i));
        }
        return res;
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
    public static TIntArrayList getPseudoIndexList(String pseudoHexaploidInfo){
        List<String> pseudoNameList=Vmap2ComplementaryVCF.getPseudoTaxonList(pseudoHexaploidInfo);
        return Vmap2ComplementaryVCF.getTaxonIndex(pseudoNameList);
    }

    /**
     * 获取六倍体index
     * @param pseudoHexaploidInfo
     * @return
     */
    public static TIntArrayList getHexaploidIndexList(String pseudoHexaploidInfo){
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
        System.out.println(DateTime.getDateTimeOfNow());
        System.out.println("Start written to "+tTestStaticsOutFile);
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
            int count=0;
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
                count++;
                if (count%5000==0){
                    System.out.println(count+" triads had been written to "+tTestStaticsOutFile);
                }
            }
            bw.flush();
            System.out.println("Total "+count+" triads had been written to "+tTestStaticsOutFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Two sample t-value had been written to "+tTestStaticsOutFile);
        System.out.println(DateTime.getDateTimeOfNow());
    }

    public void calculateStaticsValueMatrix(String pseudohexaploidInfo,
                                            String tOrZScoreOutFile,
                                            SubgenomeCombination subgenomeCombination, WheatLineage positionBySub,
                                            Statics statics){
        System.out.println(DateTime.getDateTimeOfNow());
        System.out.println("Start writing to "+tOrZScoreOutFile);
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
            double[][][] slightlyStronglyAdditiveDominance_StaticsValue;
            double staticsValue;
            TIntArrayList hexaploidTaxonIndexList=Vmap2ComplementaryVCF.getHexaploidIndexList(pseudohexaploidInfo);
            TIntArrayList pseudoTaxonIndexList=Vmap2ComplementaryVCF.getPseudoIndexList(pseudohexaploidInfo);
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                slightlyStronglyAdditiveDominance_StaticsValue=
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_StaticsValue(pseudoTaxonIndexList, hexaploidTaxonIndexList, subgenomeCombination, statics);
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
                        for (int l = 0; l < slightlyStronglyAdditiveDominance_StaticsValue[j][k].length; l++) {
                            staticsValue=slightlyStronglyAdditiveDominance_StaticsValue[j][k][l];
                            if (Double.MIN_VALUE==staticsValue){
                                sb.append("NA").append("\t");
                                continue;
                            }
                            if (Double.isInfinite(staticsValue)|| Double.isNaN(staticsValue)){
                                sb.append(staticsValue).append("\t");
                                continue;
                            }
                            sb.append(numberFormat.format(staticsValue)).append("\t");
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
        System.out.println(statics.value+" matrix had been written to "+tOrZScoreOutFile);
        System.out.println(DateTime.getDateTimeOfNow());
    }

    public void calculateStaticsValueByTaxon(String pseudohexaploidInfo, String taxaInfoFile,
                                             String oneSampleTByTaxonOutFile,
                                             SubgenomeCombination subgenomeCombination, Statics statics){
        System.out.println(DateTime.getDateTimeOfNow());
        System.out.println("Start written to "+oneSampleTByTaxonOutFile);
        List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(3);
        numberFormat.setGroupingUsed(false);
        try (BufferedWriter bw = IOTool.getWriter(oneSampleTByTaxonOutFile)) {
            StringBuilder sb=new StringBuilder();
            sb.setLength(0);
            sb.append("Taxon\tGroupBySubcontinent\tSlightlyOrStrongly\tAdditiveOrDominance\t");
            sb.append(statics.value).append("\t").append("N");
            bw.write(sb.toString());
            bw.newLine();
            List<String> hexaploidNameList=Vmap2ComplementaryVCF.getHexaploidTaxonList(pseudohexaploidInfo);
            double[][][] slightlyStronglyAdditiveDominance_StaticsValue;
            double staticsValue;
//            int[][][] slightlyStronglyAdditiveDominanceIndiv_countPositiveInf= new int[SlightlyOrStrongly.values().length][][];
//            int[][][] slightlyStronglyAdditiveDominanceIndiv_countNegativeInf= new int[SlightlyOrStrongly.values().length][][];
            int[][][] slightlyStronglyAdditiveDominanceIndiv_count= new int[SlightlyOrStrongly.values().length][][];
            double[][][] slightlyStronglyAdditiveDominanceIndiv_staticsValue=new double[SlightlyOrStrongly.values().length][][];
            for (int i = 0; i < slightlyStronglyAdditiveDominanceIndiv_staticsValue.length; i++) {
                slightlyStronglyAdditiveDominanceIndiv_staticsValue[i]=new double[AdditiveOrDominance.values().length][];
                slightlyStronglyAdditiveDominanceIndiv_count[i]=new int[AdditiveOrDominance.values().length][];
//                slightlyStronglyAdditiveDominanceIndiv_countPositiveInf[i]= new int[AdditiveOrDominance.values().length][];
//                slightlyStronglyAdditiveDominanceIndiv_countNegativeInf[i]= new int[AdditiveOrDominance.values().length][];
                for (int j = 0; j < slightlyStronglyAdditiveDominanceIndiv_staticsValue[i].length; j++) {
                    slightlyStronglyAdditiveDominanceIndiv_staticsValue[i][j]=new double[hexaploidNameList.size()];
                    slightlyStronglyAdditiveDominanceIndiv_count[i][j]=new int[hexaploidNameList.size()];
//                    slightlyStronglyAdditiveDominanceIndiv_countPositiveInf[i][j]=new int[hexaploidNameList.size()];
//                    slightlyStronglyAdditiveDominanceIndiv_countNegativeInf[i][j]=new int[hexaploidNameList.size()];
                }
            }
//            double[] hexaploidMini=new double[hexaploidNameList.size()];
//            double[] hexaploidMax=new double[hexaploidNameList.size()];
//            double max=0, mini=0;
            TIntArrayList hexaploidTaxonIndexList=Vmap2ComplementaryVCF.getHexaploidIndexList(pseudohexaploidInfo);
            TIntArrayList pseudoTaxonIndexList=Vmap2ComplementaryVCF.getPseudoIndexList(pseudohexaploidInfo);
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                slightlyStronglyAdditiveDominance_StaticsValue=
                                triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxon_StaticsValue(pseudoTaxonIndexList, hexaploidTaxonIndexList, subgenomeCombination, statics);
                for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                    for (int k = 0; k < AdditiveOrDominance.values().length; k++) {
                        for (int l = 0; l < slightlyStronglyAdditiveDominance_StaticsValue[j][k].length; l++) {
                            staticsValue=slightlyStronglyAdditiveDominance_StaticsValue[j][k][l];
                            if (Double.isNaN(staticsValue) || staticsValue==Double.MIN_VALUE || Double.isInfinite(staticsValue)) continue;
//                            max=tOrZScore > max ? tOrZScore : max;
//                            mini= tOrZScore < mini ? tOrZScore : mini;
//                            hexaploidMax[l]=max;
//                            hexaploidMini[l]=mini;
                            slightlyStronglyAdditiveDominanceIndiv_staticsValue[j][k][l]+=staticsValue;
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
            for (int i = 0; i < slightlyStronglyAdditiveDominanceIndiv_staticsValue[0][0].length; i++) {
                taxonName=hexaploidNameList.get(i);
                groupBySubcontinent=taxaGroupBySubcontinentMap.get(taxonName);
                if (groupBySubcontinent.equals("OtherHexaploid") | groupBySubcontinent.equals("NA") | groupBySubcontinent.equals("IndianDwarfWheat")) continue;
                for (int j = 0; j < SlightlyOrStrongly.values().length; j++) {
                    slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(j);
                    for (int k = 0; k < AdditiveOrDominance.values().length; k++) {
                        additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(k);
                        count=slightlyStronglyAdditiveDominanceIndiv_count[j][k][i];
                        staticsValue=slightlyStronglyAdditiveDominanceIndiv_staticsValue[j][k][i];
                        sb.setLength(0);
                        sb.append(taxonName).append("\t").append(groupBySubcontinent).append("\t");
                        sb.append(slightlyOrStrongly.getValue()).append("\t");
                        sb.append(additiveOrDominance.getValue()).append("\t");
                        sb.append(numberFormat.format(staticsValue)).append("\t");
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
        System.out.println(statics.value+" by taxon had been written to "+oneSampleTByTaxonOutFile);
        System.out.println(DateTime.getDateTimeOfNow());
    }

    public void calculateStaticsValueByTriads(String pseudohexaploidInfo, String taxaInfoFile,
                                              String staticsValueByTriadsOutFile, SubgenomeCombination subgenomeCombination,
                                              WheatLineage positionBySub, Statics statics){
        System.out.println(DateTime.getDateTimeOfNow());
        System.out.println("Start written to "+staticsValueByTriadsOutFile);
        List<TriadsBlockRecord> triadsBlockRecordList=Vmap2ComplementaryVCF.triadsBlockRecordList;
        double[][][] slightlyStronglyAdditiveDominanceTaxon_StaticsValue;
        TIntArrayList[] groupBySubcontinentIndex=Vmap2ComplementaryVCF.getGroupBySubcontinentIndexList(taxaInfoFile);
        TIntArrayList indexList;
        SlightlyOrStrongly slightlyOrStrongly;
        AdditiveOrDominance additiveOrDominance;
        TDoubleArrayList staticsValue_list;
        double[] staticsValueArray;
        TDoubleArrayList staticsValueArrayNonNANaNInfList;
        double staticsValue;
        StringBuilder sb=new StringBuilder();
        sb.setLength(0);
        String groupBySubcontinent;
        String[] groupBySubcontinentArray=Vmap2ComplementaryVCF.groupBySubcontinent;
        Arrays.sort(groupBySubcontinentArray);
        ChrRange chrRange;
        NumberFormat numberFormat=NumberFormat.getInstance();
        numberFormat.setMaximumFractionDigits(3);
        numberFormat.setGroupingUsed(false);
        try (BufferedWriter bw = IOTool.getWriter(staticsValueByTriadsOutFile)) {
            sb.append("TriadsBlock\tChr\tPosStart\tPosEnd\tGroupBySubcontinent\tSlightlyOrStrongly" +
                    "\tAdditiveOrDominance\t");
            sb.append(statics.value).append("\t").append("N");
            bw.write(sb.toString());
            bw.newLine();
            int count=0;
            TIntArrayList hexaploidTaxonIndexList=Vmap2ComplementaryVCF.getHexaploidIndexList(pseudohexaploidInfo);
            TIntArrayList pseudoTaxonIndexList=Vmap2ComplementaryVCF.getPseudoIndexList(pseudohexaploidInfo);
            for (TriadsBlockRecord triadsBlockRecord: triadsBlockRecordList){
                slightlyStronglyAdditiveDominanceTaxon_StaticsValue=
                        triadsBlockRecord.getSlightStronglyAdditiveDominanceTaxon_StaticsValue(pseudoTaxonIndexList,hexaploidTaxonIndexList,
                                subgenomeCombination, statics);
                for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                    slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(i);
                    for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                        additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(j);
                        for (int k = 0; k < groupBySubcontinentIndex.length; k++) {
                            indexList=groupBySubcontinentIndex[k];
                            groupBySubcontinent=groupBySubcontinentArray[k];
                            staticsValue_list=new TDoubleArrayList();
                            for (int l = 0; l < slightlyStronglyAdditiveDominanceTaxon_StaticsValue[i][j].length; l++) {
                                if (!indexList.contains(l)) continue;
                                staticsValue_list.add(slightlyStronglyAdditiveDominanceTaxon_StaticsValue[i][j][l]);
                            }
                            staticsValueArray=staticsValue_list.toArray();
                            staticsValueArrayNonNANaNInfList=new TDoubleArrayList();
                            for (int l = 0; l < staticsValueArray.length; l++) {
                                if (staticsValueArray[l]==Double.MIN_VALUE) continue;
                                if (Double.isNaN(staticsValueArray[l])) continue;
                                if (Double.isInfinite(staticsValueArray[l])) continue;
                                staticsValueArrayNonNANaNInfList.add(staticsValueArray[l]);
                            }
                            if (staticsValueArrayNonNANaNInfList.size() < 2) continue;
                            staticsValue= staticsValueArrayNonNANaNInfList.sum();
                            chrRange=triadsBlockRecord.getChrRange()[positionBySub.getIndex()];
                            sb.setLength(0);
                            sb.append(triadsBlockRecord.getTriadsBlockID()).append("\t");
                            sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
                            sb.append(chrRange.getEnd()-1).append("\t");
                            sb.append(groupBySubcontinent).append("\t");
                            sb.append(slightlyOrStrongly.getValue()).append("\t");
                            sb.append(additiveOrDominance.getValue()).append("\t");
                            sb.append(numberFormat.format(staticsValue)).append("\t").append(staticsValueArrayNonNANaNInfList.size());
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
                count++;
                if (count%2000==0){
                    System.out.println(count+" triads had been written to "+staticsValueByTriadsOutFile);
                }
            }
            bw.flush();
            System.out.println("Total "+count+" triads had been written to "+ staticsValueByTriadsOutFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(statics.value+" by triads had been written to "+staticsValueByTriadsOutFile);
        System.out.println(DateTime.getDateTimeOfNow());
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

    public enum Statics{
        T("T-value"), Z("Z-score"), Q("Quantile");

        String value;

        Statics(String value) {
            this.value=value;
        }

        public double getHexaploidStaticsValue(double hexaplidValue, double[] pseudoHexaploid) {
            if (this.value.equals("T-value")){
                return TestUtils.t(hexaplidValue, pseudoHexaploid);
            }else if (this.value.equals("Z-score")){
                double mean=StatUtils.mean(pseudoHexaploid);
                double sd=StatUtils.populationVariance(pseudoHexaploid, mean);
                return (hexaplidValue-mean)/sd;
            }else {
                EmpiricalDistribution empiricalDistribution=new EmpiricalDistribution();
                empiricalDistribution.load(pseudoHexaploid);
                return empiricalDistribution.cumulativeProbability(hexaplidValue);
            }
        }
    }

    public static class TriadsBlockRecord {

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
         * @param taxonIndexList
         * @param subgenomeCombination
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
         *
         * @param taxonIndexList
         * @param subgenomeCombination
         * @return -1: NA
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
         *
         * @param pseudoTaxonIndexList
         * @param hexaploidTaxonIndexList
         * @param subgenomeCombination
         * @param statics
         * @return Double.MIN_VALUE: NA
         */
        public double[][][] getSlightStronglyAdditiveDominanceTaxon_StaticsValue(TIntArrayList pseudoTaxonIndexList
                , TIntArrayList hexaploidTaxonIndexList, SubgenomeCombination subgenomeCombination, Statics statics){
            double[][][] res=new double[SlightlyOrStrongly.values().length][][];
            for (int i = 0; i < res.length; i++) {
                res[i]=new double[AdditiveOrDominance.values().length][];
                for (int j = 0; j < res[i].length; j++) {
                    res[i][j]=new double[hexaploidTaxonIndexList.size()];
                    Arrays.fill(res[i][j], Double.MIN_VALUE);
                }
            }
            double[][][] slightlyStronglyAdditiveDominancePseudoLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(pseudoTaxonIndexList,subgenomeCombination);
            double[][][] slightlyStronglyAdditiveDominanceHexaploidLoad=
                    this.getSlightStronglyAdditiveDominanceTaxonListLoad(hexaploidTaxonIndexList,subgenomeCombination);
            double[] pseudoLoad;
            TDoubleArrayList pseudoLoadRemovedNAInfNaNList;
            for (int i = 0; i < SlightlyOrStrongly.values().length; i++) {
                for (int j = 0; j < AdditiveOrDominance.values().length; j++) {
                    pseudoLoad=slightlyStronglyAdditiveDominancePseudoLoad[i][j];
                    pseudoLoadRemovedNAInfNaNList=new TDoubleArrayList();
                    for (int k = 0; k < pseudoLoad.length; k++) {
                        if (pseudoLoad[k] < 0) continue;
                        if (Double.isNaN(pseudoLoad[k])) continue;
                        if (Double.isInfinite(pseudoLoad[k])) continue;
                        pseudoLoadRemovedNAInfNaNList.add(pseudoLoad[k]);
                    }
                    if (pseudoLoadRemovedNAInfNaNList.size() < 2) continue;
                    for (int k = 0; k < hexaploidTaxonIndexList.size(); k++) {
                        if (slightlyStronglyAdditiveDominanceHexaploidLoad[i][j][k] < 0) continue;
                        res[i][j][k]= statics.getHexaploidStaticsValue(slightlyStronglyAdditiveDominanceHexaploidLoad[i][j][k], pseudoLoadRemovedNAInfNaNList.toArray());
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
            double[] subgenome_Load;
            TDoubleArrayList subgenomeCombiantion_LoadList;
            double mini, sum;
            for (int i = 0; i < slightlyStronglySubgenome_Load.length; i++) {
                subgenome_Load=slightlyStronglySubgenome_Load[i];
                subgenomeCombiantion_LoadList=new TDoubleArrayList();
                for (int j = 0; j < indexList.size(); j++) {
                    if (subgenome_Load[indexList.get(j)] < 0) continue;
                    if (Double.isNaN(subgenome_Load[indexList.get(j)])) continue;
                    subgenomeCombiantion_LoadList.add(subgenome_Load[indexList.get(j)]);
                }
                if (subgenomeCombiantion_LoadList.size() < 1) continue;
                mini=subgenomeCombiantion_LoadList.min();
                if (subgenomeCombiantion_LoadList.size() < indexList.size()){
                    slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][0]=-1;
                }else {
                    sum=subgenomeCombiantion_LoadList.sum();
                    slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][0]=sum;
                }
                slightStronglyAdditiveDominance_SubgenomeCombiantionLoad[i][1]=mini;
            }
            return slightStronglyAdditiveDominance_SubgenomeCombiantionLoad;
        }
    }
}
