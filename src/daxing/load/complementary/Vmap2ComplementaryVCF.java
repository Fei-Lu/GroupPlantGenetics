package daxing.load.complementary;

import daxing.common.ChrRange;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.stat.inference.TestUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.function.DoublePredicate;

public class Vmap2ComplementaryVCF {

    List<TriadsBlockRecord> triadsBlockRecordList;
    List<String> taxonList;

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
                    end=Integer.parseInt(t.get(1));
                    chrRange[i]=new ChrRange(chr, start, end);
                }
                tem=PStringUtils.fastSplit(temp.get(4), ",");
                for (int i = 0; i < tem.size(); i++) {
                    cdsLen[i]=Integer.parseInt(tem.get(i));
                }
                triadsBlockRecord=new TriadsBlockRecord(triadsBlockID, blockGeneNum, genotypedGeneNum, chrRange, cdsLen);
                for (int i = 5; i < temp.size(); i++) {
                    triadsBlockRecord.addLoadINFO(temp.get(i));
                }
                this.triadsBlockRecordList.add(triadsBlockRecord);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public List<String> getTaxonList() {
        return taxonList;
    }

    public List<TriadsBlockRecord> getTriadsBlockRecordList() {
        return triadsBlockRecordList;
    }

    public void calculateTwoSampleTTestStatics(String tTestStaticsOutFile, String taxaInfoDB,
                                               String pseudohexaploidInfo){
        String[] groupBySubcontinent={"LR_America","LR_Africa","LR_EU","LR_WA","LR_CSA","LR_EA","Cultivar"};
        Arrays.sort(groupBySubcontinent);
        try (BufferedWriter bw = IOTool.getWriter(tTestStaticsOutFile)) {
            Map<String, String> taxonGroupBySubcontinentMap=RowTableTool.getMap(taxaInfoDB, 0, 24);
            Map<String, String> taxonGroupPseudoMap=RowTableTool.getMap(pseudohexaploidInfo, 0, 1);
            List<String> taxonList=this.getTaxonList();
            TIntArrayList[] groupBySubcontinentIndexList=new TIntArrayList[groupBySubcontinent.length];
            TIntArrayList pseudhoIndexList=new TIntArrayList();
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
            for (int i = 0; i < taxonList.size(); i++) {
                if (!taxonGroupPseudoMap.containsKey(taxonList.get(i))) continue;
                value=taxonGroupPseudoMap.get(taxonList.get(i));
                if (!value.equals("PseudoHexaploid")) continue;
                pseudhoIndexList.add(i);
            }
            bw.write("TriadsBlockID\tGroupBySubcontinent\tGenotypedHexaploidTaxaNum" +
                    "\tGenotypedPseudoTaxaNum\tSlightlyOrStrongly\tAdditiveOrDominance\tT_notEqualVariances\tp_notEqualVariances");
            bw.newLine();
            List<TriadsBlockRecord> triadsBlockRecordList=this.getTriadsBlockRecordList();
            double[][][] slightlyStronglyAdditiveDominanceTaxonLoad, slightlyStronglyAdditiveDominancePseudoLoad;
            double[] taxonLoadGreaterOrEqualThan0, pseudoLoadGreaterOrEqualThan0;
            double t, p;
            StringBuilder sb=new StringBuilder();
            SlightlyOrStrongly slightlyOrStrongly;
            AdditiveOrDominance additiveOrDominance;
            DoublePredicate lessThan0=num -> num < 0;
            DoublePredicate infinite=num -> Double.isInfinite(num);
            DoublePredicate nan =num -> Double.isNaN(num);
            String[] na=new String[6];
            Arrays.fill(na, "NA");
            for (int i = 0; i < triadsBlockRecordList.size(); i++) {
                for (int j = 0; j < groupBySubcontinentIndexList.length; j++) {
                    slightlyStronglyAdditiveDominanceTaxonLoad=
                            triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxonListLoad(groupBySubcontinentIndexList[j]);
                    slightlyStronglyAdditiveDominancePseudoLoad=
                            triadsBlockRecordList.get(i).getSlightStronglyAdditiveDominanceTaxonListLoad(pseudhoIndexList);
                    for (int l = 0; l < SlightlyOrStrongly.values().length; l++) {
                        slightlyOrStrongly=SlightlyOrStrongly.newInstanceFromIndex(l);
                        for (int m = 0; m < AdditiveOrDominance.values().length; m++) {
                            additiveOrDominance=AdditiveOrDominance.newInstanceFromIndex(m);
                            sb.setLength(0);
                            sb.append(triadsBlockRecordList.get(i).getTriadsBlockID()).append("\t");
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
                                    Arrays.stream(slightlyStronglyAdditiveDominanceTaxonLoad[l][m]).filter(lessThan0.negate()).filter(infinite.negate()).filter(nan.negate()).toArray();
                            pseudoLoadGreaterOrEqualThan0=
                                    Arrays.stream(slightlyStronglyAdditiveDominancePseudoLoad[l][m]).filter(lessThan0.negate()).filter(infinite.negate()).filter(nan.negate()).toArray();
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

    private class TriadsBlockRecord {

        String triadsBlockID;
        byte[] blockGeneNum;
        byte[] genotypedGeneNum;
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

        boolean addLoadINFO(String loadINFOStr){
            LoadINFO loadINFO=new LoadINFO(loadINFOStr);
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
        private List<double[][]> getSlightStronglyAdditiveDominanceLoad(TIntArrayList taxonIndexList){
            List<LoadINFO> loadINFOS=this.getLoadINFO(taxonIndexList);
            List<double[][]> res=new ArrayList<>();
            for (int i = 0; i < loadINFOS.size(); i++) {
                res.add(loadINFOS.get(i).getSlightStronglyAdditiveDominanceLoad());
            }
            return res;
        }

        /**
         * -1: NA
         * @param taxonIndexList
         * @return
         */
        public double[][][] getSlightStronglyAdditiveDominanceTaxonListLoad(TIntArrayList taxonIndexList){
            List<double[][]> slightStronglyAdditiveDominanceTaxonLoad=this.getSlightStronglyAdditiveDominanceLoad(taxonIndexList);
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
    }

    private class LoadINFO{

        /**
         * -1: NA
         * dim1: slightly or strongly;
         * dim2: subgenome
         */
        double[][] slightlyStronglySubgenome_Load;

        LoadINFO(double[][] slightlyStronglySubgenome_Load){
            this.slightlyStronglySubgenome_Load = slightlyStronglySubgenome_Load;
        }

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

        /**
         * -1: all ABD is NA
         *
         * @return
         */
        double[][] getSlightStronglyAdditiveDominanceLoad(){
            double[][] slightlyStronglySubgenome_Load=this.getSlightlyStronglySubgenome_Load();
            double[][] slightStronglyAdditiveDominance_Load=new double[SlightlyOrStrongly.values().length][];
            for (int i = 0; i < slightStronglyAdditiveDominance_Load.length; i++) {
                slightStronglyAdditiveDominance_Load[i]=new double[AdditiveOrDominance.values().length];
                Arrays.fill(slightStronglyAdditiveDominance_Load[i], -1);
            }
            double[] subgenome_Load;
            double mini, sum;
            DoublePredicate lessThan0 = value -> value < 0;
            for (int i = 0; i < slightlyStronglySubgenome_Load.length; i++) {
                subgenome_Load=slightlyStronglySubgenome_Load[i];
                if (Arrays.stream(subgenome_Load).allMatch(lessThan0)) continue;
                mini=Arrays.stream(subgenome_Load).filter(lessThan0.negate()).min().getAsDouble();
                if (Arrays.stream(subgenome_Load).min().getAsDouble() < 0){
                    slightStronglyAdditiveDominance_Load[i][0]=-1;
                    slightStronglyAdditiveDominance_Load[i][1]=mini;
                }else {
                    sum=Arrays.stream(subgenome_Load).filter(lessThan0.negate()).sum();
                    slightStronglyAdditiveDominance_Load[i][0]=sum;
                    slightStronglyAdditiveDominance_Load[i][1]=mini;
                }
            }
            return slightStronglyAdditiveDominance_Load;
        }
    }
}
