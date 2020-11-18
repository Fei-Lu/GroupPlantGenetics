package daxing.individualIntrogression;

import daxing.common.*;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class PopulationIndividualFd {

    ChrRange[] chrRangeArray;
    double[][][][] dFdValueArray; // dimension 1,2,3,4: chrRange,p2,p3,d fd

    static List<IndividualFdRecord[]> individualFdRecords;
    static List<String> taxaNameList;

    static final double threshForIndividualHasNoIntrogression=0.5;

    private static class IndividualFdRecord{

        double dValue;
        double fdValue;
        String miniIBSP3; // AT001, different from P3 Ae

        public static final IndividualFdRecord DEFAULT=getDefault();

        /**
         * Double.MIN_VALUE表示该个体没有对应的ChrRange
         * @return
         */
        private static final IndividualFdRecord getDefault(){
            return new IndividualFdRecord(Double.MIN_VALUE, Double.MIN_VALUE, ".");
        }

        IndividualFdRecord(double dValue, double fdValue, String miniIBSP3){
            this.dValue=dValue;
            this.fdValue=fdValue;
            this.miniIBSP3=miniIBSP3;
        }

        public boolean contain(P3 p3){
            return miniIBSP3.substring(0,1).equals(p3.name().substring(0,1));
        }
    }

    private static class  SingleWindow2{
        ChrRange chrRange;
        double[][] p2P3_popFd;
        List<String>[][] p2P3_IntrogressedTaxonList;
        List<String>[] p2_nonIntrogressedTaxonList;
        double[][][] p2P3LoadType_introgressedDerivedTotalCount;
        double[][] p2LoadType_nonIntrogressedDerivedTotalCount;

        SingleWindow2(ChrRange chrRange, double[][] p2P3_popFd, List<String>[][] p2P3_IntrogressedTaxonList, List<String>[] p2_nonIntrogressedTaxonList){
            this.chrRange=chrRange;
            this.p2P3_popFd=p2P3_popFd;
            this.p2P3_IntrogressedTaxonList=p2P3_IntrogressedTaxonList;
            this.p2_nonIntrogressedTaxonList=p2_nonIntrogressedTaxonList;

            /**
             * 初始化并赋值为-1
             */
            this.p2P3LoadType_introgressedDerivedTotalCount= new double[P2.values().length][P3.values().length][LoadType.values().length];
            this.p2LoadType_nonIntrogressedDerivedTotalCount= new double[P2.values().length][LoadType.values().length];
            for (int i = 0; i < this.p2P3LoadType_introgressedDerivedTotalCount.length; i++) {
                for (int j = 0; j < p2P3LoadType_introgressedDerivedTotalCount[i].length; j++) {
                    Arrays.fill(p2P3LoadType_introgressedDerivedTotalCount[i][j], -1);
                }
            }
            for (int i = 0; i < p2LoadType_nonIntrogressedDerivedTotalCount.length; i++) {
                Arrays.fill(p2LoadType_nonIntrogressedDerivedTotalCount[i], -1);
            }
        }

        private TIntArrayList[][] getP2P3_IntrogressedTaxonIndex(){
            TIntArrayList[][] p2p3_IntrogressedTaxonIndex=new TIntArrayList[P2.values().length][];
            for (int i = 0; i < p2p3_IntrogressedTaxonIndex.length; i++) {
                p2p3_IntrogressedTaxonIndex[i]=new TIntArrayList[P3.values().length];
                for (int j = 0; j < p2p3_IntrogressedTaxonIndex[i].length; j++) {
                    p2p3_IntrogressedTaxonIndex[i][j]=new TIntArrayList();
                }
            }
            int index;
            for (int i = 0; i < p2P3_IntrogressedTaxonList.length; i++) {
                for (int j = 0; j < p2P3_IntrogressedTaxonList[i].length; j++) {
                    for (int k = 0; k < p2P3_IntrogressedTaxonList[i][j].size(); k++) {
                        index=Collections.binarySearch(taxaNameList, p2P3_IntrogressedTaxonList[i][j].get(k));
                        p2p3_IntrogressedTaxonIndex[i][j].add(index);
                    }
                    p2p3_IntrogressedTaxonIndex[i][j].sort();
                }
            }
            return p2p3_IntrogressedTaxonIndex;
        }

        private TIntArrayList[] getP2_NonIntrogressedTaxonIndex(){
            TIntArrayList[] p2_NonIntrogressedTaxonIndex=new TIntArrayList[P2.values().length];
            for (int i = 0; i < p2_NonIntrogressedTaxonIndex.length; i++) {
                p2_NonIntrogressedTaxonIndex[i]=new TIntArrayList();
            }
            int index;
            for (int i = 0; i < p2_nonIntrogressedTaxonList.length; i++) {
                for (int j = 0; j < p2_nonIntrogressedTaxonList[i].size(); j++) {
                    index=Collections.binarySearch(taxaNameList, p2_nonIntrogressedTaxonList[i].get(j));
                    p2_NonIntrogressedTaxonIndex[i].add(index);
                }
                p2_NonIntrogressedTaxonIndex[i].sort();
            }
            return p2_NonIntrogressedTaxonIndex;
        }

        public void allLoad(Map<ChrPos, String> chrPosLinesMap){
            TIntArrayList[][] p2P3_IntrogressedTaxonIndex=this.getP2P3_IntrogressedTaxonIndex();
            TIntArrayList[] p2_nonIntrogressedTaxonIndex=this.getP2_NonIntrogressedTaxonIndex();
            List<String> temp, tem;
            LoadType loadType;
            double[] loadType_TotalDerivedCount;
            double heter, derived;
            /**
             * one chrRange, one P2P3 combination, one LoadType -> all taxa derived count
             */
            for (int i = 0; i < p2P3_IntrogressedTaxonIndex.length; i++) {
                for (int j = 0; j < p2P3_IntrogressedTaxonIndex[i].length; j++) {
                    if (p2P3_IntrogressedTaxonIndex[i][j].size()==0) continue;
                    loadType_TotalDerivedCount=new double[3];
                    for (Map.Entry<ChrPos, String> entry : chrPosLinesMap.entrySet()){
                        temp=PStringUtils.fastSplit(entry.getValue());
                        loadType=LoadType.valueOf(temp.get(2));
                        for (int l = 0; l < p2P3_IntrogressedTaxonIndex[i][j].size(); l++) {
                            int index=p2P3_IntrogressedTaxonIndex[i][j].get(l) + 3;
                            tem=PStringUtils.fastSplit(temp.get(index), ",");
                            if (Double.parseDouble(tem.get(0)) < 0) continue; //
                            heter=Double.parseDouble(tem.get(0));
                            derived=Double.parseDouble(tem.get(1));
                            loadType_TotalDerivedCount[loadType.getIndex()]+=heter*0.5+derived;
                        }
                    }
                    p2P3LoadType_introgressedDerivedTotalCount[i][j]=loadType_TotalDerivedCount;
                }
            }
            for (int i = 0; i < p2_nonIntrogressedTaxonIndex.length; i++) {
                if (p2_nonIntrogressedTaxonIndex[i].size()==0) continue;
                loadType_TotalDerivedCount=new double[3];
                for (Map.Entry<ChrPos, String> entry : chrPosLinesMap.entrySet()){
                    temp=PStringUtils.fastSplit(entry.getValue());
                    loadType=LoadType.valueOf(temp.get(2));
                    for (int l = 0; l < p2_nonIntrogressedTaxonIndex[i].size(); l++) {
                        int index=p2_nonIntrogressedTaxonIndex[i].get(l) + 3;
                        tem=PStringUtils.fastSplit(temp.get(index), ",");
                        if (Double.parseDouble(tem.get(0)) < 0) continue;
                        heter=Double.parseDouble(tem.get(0));
                        derived=Double.parseDouble(tem.get(1));
                        loadType_TotalDerivedCount[loadType.getIndex()]+=(heter*0.5+derived);
                    }
                }
                p2LoadType_nonIntrogressedDerivedTotalCount[i]=loadType_TotalDerivedCount;
            }
        }

        public String toString(){
            StringBuilder sb=new StringBuilder();
            ChrRange chrRange;
            double[][] p2p3_popFd;
            chrRange=this.chrRange;
            p2p3_popFd=this.p2P3_popFd;
            sb.setLength(0);
            for (int i = 0; i < this.p2P3LoadType_introgressedDerivedTotalCount.length; i++) {
                for (int j = 0; j < this.p2P3LoadType_introgressedDerivedTotalCount[i].length; j++) {
                    if (chrRange.getChr().substring(1,2).equals("D")){
                        if (j < 3) continue;
                    }else{
                        if (j==3) continue;
                    }
                    sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
                    sb.append(chrRange.getEnd()-1).append("\t").append(P2.newInstanceFrom(i).name()).append("\t");
                    sb.append(P3.newInstanceFrom(j).name()).append("\t");
                    if (0 <= p2p3_popFd[i][j] && p2p3_popFd[i][j] <= 1){
                        sb.append(p2p3_popFd[i][j]).append("\t");
                    }else {
                        sb.append("NA").append("\t");
                    }
                    sb.append(p2P3_IntrogressedTaxonList[i][j].size()).append("\t");
                    sb.append(p2_nonIntrogressedTaxonList[i].size()).append("\t");
                    for (int k = 0; k < this.p2P3LoadType_introgressedDerivedTotalCount[i][j].length; k++) {
                        if(this.p2P3LoadType_introgressedDerivedTotalCount[i][j][k] < 0){ //-1
                            sb.append("NA").append("\t");
                        }else {
                            sb.append(this.p2P3LoadType_introgressedDerivedTotalCount[i][j][k]).append("\t");
                        }
                    }
                    for (int k = 0; k < this.p2LoadType_nonIntrogressedDerivedTotalCount[i].length; k++) {
                        if(this.p2LoadType_nonIntrogressedDerivedTotalCount[i][k] < 0){
                            sb.append("NA").append("\t");
                        }else {
                            sb.append(this.p2LoadType_nonIntrogressedDerivedTotalCount[i][k]).append("\t");
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                    sb.append("\n");
                }
            }
            /**
             * 一共写出6行，并去掉最后一行的换行符
             */
            sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }

        public static String getHeader(){
            StringBuilder sb=new StringBuilder();
            sb.append("Chr").append("\t").append("Start").append("\t").append("End").append("\t").append("P2");
            sb.append("\t").append("P3").append("\t").append("PopulationFd").append("\t");
            sb.append("NumAccessionIntrogressed").append("\t").append("NumAccessionNonIntrogressed");
            sb.append("\t").append("IntrogressedSynTotalCount").append("\t").append("IntrogressedNonTotalCount");
            sb.append("\t").append("IntrogressedDelTotalCount").append("\t");
            sb.append("NonIntrogressedSynTotalCount").append("\t").append("NonIntrogressedNonTotalCount").append("\t");
            sb.append("NonIntrogressedDelTotalCount");
            return sb.toString();
        }
    }

    /**
     * Double.MIN_VALUE表示群体p3 AB亚基因组D的默认值，或D亚基因组AB的默认值
     * @param popFdFile
     */
    PopulationIndividualFd(String popFdFile, String individualFdDir){
        ChrRange[] chrRanges=extractChrRange(popFdFile);
        double[][][][] dFdValueArray=new double[chrRanges.length][][][];
        for (int i = 0; i < dFdValueArray.length; i++) {
            dFdValueArray[i]=new double[2][][];
            for (int j = 0; j < dFdValueArray[i].length; j++) {
                dFdValueArray[i][j]=new double[4][];
                for (int k = 0; k < dFdValueArray[i][j].length; k++) {
                    dFdValueArray[i][j][k]=new double[2];
                    Arrays.fill(dFdValueArray[i][j][k], Double.MIN_VALUE);
                }
            }
        }
        try (BufferedReader br = IOTool.getReader(popFdFile)) {
            br.readLine();
            String line, chr, dValue, fdValue;
            List<String> temp;
            int start, end, chrRangeIndex;
            double d, fd;
            P2 p2;
            P3 p3;
            ChrRange chrRange;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                chr=temp.get(0);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                dValue=temp.get(8);
                fdValue=temp.get(9);
                d= NearestIBS.ifDChangeTo0(dValue) ? -1 : Double.parseDouble(dValue);
                fd = NearestIBS.ifFdChangeTo0(dValue, fdValue) ? -1 : Double.parseDouble(fdValue);
                p2=P2.valueOf(temp.get(11));
                p3=P3.valueOf(temp.get(12));
                chrRange=new ChrRange(chr, start, end+1);
                chrRangeIndex=Arrays.binarySearch(chrRanges, chrRange);
                dFdValueArray[chrRangeIndex][p2.index][p3.index][0]=d;
                dFdValueArray[chrRangeIndex][p2.index][p3.index][1]=fd;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.chrRangeArray=chrRanges;
        this.dFdValueArray=dFdValueArray;
        PopulationIndividualFd.individualFdRecords=new ArrayList<>();
        PopulationIndividualFd.taxaNameList=new ArrayList<>();
        List<File> individualFdFiles=IOUtils.getVisibleFileListInDir(individualFdDir);
        for (int i = 0; i < individualFdFiles.size(); i++) {
            this.addIndividual(individualFdFiles.get(i));
        }
    }

    private ChrRange[] extractChrRange(String popFdFile){
        Set<ChrRange> chrRangeSet=new HashSet<>();
        try (BufferedReader br = IOTool.getReader(popFdFile)) {
            br.readLine();
            String line, chr;
            List<String> temp;
            ChrRange chrRange;
            int start, end;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line, ",");
                chr=temp.get(0);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                chrRange=new ChrRange(chr, start, end+1);
                chrRangeSet.add(chrRange);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        List<ChrRange> chrRangeList=new ArrayList<>(chrRangeSet);
        Collections.sort(chrRangeList);
        ChrRange[] chrRanges=new ChrRange[chrRangeList.size()];
        for (int i = 0; i < chrRangeList.size(); i++) {
            chrRanges[i]=chrRangeList.get(i);
        }
        return chrRanges;
    }

    private ChrRange[] getChrRangeArray(){
        return this.chrRangeArray;
    }

    private double[][][][] getdFdValueArray() {
        return dFdValueArray;
    }

    private List<IndividualFdRecord[]> getIndividualFdRecords() {
        return individualFdRecords;
    }

    private int getChrRangeNum(){
        return this.chrRangeArray.length;
    }

    private int getTaxonNum(){
        return this.taxaNameList.size();
    }

    private int binarySearch(ChrRange chrRange){
        return Arrays.binarySearch(getChrRangeArray(), chrRange);
    }

    private void addIndividual(File individualFdFile){
        String taxaName=individualFdFile.getName().substring(0, 5);
        IndividualFdRecord[] individualFdRecords= new IndividualFdRecord[this.getChrRangeNum()];
        Arrays.fill(individualFdRecords, IndividualFdRecord.DEFAULT);
        try (BufferedReader br = IOTool.getReader(individualFdFile)) {
            br.readLine();
            String line, chr, miniIBSP3;
            List<String> temp;
            int start, end, rangeIndex;
            double d, fd;
            ChrRange chrRange;
            IndividualFdRecord individualFdRecord;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line, ",");
                chr=temp.get(0);
                start=Integer.parseInt(temp.get(1));
                end=Integer.parseInt(temp.get(2));
                d=StringTool.isNumeric(temp.get(8)) ? Double.parseDouble(temp.get(8)) : Double.NaN;
                fd=StringTool.isNumeric(temp.get(9)) ? Double.parseDouble(temp.get(9)) : Double.NaN;
                miniIBSP3=temp.get(11);
                individualFdRecord=new IndividualFdRecord(d, fd, miniIBSP3);
                chrRange=new ChrRange(chr, start, end+1);
                rangeIndex=binarySearch(chrRange);
                individualFdRecords[rangeIndex]=individualFdRecord;
            }
            PopulationIndividualFd.individualFdRecords.add(individualFdRecords);
            PopulationIndividualFd.taxaNameList.add(taxaName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private SingleWindow2 getSingleWindow2(int chrRangeIndex){
        ChrRange chrRange=this.getChrRangeArray()[chrRangeIndex];
        List<IndividualFdRecord> individualFdRecordList=new ArrayList<>(); //switch individual all ranges to one
        // range from all individual
        for (int i = 0; i < this.getIndividualFdRecords().size(); i++) {
            individualFdRecordList.add(this.getIndividualFdRecords().get(i)[chrRangeIndex]);
        }
        Set<String>[][] p2P3_IntrogressedTaxonSet=new Set[P2.values().length][];
        Set<String>[] p2_NonIntrogressedTaxonSet=new Set[P2.values().length];
        double[][] p2p3_popFdArray=new double[P2.values().length][];
        for (int i = 0; i < p2P3_IntrogressedTaxonSet.length; i++) {
            p2P3_IntrogressedTaxonSet[i]=new Set[P3.values().length];
            p2_NonIntrogressedTaxonSet[i]=new HashSet<>();
            p2p3_popFdArray[i]=new double[P3.values().length];
            Arrays.fill(p2p3_popFdArray[i], Double.MIN_VALUE);
            for (int j = 0; j < p2P3_IntrogressedTaxonSet[i].length; j++) {
                p2P3_IntrogressedTaxonSet[i][j]=new HashSet<>();
            }
        }
        double[][][] p2P3_DFdValueArray=this.getdFdValueArray()[chrRangeIndex];
        double fd;
        P2 p2;
        P3 p3;
        IndividualFdRecord individualFdRecord;
        double individualFd;
        /**
         * fori CL LR; forj WE DE FTT AT; based on indivi fd value
         * if fd(indi) > fdThreshod(0.5) -> introgression
         * p2P3_IntrogressedTaxonSet[i][j] : taxa list: num_1
         * p2_NonIntrogressedTaxonSet[i] : taxa list: num_2, so num1 + num_2 = 299
         * @Author: Aoyue
         */
        for (int i = 0; i < p2P3_DFdValueArray.length; i++) {
            p2=P2.newInstanceFrom(i);
            for (int j = 0; j < p2P3_DFdValueArray[i].length; j++) {
                p3=P3.newInstanceFrom(j);
                fd=p2P3_DFdValueArray[i][j][1];
                p2p3_popFdArray[i][j]=fd;
                for (int k = 0; k < individualFdRecordList.size(); k++) {
                    if (!(taxaNameList.get(k).substring(0,1).equals(p2.name().substring(0,1)))) continue;
                    individualFdRecord=individualFdRecordList.get(k);
                    individualFd=individualFdRecord.fdValue;
                    /**
                     * why set not list? because when we iterator 299 taxa, it will be distributed into CL-WE CL-DE
                     * CL-FTT LR-WE LR-DE LR-FTT or CL-noIntrogressed or LR-onIntrogressed
                     * E.G. CL_005在个体fd小于0.5时，遍历WE DDE FTT 都加入到p2_NonIntrogressedTaxonSet中，故用set
                     */
                    if (individualFd < threshForIndividualHasNoIntrogression){
                        p2_NonIntrogressedTaxonSet[i].add(taxaNameList.get(k));
                    }else{
                        if (!individualFdRecord.contain(p3)) continue;
                        p2P3_IntrogressedTaxonSet[i][j].add(taxaNameList.get(k));
                    }
                }
            }
        }
        List<String>[][] p2P3_IntrogressedTaxonList=new List[P2.values().length][];
        for (int i = 0; i < p2P3_IntrogressedTaxonList.length; i++) {
            p2P3_IntrogressedTaxonList[i]=new List[P3.values().length];
            for (int j = 0; j < p2P3_IntrogressedTaxonList[i].length; j++) {
                p2P3_IntrogressedTaxonList[i][j]=new ArrayList<>();
            }
        }
        List<String>[] p2_NonIntrogressedTaxonList=new List[P2.values().length];
        for (int i = 0; i < p2_NonIntrogressedTaxonList.length; i++) {
            p2_NonIntrogressedTaxonList[i]=new ArrayList<>();
        }
        for (int i = 0; i < p2P3_IntrogressedTaxonList.length; i++) {
            p2_NonIntrogressedTaxonList[i]=new ArrayList<>(p2_NonIntrogressedTaxonSet[i]);
            Collections.sort(p2_NonIntrogressedTaxonList[i]);
            for (int j = 0; j < p2P3_IntrogressedTaxonList[i].length; j++) {
                p2P3_IntrogressedTaxonList[i][j]=new ArrayList<>(p2P3_IntrogressedTaxonSet[i][j]);
                Collections.sort(p2P3_IntrogressedTaxonList[i][j]);
            }
        }
        return new SingleWindow2(chrRange, p2p3_popFdArray, p2P3_IntrogressedTaxonList, p2_NonIntrogressedTaxonList);
    }

    private void writeChrRange(String outFile){
        try (BufferedWriter bw = IOUtils.getTextWriter(outFile)) {
            bw.write("chr\tstart\tend");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            ChrRange[] chrRanges=this.getChrRangeArray();
            for (int i = 0; i < chrRanges.length; i++) {
                sb.setLength(0);
                sb.append(chrRanges[i].getChr()).append("\t").append(chrRanges[i].getStart());
                sb.append("\t").append(chrRanges[i].getEnd()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void addLoadPerWindow(String individualLoadSummaryFile, String outFile){
        try (BufferedReader br = IOTool.getReader(individualLoadSummaryFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            bw.write(SingleWindow2.getHeader());
            bw.newLine();
            br.readLine();
            int windowNum=this.getChrRangeNum();
            SingleWindow2 singleWindow;
            String line, refChr;
            List<String> temp;
            int chrID, posOnChrID, refPos;
            ChrRange chrRange;
            ChrPos chrPos; //daxing's util
            Map<ChrPos, String> tempMap=new HashMap<>(); //
            for (int i = 0; i < windowNum; i++) {
                singleWindow=this.getSingleWindow2(i);
                chrRange=singleWindow.chrRange;
                retainAll(tempMap, chrRange); //
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID=Integer.parseInt(temp.get(0));
                    posOnChrID=Integer.parseInt(temp.get(1));
                    refChr= RefV1Utils.getChromosome(chrID, posOnChrID);
                    refPos=RefV1Utils.getPosOnChromosome(chrID, posOnChrID);
                    chrPos=new ChrPos(refChr, refPos);
                    if (chrRange.contain(refChr, refPos)){
                        tempMap.put(chrPos, line);
                    }else {
                        singleWindow.allLoad(tempMap);
                        bw.write(singleWindow.toString());
                        bw.newLine();
                        tempMap.put(chrPos, line);
                        break;
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Map<ChrPos, String> retainAll(Map<ChrPos, String> chrPosLineMap, ChrRange chrRange){
        Iterator<ChrPos> iterator=chrPosLineMap.keySet().iterator();
        ChrPos chrPosKey;
        while (iterator.hasNext()){
            chrPosKey=iterator.next();
            if (!chrRange.contain(chrPosKey)){
                iterator.remove();
            }
        }
        return chrPosLineMap;
    }



    public static void start(){
        String popFdFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/002_fdResBySubspecies/002_plot_100SNPwindow_50Step/fdBySubspecies.csv";
        String individualFdDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/006_fdResByIndividual/003_fdByIndividual.newMethod";
        String individualLoadSummaryFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftPerSite/007_individualLoadFdSummary/IndividualLoadFdSummary.txt.gz";
        String outFile="/Users/xudaxing/Desktop/fdLoadBySubspecies100SNPwindow_50Step.txt";
//        String popFdFile="/Users/xudaxing/Desktop/fdLoadRelationship/001_fdBysubspecies/fdBySubspecies.subset.csv";
//        String individualFdDir="/Users/xudaxing/Desktop/fdLoadRelationship/002_fdByIndividual";
//        String individualLoadSummaryFile="/Users/xudaxing/Desktop/fdLoadRelationship/003_IndividualLoadFdSummary/IndividualLoadFdSummary.txt";
//        String outFile="//Users/xudaxing/Desktop/fdLoadRelationship/004_outFile/res2.txt";
        writeWindowSize(popFdFile, individualFdDir, individualLoadSummaryFile, outFile);
    }

    /**
     * WeakLoad and StrongLoad特殊值解释 NaN: 0/0 NA: 为初始化的默认值(即没有渗入) Infinity: 分母为0(即NonIntrogression load为0)
     * PopulationFd特殊值解释 NA: D值不在[0,1]之间或fd值不在[0,1]之间
     * @param popFdFile
     * @param individualFdDir
     * @param individualLoadSummaryFile
     * @param outFile
     */
    public static void writeWindowSize(String popFdFile, String individualFdDir,
                                       String individualLoadSummaryFile, String outFile){
        PopulationIndividualFd populationIndividualFd=new PopulationIndividualFd(popFdFile, individualFdDir);
        populationIndividualFd.addLoadPerWindow(individualLoadSummaryFile, outFile);
    }

}
