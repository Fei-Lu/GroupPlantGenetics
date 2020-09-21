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

    List<IndividualFdRecord[]> individualFdRecords;
    List<String> taxaNameList;

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

    /**
     *  循环每一个window
     *  每个window内部都有 2*3 (p2*p3)种组合(introgression list和nonintrogression list)
     */
    private static class SingleWindow {


        ChrRange chrRange;
        /**
         * dimension 1,2: p2,p3
         */
        double[][] fd;
        /**
         * dimension 1,2,3: p2,p3,Introgression_NonIntrogression
         */
        List<String>[][][] p2P3IfIntrogressionTaxonList;
        /**
         * non load and del load
         * dimension 1,2,3,4: p2,p3,Introgression_NonIntrogression,LoadType
         * default Double.MIN_VALUE
         */
        double[][][] load_P2P3LoadType;

        SingleWindow(ChrRange chrRange, double[][] fd, List<String>[][][] p2P3IfIntrogressionTaxonArray){
            this.chrRange=chrRange;
            this.fd=fd;
            this.p2P3IfIntrogressionTaxonList = p2P3IfIntrogressionTaxonArray;
            double[][][] p2P3IfIntrogressionLoadType=new double[p2P3IfIntrogressionTaxonArray.length][][];
            for (int i = 0; i < p2P3IfIntrogressionLoadType.length; i++) {
                p2P3IfIntrogressionLoadType[i]=new double[p2P3IfIntrogressionTaxonArray[i].length][];
                for (int j = 0; j < p2P3IfIntrogressionLoadType[i].length; j++) {
                    p2P3IfIntrogressionLoadType[i][j]=new double[p2P3IfIntrogressionTaxonArray[i][j].length];
                    p2P3IfIntrogressionLoadType[i][j]=new double[2];
                    Arrays.fill(p2P3IfIntrogressionLoadType[i][j], Double.MIN_VALUE);
                }
            }
            this.load_P2P3LoadType =p2P3IfIntrogressionLoadType;
        }

        public ChrRange getChrRange() {
            return chrRange;
        }

        public double[][] getFd() {
            return fd;
        }

        /**
         * 返回SingleWindow.p2P3IfIntrogressionTaxonList在 taxaNameList 中对应的index
         * @param headerList
         * @return
         */
        private TIntArrayList[][][] getP2P3IfIntrogressionTaxonIndexFrom(List<String> headerList){
            List<String>[][][] p2p3IfIntrogressionTaxonList=this.p2P3IfIntrogressionTaxonList;
            TIntArrayList[][][] p2p3IfIntrogressionTaxonIndexList=new TIntArrayList[p2p3IfIntrogressionTaxonList.length][][];
            for (int i = 0; i < p2p3IfIntrogressionTaxonIndexList.length; i++) {
                p2p3IfIntrogressionTaxonIndexList[i]=new TIntArrayList[p2p3IfIntrogressionTaxonList[i].length][];
                for (int j = 0; j < p2p3IfIntrogressionTaxonIndexList[i].length; j++) {
                    p2p3IfIntrogressionTaxonIndexList[i][j]=new TIntArrayList[p2p3IfIntrogressionTaxonList[i][j].length];
                    for (int k = 0; k < p2p3IfIntrogressionTaxonIndexList[i][j].length; k++) {
                        p2p3IfIntrogressionTaxonIndexList[i][j][k]=new TIntArrayList();
                    }
                }
            }
            int taxonIndex;
            for (int i = 0; i < p2p3IfIntrogressionTaxonList.length; i++) {
                for (int j = 0; j < p2p3IfIntrogressionTaxonList[i].length; j++) {
                    for (int k = 0; k < p2p3IfIntrogressionTaxonList[i][j].length; k++) {
                        for (int l = 0; l < p2p3IfIntrogressionTaxonList[i][j][k].size(); l++) {
                            taxonIndex=Collections.binarySearch(headerList,p2p3IfIntrogressionTaxonList[i][j][k].get(l));
                            p2p3IfIntrogressionTaxonIndexList[i][j][k].add(taxonIndex);
                        }
                    }
                }
            }
            return p2p3IfIntrogressionTaxonIndexList;
        }

        public static String getHeaderAB(){
            StringBuilder sb=new StringBuilder();
            sb.append("Chr").append("\t").append("Start").append("\t").append("End").append("\t");
            String[] weakStrong={"WeakDel", "StrongDel"};
            String[] p2={"CL","LR"};
            String[] p3={"WE","DE","FT"};
            for (int i = 0; i < p2.length; i++) {
                for (int j = 0; j < p3.length; j++) {
                    for (int k = 0; k < weakStrong.length; k++) {
                        sb.append(p2[i]).append("_").append(p3[j]).append("_").append(weakStrong[k]).append("\t");
                        if (k>0){
                            sb.append(p2[i]).append("_").append(p3[j]).append("_").append("fdInPop").append("\t");
                        }
                    }
                }

            }
            sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }

        public static String getHeaderD(){
            StringBuilder sb=new StringBuilder();
            sb.append("Chr").append("\t").append("Start").append("\t").append("End").append("\t");
            String[] weakStrong={"WeakDel", "StrongDel"};
            String[] p2={"CL","LR"};
            String[] p3={"AT"};
            for (int i = 0; i < p2.length; i++) {
                for (int j = 0; j < p3.length; j++) {
                    for (int k = 0; k < weakStrong.length; k++) {
                        sb.append(p2[i]).append("_").append(p3[j]).append("_").append(weakStrong[k]).append("\t");
                        if (k>0){
                            sb.append(p2[i]).append("_").append(p3[j]).append("_").append("fdInPop").append("\t");
                        }
                    }
                }

            }
            sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }

        public void allLoad(Map<ChrPos, String> chrPosLineMap, List<String> taxaNameList){
            TIntArrayList[][][] p2P3IfIntrogressionTaxonIndex=getP2P3IfIntrogressionTaxonIndexFrom(taxaNameList);
            List<String> temp, tem;
            LoadType loadType;
            double[] synNonDelTotalDerivedCount;
            double heter, derived;
            double[] nonLoad_IntrogressionNonIntrogression, delLoad_Introgression_NonIntrogression;
            for (int i = 0; i < p2P3IfIntrogressionTaxonIndex.length; i++) {
                for (int j = 0; j < p2P3IfIntrogressionTaxonIndex[i].length; j++) {
                    if (p2P3IfIntrogressionTaxonIndex[i][j][0].size()==0) continue;
                    if (p2P3IfIntrogressionTaxonIndex[i][j][1].size()==0) continue;
                    nonLoad_IntrogressionNonIntrogression=new double[2];
                    delLoad_Introgression_NonIntrogression=new double[2];
                    for (int k = 0; k < p2P3IfIntrogressionTaxonIndex[i][j].length; k++) {
                        synNonDelTotalDerivedCount=new double[3];
                        /**
                         * default load is zero
                         */
                        Arrays.fill(synNonDelTotalDerivedCount, 0);
                        for (Map.Entry<ChrPos, String> entry: chrPosLineMap.entrySet()){
                            temp=PStringUtils.fastSplit(entry.getValue());
                            loadType=LoadType.valueOf(temp.get(2));
                            for (int l = 0; l < p2P3IfIntrogressionTaxonIndex[i][j][k].size(); l++) {
                                int index=p2P3IfIntrogressionTaxonIndex[i][j][k].get(l);
                                tem=PStringUtils.fastSplit(temp.get(index), ",");
                                if (Double.parseDouble(tem.get(0)) < 0) continue;
                                heter=Double.parseDouble(tem.get(0));
                                derived=Double.parseDouble(tem.get(1));
                                synNonDelTotalDerivedCount[loadType.getIndex()]+=heter*0.5+derived;
                            }
                        }
                        nonLoad_IntrogressionNonIntrogression[k]=
                                (synNonDelTotalDerivedCount[1]+synNonDelTotalDerivedCount[2])/synNonDelTotalDerivedCount[0];
                        delLoad_Introgression_NonIntrogression[k]=
                                synNonDelTotalDerivedCount[2]/synNonDelTotalDerivedCount[0];
                    }
                    this.load_P2P3LoadType[i][j][0]=
                            nonLoad_IntrogressionNonIntrogression[0]/nonLoad_IntrogressionNonIntrogression[1];
                    this.load_P2P3LoadType[i][j][1]=
                            delLoad_Introgression_NonIntrogression[0]/delLoad_Introgression_NonIntrogression[1];
                }
            }
        }

        @Override
        public String toString() {
            StringBuilder sb=new StringBuilder();
            ChrRange chrRange;
            double[][] p2p3FdArray;
            chrRange=this.chrRange;
            p2p3FdArray=this.getFd();
            sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
            sb.append(chrRange.getEnd()).append("\t");
            for (int i = 0; i < this.load_P2P3LoadType.length; i++) {
                for (int j = 0; j < this.load_P2P3LoadType[i].length; j++) {
                    if (chrRange.getChr().substring(1,2).equals("D")){
                        if (j < 3) continue;
                    }else{
                        if (j==3) continue;
                    }
                    for (int k = 0; k < this.load_P2P3LoadType[i][j].length; k++) {
                        if(this.load_P2P3LoadType[i][j][k]==Double.MIN_VALUE){
                            sb.append("NA").append("\t");
                        }else if (Double.isFinite(this.load_P2P3LoadType[i][j][k])){
                            sb.append(NumberTool.format(this.load_P2P3LoadType[i][j][k], 5)).append("\t");
                        }else {
                            sb.append(this.load_P2P3LoadType[i][j][k]).append("\t");
                        }
                    }
                    if (0 <= p2p3FdArray[i][j] && p2p3FdArray[i][j] <= 1){
                        sb.append(p2p3FdArray[i][j]).append("\t");
                    }else {
                        sb.append("NA").append("\t");
                    }
                }
            }
            sb.deleteCharAt(sb.length()-1);
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
            String line, chr;
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
                d = StringTool.isNumeric(temp.get(8)) ? Double.parseDouble(temp.get(8)) : Double.NaN;
                fd= StringTool.isNumeric(temp.get(9)) ? Double.parseDouble(temp.get(9)) : Double.NaN;
                p2=P2.valueOf(temp.get(11));
                p3=P3.valueOf(temp.get(12));
                chrRange=new ChrRange(chr, start, end);
                chrRangeIndex=Arrays.binarySearch(chrRanges, chrRange);
                dFdValueArray[chrRangeIndex][p2.index][p3.index][0]=d;
                dFdValueArray[chrRangeIndex][p2.index][p3.index][1]=fd;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        this.chrRangeArray=chrRanges;
        this.dFdValueArray=dFdValueArray;
        this.individualFdRecords=new ArrayList<>();
        this.taxaNameList=new ArrayList<>();
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
                chrRange=new ChrRange(chr, start, end);
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
                chrRange=new ChrRange(chr, start, end);
                rangeIndex=binarySearch(chrRange);
                individualFdRecords[rangeIndex]=individualFdRecord;
            }
            this.individualFdRecords.add(individualFdRecords);
            this.taxaNameList.add(taxaName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private SingleWindow getSingleWindowDFd(int chrRangeIndex){
        ChrRange chrRange=this.getChrRangeArray()[chrRangeIndex];
        List<IndividualFdRecord> individualFdRecordList=new ArrayList<>();
        for (int i = 0; i < this.getIndividualFdRecords().size(); i++) {
            individualFdRecordList.add(this.getIndividualFdRecords().get(i)[chrRangeIndex]);
        }
        List<String>[][][] p2P3IfIntrogressionTaxonList=new List[2][][];
        double[][] p2p3FdArray=new double[2][];
        for (int i = 0; i < p2P3IfIntrogressionTaxonList.length; i++) {
            p2P3IfIntrogressionTaxonList[i]=new List[4][];
            p2p3FdArray[i]=new double[4];
            Arrays.fill(p2p3FdArray[i], Double.MIN_VALUE);
            for (int j = 0; j < p2P3IfIntrogressionTaxonList[i].length; j++) {
                p2P3IfIntrogressionTaxonList[i][j]=new List[2];
                for (int k = 0; k < p2P3IfIntrogressionTaxonList[i][j].length; k++) {
                    p2P3IfIntrogressionTaxonList[i][j][k]=new ArrayList<>();
                }
            }
        }
        double[][][] p2P3DFdValueArray=this.getdFdValueArray()[chrRangeIndex];
        double d, fd;
        List<String> introgressionTaxonList, nonIntrogressionTaxonList;
        P2 p2;
        P3 p3;
        for (int i = 0; i < p2P3DFdValueArray.length; i++) {
            p2=P2.newInstanceFrom(i);
            for (int j = 0; j < p2P3DFdValueArray[i].length; j++) {
                p3=P3.newInstanceFrom(j);
                d=p2P3DFdValueArray[i][j][0];
                fd=p2P3DFdValueArray[i][j][1];
                p2p3FdArray[i][j]=fd;
                if (Double.isNaN(d)) continue;
                if (Double.isNaN(fd)) continue;
                if (d==Double.MIN_VALUE) continue;
                if (fd==Double.MIN_VALUE) continue;
                if (d <= 0) continue; // include Double.MIN_VALUE
                if (d > 1) continue;
                if (fd <= 0) continue; // include Double.MIN_VALUE
                if (fd > 1) continue;
                introgressionTaxonList=new ArrayList<>();
                nonIntrogressionTaxonList=new ArrayList<>();
                for (int k = 0; k < individualFdRecordList.size(); k++) {
                    if (!(taxaNameList.get(k).substring(0,1).equals(p2.name().substring(0,1)))) continue;
                    if (individualFdRecordList.get(k).contain(p3)){
                        introgressionTaxonList.add(taxaNameList.get(k));
                    }else {
                        nonIntrogressionTaxonList.add(taxaNameList.get(k));
                    }
                }
                p2P3IfIntrogressionTaxonList[i][j][0]=introgressionTaxonList;
                p2P3IfIntrogressionTaxonList[i][j][1]=nonIntrogressionTaxonList;
            }
        }
        return new SingleWindow(chrRange, p2p3FdArray, p2P3IfIntrogressionTaxonList);
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
                sb.append("\t").append(chrRanges[i].getEnd());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public void addLoadPerWindow(String individualLoadSummaryFile, String outFileAB, String outFileD){
        try (BufferedReader br = IOTool.getReader(individualLoadSummaryFile);
             BufferedWriter bwAB =IOTool.getWriter(outFileAB);
             BufferedWriter bwD = IOTool.getWriter(outFileD)) {
            bwAB.write(SingleWindow.getHeaderAB());
            bwAB.newLine();
            bwD.write(SingleWindow.getHeaderD());
            bwD.newLine();
            String header=br.readLine();
            List<String> headerList=PStringUtils.fastSplit(header);
            int windowNum=this.getChrRangeNum();
            SingleWindow singleWindow;
            String line, refChr;
            List<String> temp;
            int chrID, posOnChrID, refPos;
            ChrRange chrRange;
            ChrPos chrPos;
            Map<ChrPos, String> tempMap=new HashMap<>();
            for (int i = 0; i < windowNum; i++) {
                singleWindow=this.getSingleWindowDFd(i);
                chrRange=singleWindow.chrRange;
                retainAll(tempMap, chrRange);
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
                        break;
                    }
                }
                singleWindow.allLoad(tempMap, headerList);
                if (chrRange.getChr().substring(1,2).equals("D")){
                    bwD.write(singleWindow.toString());
                    bwD.newLine();
                }else {
                    bwAB.write(singleWindow.toString());
                    bwAB.newLine();
                }
            }
            bwAB.flush();
            bwD.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private Map<ChrPos, String> retainAll(Map<ChrPos, String> chrPosLineMap, ChrRange chrRange){
        Iterator<ChrPos> iterator=chrPosLineMap.keySet().iterator();
        ChrPos chrPosKey;
        while (iterator.hasNext()){
            chrPosKey=iterator.next();
            if (chrRange.contain(chrPosKey)){
                iterator.remove();
            }
        }
        return chrPosLineMap;
    }



    public static void start(){
        String popFdFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/002_fdResBySubspecies/002_plot_100SNPwindow_50Step/fdBySubspecies.csv.gz";
        String individualFdDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/006_fdResByIndividual/002_fdByIndividual";
        String individualLoadSummaryFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftPerSite/007_individualLoadFdSummary/IndividualLoadFdSummary.txt";
        String outFileAB="/Users/xudaxing/Desktop/resAB.txt";
        String outFileD="/Users/xudaxing/Desktop/resD.txt";
        writeWindowSize(popFdFile, individualFdDir, individualLoadSummaryFile, outFileAB, outFileD);
    }

    public static void writeWindowSize(String popFdFile, String individualFdDir,
                                       String individualLoadSummaryFile, String outFileAB, String outFileD){
        PopulationIndividualFd populationIndividualFd=new PopulationIndividualFd(popFdFile, individualFdDir);
        populationIndividualFd.addLoadPerWindow(individualLoadSummaryFile, outFileAB, outFileD);
//        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
//            bw.write("Chr\tStart\tEnd\tCL_WE\tCL_DE\tCL_FT\tLR_WE\tLR_DE\tLR_FT");
//            bw.newLine();
//            StringBuilder sb=new StringBuilder();
//            StringBuilder size=new StringBuilder();
//            SingleWindow singleWindow;
//            ChrRange chrRange;
//            double[][] p2p3FdArray;
//            for (int i = 0; i < populationIndividualFd.getChrRangeNum(); i++) {
//                singleWindow=populationIndividualFd.getSingleWindowDFd(i);
//                chrRange=singleWindow.chrRange;
//                p2p3FdArray=singleWindow.getFd();
//                sb.setLength(0);
//                sb.append(chrRange.getChr()).append("\t").append(chrRange.getStart()).append("\t");
//                sb.append(chrRange.getEnd()).append("\t");
//                for (int j = 0; j < singleWindow.p2P3IfIntrogressionTaxonList.length; j++) {
//                    for (int k = 0; k < singleWindow.p2P3IfIntrogressionTaxonList[j].length; k++) {
//                        if (chrRange.getChr().substring(1,2).equals("D")){
//                            if (k < 3) continue;
//                        }else{
//                            if (k==3) continue;
//                        }
//                        int size0=singleWindow.p2P3IfIntrogressionTaxonList[j][k][0].size();
//                        int size1=singleWindow.p2P3IfIntrogressionTaxonList[j][k][1].size();
//                        size.setLength(0);
//                        size.append(size0).append("_").append(size1);
//                        sb.append(size.toString()).append(",").append(p2p3FdArray[j][k]).append("\t");
//                    }
//                }
//                sb.deleteCharAt(sb.length()-1);
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
    }

}
