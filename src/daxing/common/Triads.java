package daxing.common;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import daxing.load.ancestralSite.Standardization;
import daxing.load.ancestralSite.TriadGenotype;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Triads {

    List<TriadsRecord> triadsRecordList;
    /**
     * 0: sort by triadsID
     * 1: sort by sub gene name
     */
    int sortType=0;
    WheatLineage sortTypeBySubGeneName =WheatLineage.A;

    public Triads(String triadFile){
        this.triadsRecordList=new ArrayList<>(19000);
        TriadsRecord triadsRecord;
        try (BufferedReader br = IOTool.getReader(triadFile)) {
            String line, triadID;
            List<String> temp;
            br.readLine();
            String[] triadGeneNameArray;
            int[] cdsLenABD;
            int ifSyntenic,ifExpressed;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                triadID=(temp.get(0));
                triadGeneNameArray=new String[3];
                triadGeneNameArray[0]=temp.get(1);
                triadGeneNameArray[1]=temp.get(2);
                triadGeneNameArray[2]=temp.get(3);
                ifSyntenic=temp.get(4).equals("syntenic") ? 1 : 0;
                ifExpressed=temp.get(5).equals("TRUE") ? 1 : 0;
                cdsLenABD=new int[3];
                cdsLenABD[0]=StringTool.isNumeric(temp.get(6)) ? Integer.parseInt(temp.get(6)) : -1;
                cdsLenABD[1]=StringTool.isNumeric(temp.get(7)) ? Integer.parseInt(temp.get(7)) : -1;
                cdsLenABD[2]=StringTool.isNumeric(temp.get(8)) ? Integer.parseInt(temp.get(8)) : -1;
                triadsRecord=new TriadsRecord(triadID, triadGeneNameArray, ifSyntenic, ifExpressed, cdsLenABD);
                this.triadsRecordList.add(triadsRecord);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void setSortTypeBySubGeneName(WheatLineage sortTypeBySubGeneName) {
        this.sortTypeBySubGeneName = sortTypeBySubGeneName;
        this.sortType=1;
    }

    private void setSortTypeByTriadsID(){
        this.sortType=0;
    }

    public void sortByTriadsID(){
        this.setSortTypeByTriadsID();
        Collections.sort(this.triadsRecordList);
    }

    public void sortBySubGeneName(WheatLineage subgenome){
        this.setSortTypeBySubGeneName(subgenome);
        Collections.sort(this.triadsRecordList);
    }

    public List<TriadsRecord> getTriadsRecordList() {
        return triadsRecordList;
    }

    public class TriadsRecord implements Comparable<TriadsRecord>{

        String triadID;
        String[] triadGeneNameArray;
        int ifSyntenic; //0 is non-syntenic and 1 is syntenic
        int ifExpressed; // 0 is false and 1 is true
        int[] cdsLen; //-1: NA

        TriadsRecord(String triadID, String[] triadGeneNameArray, int ifSyntenic, int ifExpressed, int[] cdsLen){
            this.triadID=triadID;
            this.triadGeneNameArray=triadGeneNameArray;
            this.ifSyntenic=ifSyntenic;
            this.ifExpressed=ifExpressed;
            this.cdsLen=cdsLen;
        }

        public String getTriadID() {
            return triadID;
        }

        public String[] getTriadGeneNameArray() {
            return triadGeneNameArray;
        }

        public int getIfExpressed() {
            return ifExpressed;
        }

        public int getIfSyntenic() {
            return ifSyntenic;
        }

        public int[] getCdsLen() {
            return cdsLen;
        }

        @Override
        /**
         * sort by gene name rather than triadsID
         */
        public int compareTo(TriadsRecord o) {
            if (sortType==0){
                return this.triadID.compareTo(o.triadID);
            }else {
                return this.triadGeneNameArray[sortTypeBySubGeneName.getIndex()].compareTo(o.triadGeneNameArray[sortTypeBySubGeneName.getIndex()]);
            }
        }
    }

    /**
     *
     * @param triadsRecord
     * @return Negative number means not exist
     */
    private int getTriadsRecordIndex(TriadsRecord triadsRecord){
        return Collections.binarySearch(triadsRecordList, triadsRecord);
    }

    /**
     *
     * @param geneName
     * @return Negative number means not exist
     */
    public int getTriadsRecordIndexByGeneName(String geneName){

        String[] geneNameArray=new String[WheatLineage.values().length];
        int[] cdsLen={-1,-1,-1};
        for (int i = 0; i < geneNameArray.length; i++) {
            if(Triads.getSubgenome(geneName).getIndex()!=i) continue;
            geneNameArray[i]=geneName;
        }
        this.setSortTypeBySubGeneName(Triads.getSubgenome(geneName));
        return getTriadsRecordIndex(new TriadsRecord("", geneNameArray, Integer.MIN_VALUE, Integer.MIN_VALUE, cdsLen));
    }

    public int getTriadsRecordIndexByTriadsID(String triadsID){
        this.setSortTypeByTriadsID();
        return getTriadsRecordIndex(new TriadsRecord(triadsID, null, Integer.MIN_VALUE, Integer.MIN_VALUE,null));
    }

    public static WheatLineage getSubgenome(String geneName){
        List<String> temp=PStringUtils.fastSplit(geneName, "G");
        return WheatLineage.valueOf(temp.get(0).substring(8,9));
    }

    public TriadsRecord getTriadRecord(int triadRecordIndex){
        return this.getTriadsRecordList().get(triadRecordIndex);
    }

    public String getTraidID(String geneName){
        return this.getTriadRecord(this.getTriadsRecordIndexByGeneName(geneName)).getTriadID();
    }

    public String getTraidID(int triadRecordIndex){
        return this.getTriadRecord(triadRecordIndex).triadID;
    }

    public int[] getCDSLen(int triadRecordIndex){
        return this.getTriadRecord(triadRecordIndex).cdsLen;
    }

    public List<String> getAllGenes(){
        List<String> allGenes=new ArrayList<>();
        String[] geneName;
        for (int i = 0; i < this.getTriadRecordNum(); i++) {
            geneName=this.getTriadRecord(i).getTriadGeneNameArray();
            for (int j = 0; j < geneName.length; j++) {
                allGenes.add(geneName[j]);
            }
        }
        return allGenes;
    }

    public int getTriadRecordNum(){
        return this.getTriadsRecordList().size();
    }

    public String[] getTraidGenes(String triadsID){
        return this.getTriadRecord(this.getTriadsRecordIndexByTriadsID(triadsID)).getTriadGeneNameArray();
    }

    /**
     * 利用PGF文件、重组率文件、triad文件得到每个triad对应的pos和重组率信息
     * @param pgfFile wheat_v1.1_Lulab_geneHC.pgf
     * @param triadFile triadGenes1.1_cdsLen_geneHC.txt
     * @param recombinationFile iwgsc_refseqv1.0_recombination_rate.txt
     * @param outFile
     */
    public static void triadPosRecombinationRate(String pgfFile, String triadFile, String recombinationFile, String outFile){
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        Triads triads =new Triads(triadFile);
        RowTableTool<String> recombinationRateTable=new RowTableTool<>(recombinationFile);
        Table<String,String,String> recombiantionRate=recombinationRateTable.getTable(0,1,4);
        List<String> posABD=new ArrayList<>();
        String[] abd;
        int start, end, middle;
        int chr,refStart, refEnd;;
        String refChr;
        int geneIndex, recombinationIndex, index, nearestIndex;
        List<String> posList;
        int[] posArray;
        StringBuilder sb=new StringBuilder();
        Map<String,String> posRecombinationMap;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("TriadID\tChr\tPos\tRecombinationRate");
            bw.newLine();
            String recombination;
            int apos=-1;
            for (int i = 0; i < triads.getTriadRecordNum(); i++) {
                abd= triads.getTriadRecord(i).triadGeneNameArray;
                for (int j = 0; j < 3; j++) {
                    geneIndex=pgf.getGeneIndex(abd[j]);
                    if (geneIndex < 0) break;
                    chr=pgf.getGene(geneIndex).geneRange.chr;
                    start=pgf.getGene(geneIndex).geneRange.start;
                    end=pgf.getGene(geneIndex).geneRange.end;
                    refChr=RefV1Utils.getChromosome(chr, 1);
                    refStart=RefV1Utils.getPosOnChromosome(chr, start);
                    refEnd=RefV1Utils.getPosOnChromosome(chr, end);
                    middle=(refStart+refEnd)/2;
                    if (j==0){
                        apos=middle;
                    }
                    chr=pgf.getGene(geneIndex).geneRange.chr;
                    posRecombinationMap=recombiantionRate.row("chr"+refChr);
                    posList=new ArrayList<>(posRecombinationMap.keySet());
                    posArray=posList.stream().mapToInt(Integer::parseInt).map(s->s+5000000).toArray();
                    Arrays.sort(posArray);
                    Collections.sort(posList, Comparator.comparing(l->Integer.parseInt(l)));
                    recombinationIndex=Arrays.binarySearch(posArray, middle);
                    if (recombinationIndex==-1){
                        nearestIndex=0;
                    }else if (recombinationIndex==(0-posArray.length-1)){
                        nearestIndex=posArray.length-1;
                    }else if (recombinationIndex < -1) {
                        index = -recombinationIndex - 2;
                        nearestIndex = ((posArray[index] - middle) < (posArray[index + 1] - middle)) ? index : index + 1;
                    }else {
                        nearestIndex=recombinationIndex;
                    }
                    sb.setLength(0);
                    sb.append(triads.getTraidID(i)).append("\t");
                    sb.append(RefV1Utils.getChromosome(chr, 1)).append("\t").append(middle).append("\t");
                    recombination=posRecombinationMap.get(posList.get(nearestIndex));
                    sb.append(recombination);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 将triadPos文件按A B D pos重新赋值
     * @param triadPosFile triadPos.txt
     * TriadID	Chr	Pos	RecombinationRate
     * T000001	1A	239435	    1.72
     * T000001	6B	111691197	0.06
     * T000001	7D	74092118	0.32
     * T000002	1A	1161012	    1.72
     * T000002	1B	1205398	    0.3
     * T000002	1D	2078877	    0
     * @param outFile triadDPos.txt
     * TriadID	Chr	Pos	RecombinationRate
     * T000001	1A	74092118	1.72
     * T000001	6B	74092118	0.06
     * T000001	7D	74092118	0.32
     * T000002	1A	2078877	    1.72
     * T000002	1B	2078877	    0.3
     * T000002	1D	2078877	    0
     */
    public static void transformTriadPos(String triadPosFile, String outFile){
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line, subgenomePos;
            List<String> tempA, tempB, tempD;
            String header=br.readLine();
            bw.write(header);
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                tempA=PStringUtils.fastSplit(line);
                tempB=PStringUtils.fastSplit(br.readLine());
                tempD=PStringUtils.fastSplit(br.readLine());
                subgenomePos=tempD.get(2);
                tempA.remove(2);
                tempA.add(2, subgenomePos);
                tempB.remove(2);
                tempB.add(2, subgenomePos);
                sb.setLength(0);
                bw.write(String.join("\t", tempA));
                bw.newLine();
                bw.write(String.join("\t", tempB));
                bw.newLine();
                bw.write(String.join("\t", tempD));
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 向六倍体triad群体均值load文件添加每个triad的重组率信息
     * @param triadPosFile
     * TriadID	Chr	Pos	RecombinationRate
     * T000001	1A	239435	    1.72
     * T000001	6B	111691197	0.06
     * T000001	7D	74092118	0.32
     * T000002	1A	1161012	    1.72
     * T000002	1B	1205398	    0.3
     * T000002	1D	2078877	    0
     * @param popMeanLoadFile
     * TriadID	aabbInHexaploid	aabbdd
     * T000002	0.4630606979695424	0.3475395093062605
     * T000003	0.36036078406169736	0.24024052270779797
     * T000004	0.5448779040404048	0.5938889309764311
     * @param outFile
     * TriadID	aabbPopMeanLoadHexaploid	aabbddPopMeanLoadHexaploid	recombinationRate
     * T000002	0.4630606979695424	0.3475395093062605	0.6733333333333333
     * T000003	0.36036078406169736	0.24024052270779797	0.6733333333333333
     * T000004	0.5448779040404048	0.5938889309764311	0.6733333333333333
     */
    public static void recombinationRatePopMeanLoad(String triadPosFile, String popMeanLoadFile, String outFile){
        try (BufferedReader brTriadPos = IOTool.getReader(triadPosFile);
             BufferedReader brPopMeanLoad = IOTool.getReader(popMeanLoadFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            brTriadPos.readLine();
            bw.write("TriadID\taabbPopMeanLoadHexaploid\taabbddPopMeanLoadHexaploid\trecombinationRate");
            bw.newLine();
            String line, triadID;
            List<String> tempA, tempB, tempD, temp;
            double recombinationRateA, recombinationRateB, recombinationRateD, recombinationRateMean;
            Map<String,Double> triadRecombinationRateMap=new HashMap<>();
            while ((line=brTriadPos.readLine())!=null){
                tempA=PStringUtils.fastSplit(line);
                tempB=PStringUtils.fastSplit(brTriadPos.readLine());
                tempD=PStringUtils.fastSplit(brTriadPos.readLine());
                recombinationRateA=Double.parseDouble(tempA.get(3));
                recombinationRateB=Double.parseDouble(tempB.get(3));
                recombinationRateD=Double.parseDouble(tempD.get(3));
                recombinationRateMean=(recombinationRateA+recombinationRateB+recombinationRateD)/3;
                triadRecombinationRateMap.put(tempA.get(0), recombinationRateMean);
            }
            brPopMeanLoad.readLine();
            StringBuilder sb=new StringBuilder();
            while ((line=brPopMeanLoad.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadID=temp.get(0);
                sb.setLength(0);
                sb.append(line).append("\t").append(triadRecombinationRateMap.get(triadID));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 向triadPos重组率文件添加三元图region信息
     * @param triadPosFile
     * TriadID	Chr	Pos	RecombinationRate
     * T000001	1A	239435	    1.72
     * T000001	6B	111691197	0.06
     * T000001	7D	74092118	0.32
     * T000002	1A	1161012	    1.72
     * T000002	1B	1205398	    0.3
     * T000002	1D	2078877	    0
     * @param triadPosRegionFile
     * TriadID	 A	     B	     D	    Region
     * T000001	1.72	0.06	0.32	M100
     * T000002	1.72	0.3	    0	    M100
     */
    public static void triadRecombinationRateModel(String triadPosFile, String triadPosRegionFile){
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw =IOTool.getWriter(triadPosRegionFile)) {
            br.readLine();
            bw.write("TriadID\tA\tB\tD\tRegion");
            bw.newLine();
            String line, recombinationRateA, recombinationRateB, recombinationRateD, region;
            List<String> temp;
            double[] abdRecombinationRate;
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                recombinationRateA=temp.get(3);
                recombinationRateB=PStringUtils.fastSplit(br.readLine()).get(3);
                recombinationRateD=PStringUtils.fastSplit(br.readLine()).get(3);
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(recombinationRateA).append("\t");
                sb.append(recombinationRateB).append("\t").append(recombinationRateD).append("\t");
                abdRecombinationRate=new double[3];
                abdRecombinationRate[0]=Double.parseDouble(recombinationRateA);
                abdRecombinationRate[1]=Double.parseDouble(recombinationRateB);
                abdRecombinationRate[2]=Double.parseDouble(recombinationRateD);
                region= Standardization.getNearestPointIndex(abdRecombinationRate).getRegion();
                sb.append(region);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 将triadAPos  triadBPos triadDPos 滑窗分别得到 triadAPosSliding.txt  triadBPosSliding.txt triadDPosSliding.txt
     * @param triadPosFile
     * TriadID	Chr	Pos	RecombinationRate
     * T000001	1A	239435	1.72
     * T000001	6B	239435	0.06
     * T000001	7D	239435	0.32
     * T000002	1A	1161012	1.72
     * T000002	1B	1161012	0.3
     * T000002	1D	1161012	0
     * @param windowSize
     * @param stepSize
     * @param outFile
     * Chr	Start	End	        Value
     * 1A	0	    10000000	1.6878048780487804
     * 1A	1000000	11000000	1.6969767441860466
     * 1A	2000000	12000000	1.7482499999999999
     * 1A	3000000	13000000	1.7861702127659576
     * 1A	4000000	14000000	1.7667391304347826
     * @param triadSubgenome A or B or D
     */
    public static void triadPosSlidingWindow(String triadPosFile, int windowSize, int stepSize, String outFile,
                                              String triadSubgenome){
        try (BufferedReader br = IOTool.getReader(triadPosFile)) {
            String line;
            List<String>[] tempABD;
            StringBuilder sb=new StringBuilder();
            String chrNum;
            String header="Chr\tPos\tRecombinationRate";
            List<List<String>> cell=new ArrayList<>();
            br.readLine();
            int[] triadAPos={0,1,2};
            int[] triadBPos={1,0,2};
            int[] triadDPos={2,0,1};
            int[] triadPos=null;
            String[] bdSub={"B","D"};
            String[] adSub={"A","D"};
            String[] abSub={"A","B"};
            String[] subgenome=null;
            if (triadSubgenome.equals("A")){
                triadPos=triadAPos;
                subgenome=bdSub;
            }else if (triadSubgenome.equals("B")){
                triadPos=triadBPos;
                subgenome=adSub;
            }else if (triadSubgenome.equals("D")){
                triadPos=triadDPos;
                subgenome=abSub;
            }
            while ((line=br.readLine())!=null){
                tempABD=new List[3];
                tempABD[0]=PStringUtils.fastSplit(line);
                tempABD[1]=PStringUtils.fastSplit(br.readLine());
                tempABD[2]=PStringUtils.fastSplit(br.readLine());
                sb.setLength(0);
                chrNum=tempABD[triadPos[0]].get(1).substring(0,1);
                tempABD[triadPos[1]].set(1, chrNum+subgenome[0]);
                tempABD[triadPos[2]].set(1, chrNum+subgenome[1]);
                tempABD[0].remove(0);
                tempABD[1].remove(0);
                tempABD[2].remove(0);
                cell.add(tempABD[0]);
                cell.add(tempABD[1]);
                cell.add(tempABD[2]);
            }
            RowTableTool<String> rowTable=new RowTableTool<>(PStringUtils.fastSplit(header), cell);
            String[] tempDirName={"splitChr", "slidingWindow"};
            File[] tempDirs=new File[tempDirName.length];
            for (int i = 0; i < tempDirName.length; i++) {
                tempDirs[i]=new File(new File(outFile).getParent(), tempDirName[i]);
                tempDirs[i].mkdir();
            }
            File tempFile=new File(tempDirs[0], "temp.vcf");
            rowTable.write(tempFile, IOFileFormat.Text);
            VCF vcf=new VCF(tempFile);
            vcf.writeVcfToSplitedChr(tempDirs[0].getAbsolutePath());
            tempFile.delete();
            PlotTools.slidingWindow(tempDirs[0].getAbsolutePath(), tempDirs[1].getAbsolutePath(), 2);
            PlotTools.merge(tempDirs[1].getAbsolutePath(), outFile);

        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param nonsynGlobalIndividualFile nonsynGlobalIndividual.txt
     * @param triadRecombinationRateRegionFile
     * TriadID	A	    B	    D	    Region
     * T000001	1.72	0.06	0.32	M100
     * T000002	1.72	0.3	    0	    M100
     * T000003	1.72	0.3	    0	    M100
     * T000004	1.72	0.3	    0	    M100
     * @param triadLoadRecombinationRegionOutFile
     * TriadID	LoadA	LoadB	LoadD	LoadRegion	RecombinationRateA	RecombinationRateB	RecombinationRateD
     * RecombinationRateRegion
     * T000001	0.06701	NA	NA	NA	1.72	0.06	0.32	M100
     * T000002	0.30384	0.61895	0.11533	M110	1.72	0.3	0	M100
     * T000003	0.72134	0.0	0.0	M100	1.72	0.3	0	M100
     */
    public static void triadGlobalLoadRegionRecombinationRateRegion(String nonsynGlobalIndividualFile,
                                                                    String triadRecombinationRateRegionFile,
                                                                    String triadLoadRecombinationRegionOutFile){
        try (BufferedReader brLoad = IOTool.getReader(nonsynGlobalIndividualFile);
             BufferedReader brRecombination= IOTool.getReader(triadRecombinationRateRegionFile);
             BufferedWriter bw = IOTool.getWriter(triadLoadRecombinationRegionOutFile)) {
            brLoad.readLine();
            brRecombination.readLine();
            bw.write("TriadID\tLoadA\tLoadB\tLoadD\tLoadRegion\tRecombinationRateA\tRecombinationRateB" +
                    "\tRecombinationRateD\tRecombinationRateRegion");
            bw.newLine();
            String line, triadID;
            List<String> temp;
            Multimap<String,String> triadRecombinationMap= ArrayListMultimap.create();
            while ((line=brRecombination.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadRecombinationMap.put(temp.get(0), temp.get(1));
                triadRecombinationMap.put(temp.get(0), temp.get(2));
                triadRecombinationMap.put(temp.get(0), temp.get(3));
                triadRecombinationMap.put(temp.get(0), temp.get(4));
            }
            StringBuilder sb=new StringBuilder();
            List<String> values;
            while ((line=brLoad.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                sb.setLength(0);
                triadID=temp.get(0);
                sb.append(triadID).append("\t").append(temp.get(1)).append("\t").append(temp.get(2)).append("\t");
                sb.append(temp.get(3)).append("\t").append(temp.get(4)).append("\t");
                values=new ArrayList<>(triadRecombinationMap.get(triadID));
                sb.append(String.join("\t", values));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void triadGlobalLoadRegionRecombinationRateRegion(String nonsynGlobalIndividualFile,
                                                                    String triadRecombinationRateRegionFile,
                                                                    String triadLoadRecombinationRegionOutFile,
                                                                    String taxonName){
        try (BufferedReader brLoad = IOTool.getReader(nonsynGlobalIndividualFile);
             BufferedReader brRecombination= IOTool.getReader(triadRecombinationRateRegionFile);
             BufferedWriter bw = IOTool.getWriter(triadLoadRecombinationRegionOutFile)) {
            String header=brLoad.readLine();
            List<String> taxonNames=PStringUtils.fastSplit(header);
            int index=taxonNames.indexOf(taxonName);
            brRecombination.readLine();
            bw.write("TriadID\tLoadA\tLoadB\tLoadD\tLoadRegion\tRecombinationRateA\tRecombinationRateB" +
                    "\tRecombinationRateD\tRecombinationRateRegion");
            bw.newLine();
            String line, triadID;
            List<String> temp;
            Multimap<String,String> triadRecombinationMap= ArrayListMultimap.create();
            while ((line=brRecombination.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                triadRecombinationMap.put(temp.get(0), temp.get(1));
                triadRecombinationMap.put(temp.get(0), temp.get(2));
                triadRecombinationMap.put(temp.get(0), temp.get(3));
                triadRecombinationMap.put(temp.get(0), temp.get(4));
            }
            StringBuilder sb=new StringBuilder();
            List<String> values;
            while ((line=brLoad.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                sb.setLength(0);
                triadID=temp.get(0);
                sb.append(triadID).append("\t").append(temp.get(index)).append("\t").append(temp.get(index+1)).append(
                        "\t");
                sb.append(temp.get(index+2)).append("\t").append(temp.get(index+3)).append("\t");
                values=new ArrayList<>(triadRecombinationMap.get(triadID));
                sb.append(String.join("\t", values));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Individual Additive Model
     * @param loadThreshFile hexaploidLoadThresh3.0.txt
     * @param triadPopMeanLoadOutFile hexaploidTriadPopMeanLoad.txt
     */
    public static void triadPopMeanLoadAdditive(String loadThreshFile, String triadPopMeanLoadOutFile){
        TriadGenotype triadGenotype=new TriadGenotype(loadThreshFile);
        List<String> temp, tem, te;
        TDoubleArrayList aabb;
        TDoubleArrayList aabbdd;
        double loadA, loadB, loadD, sumAB, sumABD;
        try (BufferedWriter bw = IOTool.getWriter(triadPopMeanLoadOutFile)) {
            bw.write("TriadID\taabbInHexaploid\taabbdd");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < triadGenotype.columnTableTool.getRowNumber(); i++) {
                temp=triadGenotype.columnTableTool.getRow(i);
                aabb=new TDoubleArrayList();
                aabbdd= new TDoubleArrayList();
                for (int j = 2; j < temp.size(); j++) {
                    tem= PStringUtils.fastSplit(temp.get(j), ":");
                    te=PStringUtils.fastSplit(tem.get(0), ",");
                    if (te.get(0).equals("NA")) continue;
                    if (te.get(1).equals("NA")) continue;
                    if (te.get(2).equals("NA")) continue;
                    loadA=Double.parseDouble(te.get(0));
                    loadB=Double.parseDouble(te.get(1));
                    loadD=Double.parseDouble(te.get(2));
                    sumAB=loadA+loadB;
                    sumABD=loadA+loadB+loadD;
                    aabb.add(sumAB/2);
                    aabbdd.add(sumABD/3);
                }
                if (aabb.size() < 1 || aabbdd.size() < 1) continue;
                double meanTetraploid=aabb.sum()/aabb.size();
                double meanHexaploid=aabbdd.sum()/aabbdd.size();
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(meanTetraploid).append("\t").append(meanHexaploid);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * Individual Dominant Model
     * @param loadThreshFile hexaploidLoadThresh3.0.txt
     * @param triadPopMeanLoadOutFile hexaploidTriadPopMeanLoadDominant.txt
     */
    public static void triadPopMeanLoadDominant(String loadThreshFile, String triadPopMeanLoadOutFile){
        TriadGenotype triadGenotype=new TriadGenotype(loadThreshFile);
        List<String> temp, tem, te;
        TDoubleArrayList aabb;
        TDoubleArrayList aabbdd;
        double loadA, loadB, loadD, sumAB, sumABD;
        try (BufferedWriter bw = IOTool.getWriter(triadPopMeanLoadOutFile)) {
            bw.write("TriadID\taabbInHexaploid\taabbdd");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < triadGenotype.columnTableTool.getRowNumber(); i++) {
                temp=triadGenotype.columnTableTool.getRow(i);
                aabb=new TDoubleArrayList();
                aabbdd= new TDoubleArrayList();
                for (int j = 2; j < temp.size(); j++) {
                    tem= PStringUtils.fastSplit(temp.get(j), ":");
                    te=PStringUtils.fastSplit(tem.get(0), ",");
                    if (te.get(0).equals("NA")) continue;
                    if (te.get(1).equals("NA")) continue;
                    if (te.get(2).equals("NA")) continue;
                    loadA=Double.parseDouble(te.get(0));
                    loadB=Double.parseDouble(te.get(1));
                    loadD=Double.parseDouble(te.get(2));
                    sumAB=loadA+loadB;
                    sumABD=loadA+loadB+loadD;
                    aabb.add(sumAB/2);
                    aabbdd.add(sumABD/3);
                }
                if (aabb.size() < 1 || aabbdd.size() < 1) continue;
                double meanTetraploid=aabb.sum()/aabb.size();
                double meanHexaploid=aabbdd.sum()/aabbdd.size();
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t").append(meanTetraploid).append("\t").append(meanHexaploid);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
