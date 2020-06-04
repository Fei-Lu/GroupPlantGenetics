package daxing.common;

import com.google.common.collect.Table;
import daxing.load.ancestralSite.Standardization;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class Triad {

    List<String> triadID;
    List<String[]> triad;
    TIntArrayList ifSyntenic;  //0 is non-syntenic and 1 is syntenic
    TIntArrayList ifExpressed;  // 0 is false and 1 is true
    List<int[]> cdsList;

    public Triad(String triadFile){
        triadID=new ArrayList<>(19000);
        triad=new ArrayList<>(19000);
        ifSyntenic=new TIntArrayList();
        ifExpressed=new TIntArrayList();
        cdsList=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(triadFile)) {
            String line;
            List<String> temp;
            br.readLine();
            String[] abdGenes;
            int[] cdsLenABD;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                triadID.add(temp.get(0));
                abdGenes=new String[3];
                abdGenes[0]=temp.get(1);
                abdGenes[1]=temp.get(2);
                abdGenes[2]=temp.get(3);
                triad.add(abdGenes);
                if (temp.get(4).equals("syntenic")){
                    ifSyntenic.add(1);
                }else {
                    ifSyntenic.add(0);
                }
                if (temp.get(5).equals("TRUE")){
                    ifExpressed.add(1);
                }else {
                    ifExpressed.add(0);
                }
                cdsLenABD=new int[3];
                cdsLenABD[0]=Integer.parseInt(temp.get(6));
                cdsLenABD[1]=Integer.parseInt(temp.get(7));
                cdsLenABD[2]=Integer.parseInt(temp.get(8));
                cdsList.add(cdsLenABD);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public List<String[]> getTriad() {
        return triad;
    }

    public TIntArrayList getIfExpressed() {
        return ifExpressed;
    }

    public TIntArrayList getIfSyntenic() {
        return ifSyntenic;
    }

    public List<int[]> getCdsList() {
        return cdsList;
    }

    public List<String> getTriadID() {
        return triadID;
    }

    public WheatLineage getSubgenome(String geneName){
        List<String> temp=PStringUtils.fastSplit(geneName, "G");
        String subgenome=temp.get(0).substring(8,9);
        return WheatLineage.valueOf(subgenome);
    }

    public List<String> getSubgenomeGenes(WheatLineage subgenome){
        String[] subgenomes={"A","B","D"};
        int index= Arrays.binarySearch(subgenomes, subgenome.name());
        if (index < 0) {
            System.out.println("subgenome must be A or B or D");
            System.exit(1);
        }
        List<String[]> genes=this.getTriad();
        List<String> subgenomeGenes=new ArrayList<>();
        for (int i = 0; i < genes.size(); i++) {
            subgenomeGenes.add(genes.get(i)[index]);
        }
        return subgenomeGenes;
    }

    public List<String> getSubgenomeGens(String geneName){
        WheatLineage subgenome=this.getSubgenome(geneName);
        return this.getSubgenomeGenes(subgenome);
    }

    public List<String> getAllGenes(){
        List<String> genes=new ArrayList<>(57000);
        List<String[]> data=this.getTriad();
        for (int i = 0; i < data.size(); i++) {
            for (int j = 0; j < data.get(i).length; j++) {
                genes.add(data.get(i)[j]);
            }
        }
        return genes;
    }

    public String[] getTriad(int geneIndex){
        return this.getTriad().get(geneIndex);
    }

    public int[] getCDSLen(int geneIndex){
        return this.cdsList.get(geneIndex);
    }

    public String[] getTriad(String geneName){
        int geneIndex=this.getGeneIndex(geneName);
        return getTriad(geneIndex);
    }

    public String[] getTraidGenes(String triadID){
        List<String> triadNames=this.triadID;
        int index=triadNames.indexOf(triadID);
        return this.getTriad().get(index);
    }

    public String getTraidID(String geneName){
        int geneIndex=this.getGeneIndex(geneName);
        return this.getTriadID().get(geneIndex);
    }

    public String getTraidID(int index){
        return this.getTriadID().get(index);
    }

    public int getGeneIndex(String geneName){
        String subgenome=this.getSubgenome(geneName).name();
        List<String> subgenomeGenes=this.getSubgenomeGenes(WheatLineage.valueOf(subgenome));
        return subgenomeGenes.indexOf(geneName);
    }

    public boolean contain(String geneName){
        String subgenome=this.getSubgenome(geneName).name();
        List<String> subgenomeGenes=this.getSubgenomeGenes(WheatLineage.valueOf(subgenome));
        return subgenomeGenes.contains(geneName);
    }

    public int getRowNum(){
        return this.getTriadID().size();
    }

    public static void triadPosRecombinationRate(String pgfFile, String triadFile, String recombinationFile, String outFile){
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        Triad triad=new Triad(triadFile);
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
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            bw.write("TriadID\tChr\tPos\tRecombinationRate");
            bw.newLine();
            String recombination;
            int apos=-1;
            for (int i = 0; i < triad.getRowNum(); i++) {
                abd=triad.getTriad(i);
                for (int j = 0; j < 3; j++) {
                    geneIndex=pgf.getGeneIndex(abd[j]);
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
                    sb.append(triad.getTraidID(i)).append("\t");
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

    public static void transformTriadPos(String triadPosFile, String outFile){
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw = IOTool.getTextWriter(outFile)) {
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

    public static void recombinationRatePopMeanLoad(String triadPosFile, String popMeanLoadFile, String outFile){
        try (BufferedReader brTriadPos = IOTool.getReader(triadPosFile);
             BufferedReader brPopMeanLoad = IOTool.getReader(popMeanLoadFile);
             BufferedWriter bw = IOTool.getTextWriter(outFile)) {
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

    public static void triadModel(String triadPosFile, String triadPosRegionFile){
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw =IOTool.getTextWriter(triadPosRegionFile)) {
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
     *
     * @param triadPosFile
     * @param windowSize
     * @param stepSize
     * @param outFile
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
}
