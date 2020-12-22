package daxing.load.ancestralSite.complementary;

import daxing.common.*;
import daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal.IndividualLoadComplementary;
import daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal.SlidingWindowForLoadComplement;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.math.BigInteger;
import java.util.*;
import java.util.stream.IntStream;

/**
 * 该类用于存储每个个体在亚基因水平上syn non del在cds区域为derived的概率
 * exon vcf
 * 基因cds区域
 */
public class IndividualTaxonDerivedProbability {

    /**
     * only TreeValidatedPloidy hexaploid
     * @param inputDir retainTriadHexaploid
     * @param taxa_InfoDBFile
     * @param outFile
     */
    public static void getIndividualTaxonDerivedProbability(String inputDir, String taxa_InfoDBFile, String outFile){
        List<File> fileList= IOUtils.getVisibleFileListInDir(inputDir);
        BufferedReader br;
        Map<String, String> taxonPloidyMap=RowTableTool.getMap(taxa_InfoDBFile, 0, 3);
        Map<String,String> taxonTreeValidatedPloidyMap= RowTableTool.getMap(taxa_InfoDBFile,0, 14);
        Map<String, String> taxonSubspeciesMap=RowTableTool.getMap(taxa_InfoDBFile, 0, 15);
        Map<String, String> taxonFdBySubContinentMap=RowTableTool.getMap(taxa_InfoDBFile,0, 24);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Taxon\tSub\tSynCount\tNonCount\tDelCount\tDerivedSynCount" +
                    "\tDerivedNonCount\tDerivedDelCount\tHeterSynCount\tHeterNonCount" +
                    "\tHeterDelCount\tTreeValidatedPloidy" +
                    "\tSubspecies\tfdBySubContinent");
            bw.newLine();
            String line, taxonName, treeValidatedPloid, subspecies,fdBySubContinent;
            String subgenome;
            List<String> temp;
            int numSyn, numNonsyn, numDel, numDerivedSyn, numbDerivedNon, numDerivedDel;
            int numHeterInSyn, numHeterInNon, numHeterInDel;
            int subgenomeIndex;
            int[] totalNumDerivedSyn, totalNumDerivedNon, totalNumDerivedDel;
            int[] totalNumSyn, totalNumNon, totalNumDel;
            int[] totalNumHeterSyn, totalNumHeterNon, totalNumHeterDel;
            Ploidy ploidy;
            StringBuilder sb=new StringBuilder();
            String[] subArray={"A","B","D"};
            String sub="D";
            for (int i = 0; i < fileList.size(); i++) {
                br=IOTool.getReader(fileList.get(i));
                br.readLine();
                taxonName=PStringUtils.fastSplit(fileList.get(i).getName(), ".").get(0);
                ploidy=Ploidy.newInstanceFromSubChar(taxonPloidyMap.get(taxonName));
                treeValidatedPloid=taxonTreeValidatedPloidyMap.get(taxonName);
                subspecies=taxonSubspeciesMap.get(taxonName);
                fdBySubContinent=taxonFdBySubContinentMap.get(taxonName);
                int subNum=ploidy.getSubgenomewNum();
                totalNumSyn=new int[subNum];
                totalNumNon=new int[subNum];
                totalNumDel=new int[subNum];
                totalNumDerivedSyn=new int[subNum];
                totalNumDerivedNon=new int[subNum];
                totalNumDerivedDel=new int[subNum];
                totalNumHeterSyn=new int[subNum];
                totalNumHeterNon=new int[subNum];
                totalNumHeterDel=new int[subNum];
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    subgenome=temp.get(0).substring(8,9);
                    subgenomeIndex= Arrays.binarySearch(subArray, subgenome);
                    subgenomeIndex= subNum==1 ? subgenomeIndex-2 : subgenomeIndex;
                    numSyn=Integer.parseInt(temp.get(1));
                    numDerivedSyn=Integer.parseInt(temp.get(2));
                    numHeterInSyn=Integer.parseInt(temp.get(3));
                    numNonsyn=Integer.parseInt(temp.get(4));
                    numbDerivedNon=Integer.parseInt(temp.get(5));
                    numHeterInNon=Integer.parseInt(temp.get(6));
                    numDel=Integer.parseInt(temp.get(7));
                    numDerivedDel=Integer.parseInt(temp.get(8));
                    numHeterInDel=Integer.parseInt(temp.get(9));
                    totalNumSyn[subgenomeIndex]+=numSyn;
                    totalNumDerivedSyn[subgenomeIndex]+=numDerivedSyn;
                    totalNumHeterSyn[subgenomeIndex]+=numHeterInSyn;
                    totalNumNon[subgenomeIndex]+=numNonsyn;
                    totalNumDerivedNon[subgenomeIndex]+=numbDerivedNon;
                    totalNumHeterNon[subgenomeIndex]+=numHeterInNon;
                    totalNumDel[subgenomeIndex]+=numDel;
                    totalNumDerivedDel[subgenomeIndex]+=numDerivedDel;
                    totalNumHeterDel[subgenomeIndex]+=numHeterInDel;
                }
                br.close();
                for (int j = 0; j < subNum; j++) {
                    sb.setLength(0);
                    sb.append(taxonName).append("\t");
                    sub= subNum==1 ? sub : subArray[j];
                    sb.append(sub);
                    sb.append("\t").append(totalNumSyn[j]).append("\t");
                    sb.append(totalNumNon[j]).append("\t").append(totalNumDel[j]).append("\t");
                    sb.append(totalNumDerivedSyn[j]).append("\t").append(totalNumDerivedNon[j]).append("\t");
                    sb.append(totalNumDerivedDel[j]).append("\t");
                    sb.append(totalNumHeterSyn[j]).append("\t").append(totalNumHeterNon[j]).append("\t");
                    sb.append(totalNumHeterDel[j]).append("\t");
                    sb.append(treeValidatedPloid).append("\t").append(subspecies).append("\t").append(fdBySubContinent);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public enum LoadType{
        Syn(0), Non(1), Del(2);

        int index;
        LoadType(int i) {
            this.index=i;
        }

        public int getIndex() {
            return index;
        }
    }

    public enum ABD{
        A(0),B(1),D(2);

        int index;
        ABD(int index) {
            this.index=index;
        }

        public int getIndex() {
            return index;
        }
    }

    /**
     *
     * @param triadPosFile
     * @param window_triadIDNum
     * @param step_triadIDNum
     * @param p_ABD_ancestral
     * @param outFile
     */
    public static void calculateIndividualTaxonLocalDerivedProbability(String triadPosFile, String triadGeneFile,
                                                                       ABD triadPosSub, LoadType loadType,
                                                                       int window_triadIDNum,
                                                                       int step_triadIDNum,
                                                                       double[] p_ABD_ancestral, String outFile){
        Triads triads =new Triads(triadGeneFile);
        RowTableTool<String> triadPosFileTable=new RowTableTool<>(triadPosFile);
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            int[] subgenomeAIDArray={1,2,3,4,5,6,7};
            br.readLine();
            bw.write("TriadPosASubID\tWindow_Start\tWindow_End\tTriadIDNum\tTriadChrA\tWindow_Mid.TriadPosA\tObservedA" +
                    "\tObservedB\tObservedD\tObservedHomoAncestral\tExpectedA\tExpectedB\tExpectedD\tExpectedHomoAncestral" +
                    "\tp_value");
            bw.newLine();
            String line, triadID, chr, triadPosChr;
            int subgenomeID, subIDIndex=-1;
            List<String> temp[]=new List[subgenomeAIDArray.length];
            for (int i = 0; i < temp.length; i++) {
                temp[i]=new ArrayList<>();
            }
            while ((line=br.readLine())!=null){
                triadID=PStringUtils.fastSplit(line).get(0);
                chr=PStringUtils.fastSplit(line).get(2).substring(7,9);
                triadPosChr= IndividualLoadComplementary.getTriadPosChr(triads.getTraidGenes(triadID), triadPosSub.name(),chr);
                subgenomeID=Integer.parseInt(triadPosChr.substring(0,1));
                subIDIndex=Arrays.binarySearch(subgenomeAIDArray, subgenomeID);
                temp[subIDIndex].add(line);
            }
            SlidingWindowForLoadComplement windowTotal, windowDerived;
            Comparator<String> comparator= Comparator.comparing(s -> PStringUtils.fastSplit(s).get(5+triadPosSub.getIndex()*2));
            comparator= comparator.thenComparing(s -> Integer.parseInt(PStringUtils.fastSplit(s).get(6+triadPosSub.getIndex()*2)));
            comparator=comparator.thenComparing(s -> PStringUtils.fastSplit(s).get(20));
            double[] valueABD;
            List<String> tempA, tempB, tempD;
            TDoubleArrayList valueAList, valueBList, valueDList;
            TDoubleArrayList valueDerivedAList, valueDerivedBList, valueDerivedDList;
            StringBuilder sb=new StringBuilder();
            int[] observations;
            int cumulativeSumOfRowIndex=0;
            for (int i = 0; i < temp.length; i++) {
                Collections.sort(temp[i], comparator);
                windowTotal=new SlidingWindowForLoadComplement(temp[i].size(),window_triadIDNum*3, step_triadIDNum*3);
                windowDerived=new SlidingWindowForLoadComplement(temp[i].size(), window_triadIDNum*3, step_triadIDNum*3);
                for (int j = 0; j < temp[i].size(); j=j+3) {
                    valueABD=new double[3];
                    tempA=PStringUtils.fastSplit(temp[i].get(j));
                    tempB=PStringUtils.fastSplit(temp[i].get(j+1));
                    tempD=PStringUtils.fastSplit(temp[i].get(j+2));
                    valueABD[0]= Double.parseDouble(tempA.get(11+loadType.getIndex()*3))-Double.parseDouble(tempA.get(13+loadType.getIndex()*3));
                    valueABD[1]=Double.parseDouble(tempB.get(11+loadType.getIndex()*3))-Double.parseDouble(tempB.get(13+loadType.getIndex()*3));
                    valueABD[2]=Double.parseDouble(tempD.get(11+loadType.getIndex()*3))-Double.parseDouble(tempD.get(13+loadType.getIndex()*3));
                    windowTotal.addValue(j+1, valueABD);
                    valueABD=new double[3];
                    valueABD[0]=Double.parseDouble(tempA.get(12+loadType.getIndex()*3));
                    valueABD[1]=Double.parseDouble(tempB.get(12+loadType.getIndex()*3));
                    valueABD[2]=Double.parseDouble(tempD.get(12+loadType.getIndex()*3));
                    windowDerived.addValue(j+1, valueABD);
                }
                if (i > 0){
                    cumulativeSumOfRowIndex=cumulativeSumOfRowIndex+temp[i-1].size();
                }
                for (int j = 0; j < windowTotal.getWindowNum(); j++) {
                    valueAList=windowTotal.getWindow1Value(j);
                    valueBList=windowTotal.getWindow2Value(j);
                    valueDList=windowTotal.getWindow3Value(j);
                    observations=new int[4];
                    valueDerivedAList=windowDerived.getWindow1Value(j);
                    valueDerivedBList=windowDerived.getWindow2Value(j);
                    valueDerivedDList=windowDerived.getWindow3Value(j);
                    double total=valueAList.sum()+valueBList.sum()+valueDList.sum();
                    double homoAncestralNum=total-(valueDerivedAList.sum()+valueDerivedBList.sum()+valueDerivedDList.sum());
                    sb.setLength(0);
                    sb.append(subgenomeAIDArray[i]).append("\t").append(windowTotal.getWindowStarts()[j]).append("\t");
                    sb.append(windowTotal.getWindowEnds()[j]).append("\t").append(windowTotal.getCountInWindow(j)).append("\t");
                    int rowIndexOfTriadPos= windowTotal.getWindowStarts()[j]+(windowTotal.getWindowEnds()[j]-windowTotal.getWindowStarts()[j])/2;
                    int rowIndex=rowIndexOfTriadPos-1+cumulativeSumOfRowIndex;
                    sb.append(triadPosFileTable.getCell(rowIndex, 5+triadPosSub.getIndex()*2)).append("\t");
                    sb.append(triadPosFileTable.getCell(rowIndex, 6+triadPosSub.getIndex()*2)).append("\t");
                    sb.append(valueDerivedAList.sum()).append("\t").append(valueDerivedBList.sum()).append("\t");
                    sb.append(valueDerivedDList.sum()).append("\t").append(homoAncestralNum).append("\t");
                    double expectedA= total*p_ABD_ancestral[0];
                    double expectedB=total*p_ABD_ancestral[1];
                    double expectedD=total*p_ABD_ancestral[2];
                    double expectedHomoAncestral=total*(1-p_ABD_ancestral[0]-p_ABD_ancestral[1]-p_ABD_ancestral[2]);
                    sb.append(expectedA).append("\t").append(expectedB).append("\t");
                    sb.append(expectedD).append("\t").append(expectedHomoAncestral).append("\t");
                    observations[0]=(int) valueDerivedAList.sum();
                    observations[1]=(int) valueDerivedBList.sum();
                    observations[2]=(int) valueDerivedDList.sum();
                    observations[3]=(int) homoAncestralNum;
                    double p =getMultinomialDistribution_p(observations, p_ABD_ancestral);
                    sb.append(p);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static double getMultinomialDistribution_p(int[] observed, double[] p){
        int n= Arrays.stream(observed).sum();
        BigInteger a=factorialHavingLargeResult(n);
        BigInteger b=BigInteger.ONE;
        double c=1d;
        for (int i = 0; i < observed.length; i++) {
            b=b.multiply(factorialHavingLargeResult(observed[i]));
        }
        for (int i = 0; i < p.length; i++) {
            c=c*Math.pow(p[i], observed[i]);
        }
        BigInteger res= a.divide(b);
        double resD=res.doubleValue();
        return resD*c;
    }

    public static BigInteger factorialHavingLargeResult(int n){
        BigInteger result = BigInteger.ONE;
        for (int i = n; i > 0; i--) {
            result=result.multiply(BigInteger.valueOf(i));
        }
        return result;
    }

    public static void addTriadChrPos(String inputDir, String pgfFile, String triadGeneFile, String outDir){
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir);
        PGF pgf=new PGF(pgfFile);
        String[] outName= files.stream().map(File::getName).map(s->s.replaceAll(".txt", ".chr.pos.txt")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->addTriadChrPos(files.get(e), pgf, triadGeneFile, new File(outDir,
                outName[e])));
    }

    /**
     * triadPos 为基因中点
     * @param inputFile
     * @param pgf
     * @param outFile
     */
    private static void addTriadChrPos(File inputFile, PGF pgf, String triadGeneFile, File outFile){
        pgf.sortGeneByName();
        Triads triads =new Triads(triadGeneFile);
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            String[] geneNames;
            List<String> temp;
            List<String>[] tem;
            temp= PStringUtils.fastSplit(br.readLine());
            temp.add(3,"Chr\tPos\tTriadChrA\tTriadPosA\tTriadChrB\tTriadPosB\tTriadChrD\tTriadPosD");
            bw.write(String.join("\t",temp));
            bw.newLine();
            int[] geneIndex, posOnChrID, posOnChr, chrIDStart, chrIDEnd, chrID;
            int[] pos;
            String[] chrs, triadPosChr, triadID;
            StringBuilder sbA=new StringBuilder();
            StringBuilder sbB=new StringBuilder();
            StringBuilder sbD=new StringBuilder();
            String[] abd={"A","B","D"};
            while ((line=br.readLine())!=null){
                tem=new List[3];
                tem[0]=PStringUtils.fastSplit(line);
                tem[1]=PStringUtils.fastSplit(br.readLine());
                tem[2]=PStringUtils.fastSplit(br.readLine());
                geneNames=new String[3];
                triadID=new String[3];
                for (int i = 0; i < 3; i++) {
                    geneNames[i]=tem[i].get(2);
                    triadID[i]=tem[i].get(0);
                }
                geneIndex=new int[3];
                for (int i = 0; i < 3; i++) {
                    geneIndex[i]=pgf.getGeneIndex(geneNames[i]);
                }
                chrID = new int[3];
                chrIDStart=new int[3];
                chrIDEnd=new int[3];
                posOnChrID=new int[3];
                posOnChr=new int[3];
                pos=new int[3];
                chrs=new String[3];
                triadPosChr=new String[3];
                for (int i = 0; i < 3; i++) {
                    chrID[i]=pgf.getGene(geneIndex[i]).getGeneRange().chr;
                    chrIDStart[i]=pgf.getGene(geneIndex[i]).getGeneRange().start;
                    chrIDEnd[i]=pgf.getGene(geneIndex[i]).getGeneRange().end;
                    posOnChrID[i]=(chrIDStart[i]+chrIDEnd[i])/2;
                    posOnChr[i]= RefV1Utils.getPosOnChromosome(chrID[i], posOnChrID[i]);
                    pos[i]=posOnChr[i];
                    chrs[i]=RefV1Utils.getChromosome(chrID[i], posOnChrID[i]);
                    triadPosChr[i]=IndividualLoadComplementary.getTriadPosChr(triads.getTraidGenes(triadID[i]), abd[i], chrs[i]);
                }
                sbA.setLength(0);
                sbB.setLength(0);
                sbD.setLength(0);
//                temp.add(3,"Pos\tTriadPosA\tTriadPosB\tTriadPosD");
                sbA.append(chrs[0]).append("\t").append(pos[0]).append("\t").append(triadPosChr[0]).append("\t").append(pos[0]).append("\t").append(triadPosChr[1]).append("\t").append(pos[1]).append("\t").append(triadPosChr[2]).append("\t").append(pos[2]);
                sbB.append(chrs[1]).append("\t").append(pos[1]).append("\t").append(triadPosChr[0]).append("\t").append(pos[0]).append("\t").append(triadPosChr[1]).append("\t").append(pos[1]).append("\t").append(triadPosChr[2]).append("\t").append(pos[2]);
                sbD.append(chrs[2]).append("\t").append(pos[2]).append("\t").append(triadPosChr[0]).append("\t").append(pos[0]).append("\t").append(triadPosChr[1]).append("\t").append(pos[1]).append("\t").append(triadPosChr[2]).append("\t").append(pos[2]);
                tem[0].add(3, sbA.toString());
                tem[1].add(3, sbB.toString());
                tem[2].add(3, sbD.toString());
                for (int i = 0; i < 3; i++) {
                    bw.write(String.join("\t",tem[i]));
                    bw.newLine();
                }
                System.out.println();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * 不包含杂合子
     * @param individualTaxonDerivedProbabilityFile
     * @param taxonName hexaploid
     * @param loadType
     * @return
     */
    public static double[] getIndividualDerivedProbability(String individualTaxonDerivedProbabilityFile,
                                                           String taxa_InfoFile,
                                                           String taxonName, LoadType loadType){
        List<String> taxonList=RowTableTool.getColumnList(individualTaxonDerivedProbabilityFile, 0);
        Map<String, String> taxonPloidyMap=RowTableTool.getMap(taxa_InfoFile, 0, 3);
        Ploidy ploidy=Ploidy.newInstanceFromSubChar(taxonPloidyMap.get(taxonName));
        if (!ploidy.equals(Ploidy.HEXAPLOID)){
            System.out.println("Only support hexaploid taxon, program will quit");
            System.exit(1);
        }
        int rowIndexA=taxonList.indexOf(taxonName);
        RowTableTool<String> table=new RowTableTool<>(individualTaxonDerivedProbabilityFile);
        List<String>[] lineListABD=new List[3];
        for (int i = 0; i < lineListABD.length; i++) {
            lineListABD[i]=table.getRow(rowIndexA+i);
        }
        int[] totalCountABD=new int[3];
        int[] derivedCountABD=new int[3];
        int[] heterCountABD=new int[3];
        for (int i = 0; i < 3; i++) {
            totalCountABD[i]=Integer.parseInt(lineListABD[i].get(2+loadType.getIndex()));
            derivedCountABD[i]=Integer.parseInt(lineListABD[i].get(5+loadType.getIndex()));
            heterCountABD[i]=Integer.parseInt(lineListABD[i].get(8+loadType.getIndex()));
        }
        double[] p_abd_ancestral=new double[4];
        double total= totalCountABD[0]+totalCountABD[1]+totalCountABD[2]-heterCountABD[0]-heterCountABD[1]-heterCountABD[2];
        p_abd_ancestral[0]=derivedCountABD[0]/total;
        p_abd_ancestral[1]=derivedCountABD[1]/total;
        p_abd_ancestral[2]=derivedCountABD[2]/total;
        double ancestralNum=total-derivedCountABD[0]-derivedCountABD[1]-derivedCountABD[2];
        p_abd_ancestral[3]=ancestralNum/total;
        return p_abd_ancestral;
    }



}
