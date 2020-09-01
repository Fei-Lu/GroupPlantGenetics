package daxing.load.ancestralSite.complementary;

import daxing.common.IOTool;
import daxing.common.Ploidy;
import daxing.common.RowTableTool;
import daxing.common.Triad;
import daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal.IndividualLoadComplementary;
import daxing.load.ancestralSite.complementary.loadComplementaryGlobalLocal.SlidingWindowForLoadComplement;
import gnu.trove.list.array.TDoubleArrayList;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import umontreal.iro.lecuyer.probdistmulti.MultinomialDist;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;

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
        Syn, Non, Del
    }

    /**
     *
     * @param triadPosFile
     * @param window_triadIDNum
     * @param step_triadIDNum
     * @param p_ABD
     * @param outFile
     */
    public static void calculateIndividualTaxonLocalDerivedProbability(String triadPosFile, String triadGeneFile,
                                                                       int window_triadIDNum,
                                                                       int step_triadIDNum,
                                                                       double[] p_ABD, String outFile){
        Triad triad=new Triad(triadGeneFile);
        try (BufferedReader br = IOTool.getReader(triadPosFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            int[] subgenomeAIDArray={1,2,3,4,5,6,7};
            br.readLine();
            bw.write("TriadPosASubID\tWindow_Start\tWindow_End\tTriadIDNum\tObservedA\tObservedB\t" +
                    "ObservedD\tObservedHomoAncestral\tExpectedA\tExpectedB\tExpectedD\tExpectedHomoAncestral" +
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
                chr=PStringUtils.fastSplit(line).get(2).substring(8,10);
                triadPosChr= IndividualLoadComplementary.getTriadPosChr(triad.getTraidGenes(triadID), "A",chr);
                subgenomeID=Integer.parseInt(triadPosChr.substring(0,1));
                subIDIndex=Arrays.binarySearch(subgenomeAIDArray, subgenomeID);
                temp[subIDIndex].add(line);
            }
            SlidingWindowForLoadComplement window, windowDerived;
            Comparator<String> comparator=Comparator.comparing(s -> Integer.parseInt(PStringUtils.fastSplit(s).get(4)));
            comparator=comparator.thenComparing(s -> PStringUtils.fastSplit(s).get(16));
            double[] valueABD;
            List<String> tempA, tempB, tempD;
            TDoubleArrayList valueAList, valueBList, valueDList;
            TDoubleArrayList valueDerivedAList, valueDerivedBList, valueDerivedDList;
            StringBuilder sb=new StringBuilder();
            MultinomialDist multinomialDist;
            double[] pArray=new double[4];
            int[] observations;
            for (int i = 0; i < p_ABD.length; i++) {
                pArray[i]=p_ABD[i];
            }
            pArray[3]=1-p_ABD[0]-p_ABD[1]-p_ABD[2];
            for (int i = 0; i < temp.length; i++) {
                Collections.sort(temp[i], comparator);
                window=new SlidingWindowForLoadComplement(temp[i].size(),window_triadIDNum*3, step_triadIDNum*3);
                windowDerived=new SlidingWindowForLoadComplement(temp[i].size(), window_triadIDNum*3, step_triadIDNum*3);
                for (int j = 0; j < temp[i].size(); j=j+3) {
                    valueABD=new double[3];
                    tempA=PStringUtils.fastSplit(temp[i].get(j));
                    tempB=PStringUtils.fastSplit(temp[i].get(j+1));
                    tempD=PStringUtils.fastSplit(temp[i].get(j+2));
                    valueABD[0]=Double.parseDouble(tempA.get(13));
                    valueABD[1]=Double.parseDouble(tempB.get(13));
                    valueABD[2]=Double.parseDouble(tempD.get(13));
                    window.addValue(j+1, valueABD);
                    valueABD=new double[3];
                    valueABD[0]=Double.parseDouble(tempA.get(14));
                    valueABD[1]=Double.parseDouble(tempB.get(14));
                    valueABD[2]=Double.parseDouble(tempD.get(14));
                    windowDerived.addValue(j+1, valueABD);
                }
                for (int j = 0; j < window.getWindowNum(); j++) {
                    valueAList=window.getWindow1Value(j);
                    valueBList=window.getWindow2Value(j);
                    valueDList=window.getWindow3Value(j);
                    observations=new int[4];
                    valueDerivedAList=windowDerived.getWindow1Value(j);
                    valueDerivedBList=windowDerived.getWindow2Value(j);
                    valueDerivedDList=windowDerived.getWindow3Value(j);
                    double total=valueAList.sum()+valueBList.sum()+valueDList.sum();
                    double homoAncestralNum=total-(valueDerivedAList.sum()+valueDerivedBList.sum()+valueDerivedDList.sum());
                    sb.setLength(0);
                    sb.append(subgenomeAIDArray[i]).append("\t").append(window.getWindowStarts()[j]).append("\t");
                    sb.append(window.getWindowEnds()[j]).append("\t").append(window.getCountInWindow(j)).append("\t");
                    sb.append(valueDerivedAList.sum()).append("\t").append(valueDerivedBList.sum()).append("\t");
                    sb.append(valueDerivedDList.sum()).append("\t").append(homoAncestralNum).append("\t");
                    double expectedA=total*p_ABD[0];
                    double expectedB=total*p_ABD[1];
                    double expectedD=total*p_ABD[2];
                    double expectedHomoAncestral=total*(1-p_ABD[0]-p_ABD[1]-p_ABD[2]);
                    sb.append(expectedA).append("\t").append(expectedB).append("\t");
                    sb.append(expectedD).append("\t").append(expectedHomoAncestral).append("\t");
                    multinomialDist=new MultinomialDist((int)total, pArray);
                    observations[0]=(int) valueDerivedAList.sum();
                    observations[1]=(int) valueDerivedBList.sum();
                    observations[2]=(int) valueDerivedDList.sum();
                    observations[3]=(int) homoAncestralNum;
                    double p =multinomialDist.cdf(observations);
                    double[] mean=multinomialDist.getMean();
                    sb.append(p).append("\t").append(mean[0]).append("\t").append(mean[1]).append("\t");
                    sb.append(mean[2]).append("\t").append(mean[3]);
                    bw.write(sb.toString());
                    bw.newLine();
                }

            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
