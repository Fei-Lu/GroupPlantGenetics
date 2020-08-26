package daxing.load.ancestralSite.complementary;

import daxing.common.IOTool;
import daxing.common.Ploidy;
import daxing.common.RowTableTool;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

/**
 * 该类用于存储每个个体在亚基因水平上syn non del在cds区域为derived的概率
 * exon vcf
 * 基因cds区域
 */
public class IndividualTaxonDerivedProbability {

    /**
     *
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
        try (BufferedWriter bw = IOTool.getTextWriter(outFile)) {
            bw.write("Taxon\tSub\tSynCount\tNonCount\tDelCount\tDerivedSynCount" +
                    "\tDerivedNonCount\tDerivedDelCount\tTreeValidatedPloidy\tSubspecies\tfdBySubContinent");
            bw.newLine();
            String line, taxonName, treeValidatedPloid, subspecies,fdBySubContinent;
            List<String>[] temp;
            int numSyn, numNonsyn, numDel, numDerivedSyn, numbDerivedNon, numDerivedDel;
            int[] totalNumDerivedSyn, totalNumDerivedNon, totalNumDerivedDel;
            int[] totalNumSyn, totalNumNon, totalNumDel;
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
                while ((line=br.readLine())!=null){
                    temp=new List[subNum];
                    temp[0]=PStringUtils.fastSplit(line);
                    for (int j = 1; j < temp.length; j++) {
                        temp[j]=PStringUtils.fastSplit(br.readLine());
                    }
                    for (int j = 0; j < temp.length; j++) {
                        numSyn=Integer.parseInt(temp[j].get(3));
                        numDerivedSyn=Integer.parseInt(temp[j].get(4));
                        numNonsyn=Integer.parseInt(temp[j].get(5));
                        numbDerivedNon=Integer.parseInt(temp[j].get(6));
                        numDel=Integer.parseInt(temp[j].get(7));
                        numDerivedDel=Integer.parseInt(temp[j].get(8));
                        totalNumSyn[j]+=numSyn;
                        totalNumDerivedSyn[j]+=numDerivedSyn;
                        totalNumNon[j]+=numNonsyn;
                        totalNumDerivedNon[j]+=numbDerivedNon;
                        totalNumDel[j]+=numDel;
                        totalNumDerivedDel[j]+=numDerivedDel;
                    }
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

}
