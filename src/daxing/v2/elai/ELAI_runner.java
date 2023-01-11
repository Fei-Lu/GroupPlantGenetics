package daxing.v2.elai;

import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.GenotypeTable;
import daxing.v2.localAncestryInfer.TaxaInfo;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class ELAI_runner {

    String eLaiPath;

    String[] genotypeID;

    String[] genotypePath;

    TaxaInfo taxaInfo;

    String logFilePath;

    String outDir_multipleRun;

    int[] nWayAdmixture;
    List<String>[] referencePopList;
    String[] admixedPop;

    int expectationMaximizationSteps;

    int threadsNum;

    File[] subDir;

    public ELAI_runner(String parameterFile){
        this.initialize(parameterFile);
        this.makeSubDir();
        this.prepareFile();
    }

    private void initialize(String parameterFile){
        List<String> genotypeIDList = new ArrayList<>();
        List<String> genotypePathList=new ArrayList<>();
        IntList nWayAdmixtureList = new IntArrayList();
        List<String> admixedPopList = new ArrayList<>();
        List<List<String>> refPopList = new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line, line_genotype;
            List<String> temp, temp_genotype, tem;
            BufferedReader brGenotype;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("ELAI")){
                    this.eLaiPath=temp.get(1);
                    continue;
                }
                if (line.startsWith("GenotypePath")){
                    brGenotype = IOTool.getReader(temp.get(1));
                    brGenotype.readLine();
                    while ((line_genotype=brGenotype.readLine())!=null){
                        temp_genotype = PStringUtils.fastSplit(line_genotype);
                        genotypeIDList.add(temp_genotype.get(0));
                        genotypePathList.add(temp_genotype.get(1));
                        nWayAdmixtureList.add(Integer.parseInt(temp_genotype.get(2)));
                        admixedPopList.add(temp_genotype.get(3));
                        tem = PStringUtils.fastSplit(temp_genotype.get(4), ",");
                        refPopList.add(tem);
                    }
                    brGenotype.close();
                    this.genotypeID = genotypeIDList.toArray(new String[0]);
                    this.genotypePath = genotypePathList.toArray(new String[0]);
                    this.nWayAdmixture=nWayAdmixtureList.toIntArray();
                    this.admixedPop = admixedPopList.toArray(new String[0]);
                    this.referencePopList = refPopList.toArray(new List[0]);
                    continue;
                }
                if (line.startsWith("TaxaInfoPath")){
                    this.taxaInfo = new TaxaInfo(temp.get(1));
                    continue;
                }
                if (line.startsWith("LogFilePath")){
                    this.logFilePath = temp.get(1);
                    continue;
                }
                if (line.startsWith("OutDir")){
                    this.outDir_multipleRun = temp.get(1);
                    continue;
                }
                if (line.startsWith("EMSteps")){
                    this.expectationMaximizationSteps=Integer.parseInt(temp.get(1));
                    continue;
                }
                if (line.startsWith("threadsNum")){
                    this.threadsNum=Integer.parseInt(temp.get(1));
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void makeSubDir(){
        this.subDir = new File[this.genotypeID.length];
        for (int i = 0; i < this.genotypeID.length; i++) {
            this.subDir[i] = new File(this.outDir_multipleRun, this.genotypeID[i]);
            this.subDir[i].mkdir();
        }
    }

    private void prepareFile(){
        for (int i = 0; i < this.genotypePath.length; i++) {
            ELAI_runner.prepareInputFile(genotypePath[i], genotypeID[i], this.taxaInfo, this.admixedPop[i],
                    referencePopList[i], subDir[i], new File(subDir[i], this.genotypeID[i]+".pos.txt"));
        }
    }

    public static void prepareInputFile(String genotypeFile, String genotypeID, TaxaInfo taxaInfo, String admixedPop,
                                        List<String> referencePopList, File outDirBimBam,
                                        File posFile){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
        String[] allPopOutName = new String[referencePopList.size()+1];
        allPopOutName[0] = admixedPop;
        for (int i = 0; i < referencePopList.size(); i++) {
            allPopOutName[i+1] = referencePopList.get(i);
        }
        int[][] taxaIndex_allPop = new int[allPopOutName.length][];
        List<String> pop_sampleList;
        for (int i = 0; i < taxaIndex_allPop.length; i++) {
            pop_sampleList = taxaInfo.getTaxaListOf(allPopOutName[i]);
            taxaIndex_allPop[i] = new int[pop_sampleList.size()];
            for (int j = 0; j < taxaIndex_allPop[i].length; j++) {
                taxaIndex_allPop[i][j] = genotypeTable.getTaxonIndex(pop_sampleList.get(j));
            }
        }
        GenotypeTable[] genotypeTables = new GenotypeTable[allPopOutName.length];
        for (int i = 0; i < genotypeTables.length; i++) {
            genotypeTables[i] = genotypeTable.getSubGenotypeTableByTaxa(taxaIndex_allPop[i]);
        }
        BufferedWriter[] bwBimBams = new BufferedWriter[allPopOutName.length];
        for (int i = 0; i < bwBimBams.length; i++) {
            bwBimBams[i] = IOTool.getWriter(new File(outDirBimBam, genotypeID+"."+allPopOutName[i]+".inp"));
        }
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOTool.getWriter(posFile);
            for (int i = 0; i < genotypeTable.getSiteNumber(); i++) {
                sb.setLength(0);
                sb.append(genotypeTable.getSnps()[i].getSnpID()).append(" ");
                sb.append(genotypeTable.getSnps()[i].getPos()).append(" ");
                sb.append(genotypeTable.getSnps()[i].getChr());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            int taxonIndex;
            char alleleBase;
            String[] subGenotypeTaxa;
            for (int i = 0; i < genotypeTables.length; i++) {
                pop_sampleList = taxaInfo.getTaxaListOf(allPopOutName[i]);
                sb.setLength(0);
                sb.append(pop_sampleList.size()).append(" =  \n");
                sb.append(genotypeTables[i].getSiteNumber()).append("\n");
                subGenotypeTaxa = genotypeTables[i].getTaxa();
                sb.append("IND,");
                sb.append(String.join(",", subGenotypeTaxa));
                bwBimBams[i].write(sb.toString());
                bwBimBams[i].newLine();
                for (int j = 0; j < genotypeTables[i].getSiteNumber(); j++) {
                    sb.setLength(0);
                    sb.append(genotypeTables[i].getSnps()[j].getSnpID()).append(",");
                    for (int k = 0; k < pop_sampleList.size(); k++) {
                        taxonIndex = genotypeTables[i].getTaxonIndex(pop_sampleList.get(k));
                        alleleBase = genotypeTables[i].getAlleleBase(j, taxonIndex);
                        sb.append(alleleBase).append(alleleBase).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bwBimBams[i].write(sb.toString());
                    bwBimBams[i].newLine();
                }
                bwBimBams[i].flush();
                bwBimBams[i].close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
