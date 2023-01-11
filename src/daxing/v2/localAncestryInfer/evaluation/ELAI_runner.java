package daxing.v2.localAncestryInfer.evaluation;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.GenotypeTable;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class ELAI_runner implements LocalAncestry {

    /**
     * Soft path
     */
    String softPath;


    /**
     * genotypePath file
     */
    String[] genotypeID;

    String[] genotypePath;

    int[] nWayAdmixture;

    String[] admixedPop;

    List<String>[] referencePopList;

    int[] timeSinceAdmixture;   // -1 means unknown


    /**
     * taxaInfo file
     */
    TaxaInfo taxaInfo;


    /**
     * logFile
     */
    String logFilePath;


    /**
     * outDir for multiple runs
     */
    String outDir;


    /**
     * soft parameter
     */
    int expectationMaximizationSteps;


    /**
     * Other parameter
     */
    int threadsNum;

    /**
     * working dir for single run
     */
    File[] workingDir;

    public ELAI_runner(String parameterFile){
        this.initialize(parameterFile);
        this.prepareFile();
        this.runELAI();
    }

    private void initialize(String parameterFile){
        List<String> genotypeIDList = new ArrayList<>();
        List<String> genotypePathList=new ArrayList<>();
        IntList nWayAdmixtureList = new IntArrayList();
        List<String> admixedPopList = new ArrayList<>();
        List<List<String>> refPopList = new ArrayList<>();
        IntList timeSinceAdmixture=new IntArrayList();
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line, line_genotype;
            List<String> temp, temp_genotype, tem;
            BufferedReader brGenotype;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("ELAI")){
                    this.softPath =temp.get(1);
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
                        int time = Integer.parseInt(temp_genotype.get(5));
                        timeSinceAdmixture.add(time < 0 ? -1 : time);
                    }
                    brGenotype.close();
                    this.genotypeID = genotypeIDList.toArray(new String[0]);
                    this.genotypePath = genotypePathList.toArray(new String[0]);
                    this.nWayAdmixture=nWayAdmixtureList.toIntArray();
                    this.admixedPop = admixedPopList.toArray(new String[0]);
                    this.referencePopList = refPopList.toArray(new List[0]);
                    this.timeSinceAdmixture = timeSinceAdmixture.toIntArray();
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
                    this.outDir = temp.get(1);
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
        this.makeSubDir();
    }

    private void makeSubDir(){
        this.workingDir = new File[this.genotypeID.length];
        for (int i = 0; i < this.genotypeID.length; i++) {
            this.workingDir[i] = new File(this.outDir, this.genotypeID[i]);
            this.workingDir[i].mkdir();
        }
    }

    private void prepareFile(){
        for (int i = 0; i < this.genotypePath.length; i++) {
            ELAI_runner.prepareInputFile(genotypePath[i], genotypeID[i], this.taxaInfo, this.admixedPop[i],
                    referencePopList[i], workingDir[i], new File(workingDir[i], this.genotypeID[i]+".pos.txt"));
        }
    }

    private static void prepareInputFile(String genotypeFile, String genotypeID, TaxaInfo taxaInfo, String admixedPop,
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

    private void runELAI(){
        new File(logFilePath).delete();
        List<Callable<Integer>> callableList = new ArrayList<>();
        for (int i = 0; i < this.genotypeID.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append(this.softPath).append(" ");
            int refPopBaseNum = 10;
            for (int j = 0; j < this.referencePopList[i].size(); j++) {
                sb.append("-g ").append(new File(workingDir[i], genotypeID[i])+"."+referencePopList[i].get(j)+".inp ");
                sb.append("-p ").append(refPopBaseNum).append(" ");
                refPopBaseNum++;
            }
            sb.append("-g ").append(new File(workingDir[i], genotypeID[i]+"."+admixedPop[i]+".inp")).append(" ");
            sb.append("-p 1 -pos ");
            sb.append(new File(workingDir[i], genotypeID[i]+".pos.txt "));
            sb.append("-s ").append(expectationMaximizationSteps).append(" ");
            sb.append("-o ").append(genotypeID[i]).append(" ");
            sb.append("-C ").append(nWayAdmixture[i]).append(" ");
            sb.append("-c ").append(nWayAdmixture[i]*5).append(" ");
            if (timeSinceAdmixture[i] > 0){
                sb.append("-mg ").append(timeSinceAdmixture[i]);
            }
            int finalI = i;
            callableList.add(()-> CommandUtils.runOneCommand(sb.toString(), workingDir[finalI].getAbsolutePath(), new File(logFilePath)));
        }

        List<Integer> results = CommandUtils.run_commands(callableList, threadsNum);
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i)!=0){
                System.out.println(results.get(i));
                System.out.println(genotypeID[i]+" run failed");
            }
        }
    }

    public double[][][][] extractLocalAncestry(){
        double[][][][] localAncestry = new double[genotypeID.length][][][];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new double[this.taxaInfo.getPopSampleSize(admixedPop[i])][][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new double[this.nWayAdmixture[i]][];
            }
        }
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            int ancestryPopIndex, snpIndex;
            File outputFile;
            DoubleList[][][] localAnc = new DoubleList[genotypeID.length][][];
            for (int i = 0; i < localAnc.length; i++) {
                localAnc[i] = new DoubleList[taxaInfo.getPopSampleSize(admixedPop[i])][];
                for (int j = 0; j < localAncestry[i].length; j++) {
                    localAnc[i][j] = new DoubleList[nWayAdmixture[i]];
                    for (int k = 0; k < localAnc[i][j].length; k++) {
                        localAnc[i][j][k] = new DoubleArrayList();
                    }
                }
            }
            for (int i = 0; i < this.genotypeID.length; i++) {
                outputFile = new File(workingDir[i], "output");
                br = IOTool.getReader(new File(outputFile, genotypeID[i]+".ps21.txt"));
                int haplotypeIndex=0;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    for (int j = 0; j < temp.size()-1; j++) {
                        ancestryPopIndex = j % this.nWayAdmixture[i];
                        snpIndex  = j / this.nWayAdmixture[i];
                        localAnc[i][haplotypeIndex][ancestryPopIndex].add(Double.parseDouble(temp.get(j)));
                    }
                    br.readLine();
                    haplotypeIndex ++;
                }
            }
            for (int i = 0; i < localAnc.length; i++) {
                for (int j = 0; j < localAnc[i].length; j++) {
                    for (int k = 0; k < localAnc[i][j].length; k++) {
                        localAncestry[i][j][k] = localAnc[i][j][k].toDoubleArray();
                    }
                }
            }
        } catch (Exception e) {
           e.printStackTrace();
        }
        return localAncestry;
    }
}
