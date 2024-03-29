package daxing.v2.localAncestryInfer.runner;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
import daxing.v2.localAncestryInfer.laidp.GenotypeTable;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.Callable;

public class ELAI_runner extends LocalAncestry {

    /**
     * Soft path
     */
    String softPath;

    GenotypeMetaData genotypeMetaData;

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

    public ELAI_runner(Builder builder){
        this.softPath=builder.softPath;
        this.genotypeMetaData=builder.genotypeMetaData;
        this.logFilePath=builder.logFilePath;
        this.outDir=builder.outDir;
        this.expectationMaximizationSteps=builder.expectationMaximizationSteps;
        this.threadsNum=builder.threadsNum;
        this.makeSubDir();
    }

    public ELAI_runner(String parameterFile){
        this.initialize(parameterFile);
        this.startRun();
    }

    public void startRun(){
        this.prepareFile();
        this.runELAI();
    }

    private void initialize(String parameterFile){
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("ELAI")){
                    this.softPath =temp.get(1);
                    continue;
                }
                if (line.startsWith("GenotypeMetaDataPath")){
                    this.genotypeMetaData=new GenotypeMetaData(temp.get(1));
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
        this.workingDir = new File[genotypeMetaData.genotypeID.length];
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            this.workingDir[i] = new File(this.outDir, genotypeMetaData.genotypeID[i]);
            this.workingDir[i].mkdir();
        }
    }

    public static class Builder{

        /**
         * Required parameters
         */
        GenotypeMetaData genotypeMetaData;
        String logFilePath;
        String outDir;

        /**
         * Optional parameters
         */
        String softPath = "/Users/xudaxing/Software/ELAI/elai-mac";
        int expectationMaximizationSteps = 2;
        int threadsNum = 2;

        public Builder(GenotypeMetaData genotypeMetaData, String logFilePath, String outDir){
            this.genotypeMetaData=genotypeMetaData;
            this.logFilePath=logFilePath;
            this.outDir=outDir;
        }

        public Builder softPath(String softPath){
            this.softPath=softPath;
            return this;
        }

        public Builder expectationMaximizationSteps(int expectationMaximizationSteps){
            this.expectationMaximizationSteps=expectationMaximizationSteps;
            return this;
        }

        public Builder threadsNum(int threadsNum){
            this.threadsNum=threadsNum;
            return this;
        }

        public ELAI_runner build(){
            return new ELAI_runner(this);
        }
    }

    private void prepareFile(){
        for (int i = 0; i < genotypeMetaData.genotypePath.length; i++) {
            ELAI_runner.prepareInputFile(genotypeMetaData.genotypePath[i], genotypeMetaData.genotypeID[i],
                    genotypeMetaData.taxaInfoPath[i],
                    genotypeMetaData.admixedPop[i], genotypeMetaData.referencePopList[i], workingDir[i],
                    new File(workingDir[i], genotypeMetaData.genotypeID[i]+".pos.txt"));
        }
    }

    private static void prepareInputFile(String genotypeFile, String genotypeID, String taxaInfoFile, String admixedPop,
                                        List<String> referencePopList, File outDirBimBam,
                                        File posFile){
        TaxaInfo taxaInfo = new TaxaInfo(taxaInfoFile);
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
                    for (String s : pop_sampleList) {
                        taxonIndex = genotypeTables[i].getTaxonIndex(s);
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
        long start = System.nanoTime();
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append(this.softPath).append(" ");
            int refPopBaseNum = 10;
            for (int j = 0; j < genotypeMetaData.referencePopList[i].size(); j++) {
                sb.append("-g ").append(new File(workingDir[i], genotypeMetaData.genotypeID[i])).append(".").append(genotypeMetaData.referencePopList[i].get(j)).append(".inp ");
                sb.append("-p ").append(refPopBaseNum).append(" ");
                refPopBaseNum++;
            }
            sb.append("-g ").append(new File(workingDir[i], genotypeMetaData.genotypeID[i]+"."+genotypeMetaData.admixedPop[i]+".inp")).append(" ");
            sb.append("-p 1 -pos ");
            sb.append(new File(workingDir[i], genotypeMetaData.genotypeID[i]+".pos.txt "));
            sb.append("-s ").append(expectationMaximizationSteps).append(" ");
            sb.append("-o ").append(genotypeMetaData.genotypeID[i]).append(" ");
            sb.append("-C ").append(genotypeMetaData.nWayAdmixture[i]).append(" ");
            sb.append("-c ").append(genotypeMetaData.nWayAdmixture[i]*5).append(" ");
            if (genotypeMetaData.timeSinceAdmixture[i] > 0){
                sb.append("-mg ").append(genotypeMetaData.timeSinceAdmixture[i]);
            }
            int finalI = i;
            callableList.add(()-> CommandUtils.runOneCommand(sb.toString(), workingDir[finalI].getAbsolutePath(), new File(logFilePath)));
        }

        List<Integer> results = CommandUtils.run_commands(callableList, threadsNum);
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i)!=0){
                System.out.println(results.get(i));
                System.out.println(genotypeMetaData.genotypeID[i]+" run failed");
            }
        }
        System.out.println("elaiRunner completed in "+ Benchmark.getTimeSpanSeconds(start)+" seconds");
    }

    @Override
    public double[][][][] extractLocalAncestry(){
        double[][][][] localAncestry = new double[genotypeMetaData.genotypeID.length][][][];
        TaxaInfo taxaInfo;
        for (int i = 0; i < localAncestry.length; i++) {
            taxaInfo = new TaxaInfo(genotypeMetaData.taxaInfoPath[i]);
            localAncestry[i] = new double[taxaInfo.getPopSampleSize(genotypeMetaData.admixedPop[i])][][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new double[genotypeMetaData.nWayAdmixture[i]][];
            }
        }
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            int ancestryPopIndex, snpIndex;
            File outputFile;
            DoubleList[][][] localAnc = new DoubleList[genotypeMetaData.genotypeID.length][][];
            for (int i = 0; i < localAnc.length; i++) {
                taxaInfo = new TaxaInfo(genotypeMetaData.taxaInfoPath[i]);
                localAnc[i] = new DoubleList[taxaInfo.getPopSampleSize(genotypeMetaData.admixedPop[i])][];
                for (int j = 0; j < localAncestry[i].length; j++) {
                    localAnc[i][j] = new DoubleList[genotypeMetaData.nWayAdmixture[i]];
                    for (int k = 0; k < localAnc[i][j].length; k++) {
                        localAnc[i][j][k] = new DoubleArrayList();
                    }
                }
            }
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                outputFile = new File(workingDir[i], "output");
                br = IOTool.getReader(new File(outputFile, genotypeMetaData.genotypeID[i]+".ps21.txt"));
                int haplotypeIndex=0;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    for (int j = 0; j < temp.size()-1; j++) {
                        ancestryPopIndex = j % genotypeMetaData.nWayAdmixture[i];
                        snpIndex  = j / genotypeMetaData.nWayAdmixture[i];
                        localAnc[i][haplotypeIndex][ancestryPopIndex].add(Double.parseDouble(temp.get(j)));
                    }
                    br.readLine();
                    haplotypeIndex ++;
                }
                br.close();
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

    @Override
    public BitSet[][][] extractLocalAncestry_bitset(){
        BitSet[][][] localAncestry = new BitSet[genotypeMetaData.genotypeID.length][][];
        TaxaInfo taxaInfo;
        for (int i = 0; i < localAncestry.length; i++) {
            taxaInfo = new TaxaInfo(genotypeMetaData.taxaInfoPath[i]);
            localAncestry[i] = new BitSet[taxaInfo.getPopSampleSize(genotypeMetaData.admixedPop[i])][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new BitSet[genotypeMetaData.nWayAdmixture[i]];
                for (int k = 0; k < localAncestry[i][j].length; k++) {
                    localAncestry[i][j][k] = new BitSet();
                }
            }
        }
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            int ancestryPopIndex, snpIndex;
            File outputFile;
            Double inferredValue;
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                outputFile = new File(workingDir[i], "output");
                br = IOTool.getReader(new File(outputFile, genotypeMetaData.genotypeID[i]+".ps21.txt"));
                int haplotypeIndex=0;
                boolean ancestryValue;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    for (int j = 0; j < temp.size()-1; j++) {
                        ancestryPopIndex = j % genotypeMetaData.nWayAdmixture[i];
                        snpIndex  = j / genotypeMetaData.nWayAdmixture[i];
                        inferredValue = Double.parseDouble(temp.get(j));
                        ancestryValue = inferredValue > 0.5 ? true : false;
                        localAncestry[i][haplotypeIndex][ancestryPopIndex].set(snpIndex, ancestryValue);
                    }
                    br.readLine();
                    haplotypeIndex ++;
                }
                br.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return localAncestry;
    }

}
