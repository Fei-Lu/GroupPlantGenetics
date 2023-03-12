package daxing.v2.localAncestryInfer.runner;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
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

public class LAIDP_runner extends LocalAncestry {

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
    int windowSize;
    int stepSize;
    String ancestralAllele;
    double switchCostScore;
    int conjunctionNum;


    /**
     * Other parameter
     */
    int threadsNum;

    /**
     * working dir for single run
     */
    File[] workingDir;

    public LAIDP_runner(Builder builder){
        this.softPath=builder.softPath;
        this.genotypeMetaData=builder.genotypeMetaData;
        this.logFilePath=builder.logFilePath;
        this.outDir=builder.outDir;
        this.windowSize=builder.windowSize;
        this.stepSize=builder.stepSize;
        this.ancestralAllele=builder.ancestralAllele;
        this.switchCostScore=builder.switchCostScore;
        this.conjunctionNum=builder.conjunctionNum;
        this.threadsNum=builder.threadsNum;
        this.makeSubDir();
    }

    private void makeSubDir(){
        this.workingDir = new File[genotypeMetaData.genotypeID.length];
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            this.workingDir[i] = new File(this.outDir, genotypeMetaData.genotypeID[i]);
            this.workingDir[i].mkdir();
        }
    }

    public void startRun(){
        this.prepareFile();
        this.runLAIDP();
    }

    private void prepareFile(){
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            this.prepareFile(genotypeMetaData.taxaInfoPath[i], genotypeMetaData.admixedPop[i],
                    genotypeMetaData.nativePop[i], genotypeMetaData.introgressedPop[i],
                    new File(workingDir[i], genotypeMetaData.genotypeID[i]+".taxaGroup.txt"));
        }
    }

    private void prepareFile(String taxaInfoFile, String admixedPop,
                             String nativePop, List<String> introgressedPopList, File taxaGroupFile){
        TaxaInfo taxaInfo = new TaxaInfo(taxaInfoFile);
        try (BufferedWriter bw = IOTool.getWriter(taxaGroupFile)) {
            List<String> pop_sampleList;
            StringBuilder sb = new StringBuilder();
            sb.append("Taxon\tPopulationName\tSource");
            bw.write(sb.toString());
            bw.newLine();
            pop_sampleList = taxaInfo.getTaxaListOf(admixedPop);
            for (String taxon : pop_sampleList){
                sb.setLength(0);
                sb.append(taxon).append("\t").append(admixedPop).append("\t").append("ADMIXED");
                bw.write(sb.toString());
                bw.newLine();
            }
            pop_sampleList = taxaInfo.getTaxaListOf(nativePop);
            for (String taxon : pop_sampleList){
                sb.setLength(0);
                sb.append(taxon).append("\t").append(nativePop).append("\t").append("NATIVE");
                bw.write(sb.toString());
                bw.newLine();
            }
            int introgressedIndex = 1;
            for (String pop : introgressedPopList){
                pop_sampleList = taxaInfo.getTaxaListOf(pop);
                for (String taxon : pop_sampleList){
                    sb.setLength(0);
                    sb.append(taxon).append("\t").append(pop).append("\t");
                    sb.append("INTROGRESSED").append("_").append(introgressedIndex);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                introgressedIndex++;
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void runLAIDP(){
        new File(logFilePath).delete();
        List<Callable<Integer>> callableList = new ArrayList<>();
        long start = System.nanoTime();
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append("java -jar ").append(this.softPath).append(" ");
            sb.append(genotypeMetaData.genotypePath[i]).append(" ");
            sb.append(this.windowSize).append(" ");
            sb.append(this.stepSize).append(" ");
            sb.append(new File(workingDir[i], genotypeMetaData.genotypeID[i]+".taxaGroup.txt").getAbsolutePath()).append(" ");
            sb.append(this.ancestralAllele).append(" ");
            sb.append(this.conjunctionNum).append(" ");
            sb.append(this.switchCostScore).append(" ");
            sb.append(new File(workingDir[i], genotypeMetaData.genotypeID[i]+".localAnc.txt")).append(" ");
            sb.append(this.threadsNum);
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
        System.out.println("laidpRunner completed in "+ Benchmark.getTimeSpanSeconds(start)+" seconds");
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
        String softPath= "/Users/xudaxing/Software/LAIDP/LAIDP_single.jar";
        int windowSize = 200;
        int stepSize = 100;
        String ancestralAllele = "simulation";
        double switchCostScore = 1.5;
        int conjunctionNum= 2;
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

        public Builder windowSize(int windowSize){
            this.windowSize=windowSize;
            return this;
        }

        public Builder stepSize(int stepSize){
            this.stepSize=stepSize;
            return this;
        }

        public Builder ancestralAllele(String ancestralAllele){
            this.ancestralAllele=ancestralAllele;
            return this;
        }

        public Builder switchCostScore(double switchCostScore){
            this.switchCostScore=switchCostScore;
            return this;
        }

        public Builder conjunctionNum(int conjunctionNum){
            this.conjunctionNum=conjunctionNum;
            return this;
        }

        public Builder threadsNum(int threadsNum){
            this.threadsNum=threadsNum;
            return this;
        }

        public LAIDP_runner build(){
            return new LAIDP_runner(this);
        }
    }

    @Override
    public double[][][][] extractLocalAncestry() {
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
            List<String> temp, tem;
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
                br = IOTool.getReader(new File(workingDir[i], genotypeMetaData.genotypeID[i]+".localAnc.txt"));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    for (int j = 1; j < temp.size(); j++) {
                        tem = PStringUtils.fastSplit(temp.get(j), ",");
                        for (int k = 0; k < tem.size(); k++) {
                            if (k == 0){
                                // k=0, native ancestry
                                localAnc[i][j-1][tem.size()-1].add(Integer.parseInt(tem.get(k)));
                            }else {
                                // k > 0, introgressed ancestry
                                localAnc[i][j-1][k-1].add(Integer.parseInt(tem.get(k)));
                            }
                        }
                    }
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
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return localAncestry;
    }

    @Override
    public BitSet[][][] extractLocalAncestry_bitset() {
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
            List<String> temp, tem;
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                int variantIndex = 0;
                boolean ancestryValue;
                br = IOTool.getReader(new File(workingDir[i], genotypeMetaData.genotypeID[i]+".localAnc.txt"));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    for (int j = 1; j < temp.size(); j++) {
                        tem = PStringUtils.fastSplit(temp.get(j), ",");
                        for (int k = 0; k < tem.size(); k++) {
                            ancestryValue = Integer.parseInt(tem.get(k)) > 0.5 ? true : false;
                            if (k == 0){
                                // k=0, native ancestry
                                localAncestry[i][j-1][tem.size()-1].set(variantIndex, ancestryValue);
                            }else {
                                // k > 0, introgressed ancestry
                                localAncestry[i][j-1][k-1].set(variantIndex, ancestryValue);
                            }
                        }
                    }
                    variantIndex++;
                }
                br.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return localAncestry;
    }
}
