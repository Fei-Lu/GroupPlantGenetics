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

public class Loter_runner extends LocalAncestry {

    /**
     * soft path
     */
    String condaActivatePath;
    String vcftools;

    /**
     * genotypePath file
     */
    GenotypeMetaData genotypeMetaData;

    /**
     * logFile
     */
    String logFilePath; // all log will be append to this logFile

    /**
     * outDir for multiple runs
     */
    String outDir; // outDir

    /**
     * Other parameter
     */
    int threadsNum;

    String[] subDir = {"loterGroup", "loterGroupVCF", "LAI"};
    File[] subDirFile;

    public Loter_runner(Builder builder){
        this.genotypeMetaData=builder.genotypeMetaData;
        this.logFilePath=builder.logFilePath;
        this.outDir=builder.outDir;
        this.condaActivatePath=builder.condaActivatePath;
        this.vcftools=builder.vcftools;
        this.threadsNum=builder.threadsNum;
        this.subDir=builder.subDir;
        this.makeSubDir();
    }

    public Loter_runner(String parameterFile){
        this.initialize(parameterFile);
        this.startRun();
    }

    public void startRun(){
        this.preparePop();
        this.splitPopGenotype();
        this.run_loter();
    }

    private void initialize(String parameterFile){
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("CondaActivate")){
                    this.condaActivatePath =temp.get(1);
                    continue;
                }
                if (line.startsWith("vcftools")){
                    this.vcftools=temp.get(1);
                    continue;
                }
                if (line.startsWith("GenotypeMetaDataPath")){
                    this.genotypeMetaData = new GenotypeMetaData(temp.get(1));
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
                if (line.startsWith("threadsNum")){
                    this.threadsNum=Integer.parseInt(temp.get(1));
                }
            }

        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        this.makeSubDir();
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
        String condaActivatePath = "/Users/xudaxing/anaconda3/bin/activate";
        String vcftools = "/Users/xudaxing/Software/vcftools_0.1.13/bin/vcftools";
        int threadsNum = 2;
        String[] subDir = {"loterGroup", "loterGroupVCF", "LAI"};

        public Builder(GenotypeMetaData genotypeMetaData, String logFilePath, String outDir){
            this.genotypeMetaData= genotypeMetaData;
            this.logFilePath=logFilePath;
            this.outDir=outDir;
        }

        public Builder condaActivatePath(String condaActivatePath){
            this.condaActivatePath=condaActivatePath;
            return this;
        }

        public Builder vcftools(String vcftools){
            this.vcftools=vcftools;
            return this;
        }

        public Builder threadsNum(int threadsNum){
            this.threadsNum=threadsNum;
            return this;
        }

        public Builder subDir(String[] subDir){
            this.subDir=subDir;
            return this;
        }

        public Loter_runner build(){
            return new Loter_runner(this);
        }

    }

    private void makeSubDir(){
        this.subDirFile = new File[subDir.length];
        for (int i = 0; i < this.subDir.length; i++) {
            this.subDirFile[i] = new File(outDir, this.subDir[i]);
            this.subDirFile[i].mkdir();
        }
    }

    private void preparePop(){
        List<String> referencePopList, taxaList;
        String admixedPop;
        BufferedWriter bw;
        try {
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                referencePopList = genotypeMetaData.referencePopList[i];
                for (String refPop:referencePopList){
                    taxaList = genotypeMetaData.getTaxaInfo(i).getTaxaListOf(refPop);
                    bw = IOTool.getWriter(new File(subDirFile[0], genotypeMetaData.genotypeID[i]+"."+refPop+".txt"));
                    bw.write(String.join("\n", taxaList));
                    bw.newLine();
                    bw.flush();
                    bw.close();
                }
                admixedPop = genotypeMetaData.admixedPop[i];
                bw = IOTool.getWriter(new File(subDirFile[0], genotypeMetaData.genotypeID[i]+"."+admixedPop+".txt"));
                taxaList = genotypeMetaData.getTaxaInfo(i).getTaxaListOf(admixedPop);
                bw.write(String.join("\n", taxaList));
                bw.newLine();
                bw.flush();
                bw.close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    private void splitPopGenotype(){
        new File(this.logFilePath).delete();
        StringBuilder sb = new StringBuilder();
        List<String> commandList = new ArrayList<>();
        List<String> refPop_admixedPop;
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            refPop_admixedPop = new ArrayList<>(genotypeMetaData.referencePopList[i]);
            refPop_admixedPop.add(genotypeMetaData.admixedPop[i]);
            for (String pop:refPop_admixedPop){
                sb.setLength(0);
                sb.append(vcftools).append(" --gzvcf ").append(genotypeMetaData.genotypePath[i]).append(" ");
                sb.append("--keep ").append(new File(this.subDirFile[0], genotypeMetaData.genotypeID[i]+"."+pop+".txt").getAbsolutePath());
                sb.append(" ");
                sb.append("--recode --out ").append(new File(this.subDirFile[1], genotypeMetaData.genotypeID[i]+"."+pop).getAbsolutePath());
                commandList.add(sb.toString());
            }

        }
        int exitCode;
        for (String command: commandList){
            exitCode = CommandUtils.runOneCommand(command, this.outDir, new File(this.logFilePath));
            assert exitCode == 0 : command+" run failed";
        }
    }

    private void run_loter(){
        long start = System.nanoTime();
        List<Callable<Integer>> callableList = new ArrayList<>();
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            StringBuilder sb = new StringBuilder();
            StringBuilder sbRef = new StringBuilder();
            sb.setLength(0);
            sb.append("source ").append(condaActivatePath).append(" loter && ");
            sb.append("loter_cli -r ");
            sbRef.setLength(0);
            for (String refPop : genotypeMetaData.referencePopList[i]){
                sbRef.append(new File(this.subDirFile[1], genotypeMetaData.genotypeID[i]+"."+refPop+".recode.vcf").getAbsolutePath());
                sbRef.append(" ");
            }
            sb.append(sbRef).append(" -a ").append(new File(this.subDirFile[1], genotypeMetaData.genotypeID[i]+"."+genotypeMetaData.admixedPop[i]+
                    ".recode.vcf").getAbsolutePath());
            sb.append(" -f vcf -o ").append(new File(this.subDirFile[2], genotypeMetaData.genotypeID[i]+".lai.txt"));
            callableList.add(()->CommandUtils.runMultipleCommands(sb.toString(), this.outDir, new File(this.logFilePath)));
        }

        List<Integer> results = CommandUtils.run_commands(callableList, threadsNum);
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i)!=0){
                System.out.println(results.get(i));
                System.out.println(genotypeMetaData.genotypeID[i]+" run failed");
            }
        }
        System.out.println("loter_runner completed in "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
    }

    @Override
    public double[][][][] extractLocalAncestry(){
        double[][][][] localAncestry = new double[genotypeMetaData.genotypeID.length][][][];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new double[genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new double[genotypeMetaData.nWayAdmixture[i]][];
            }
        }
        BufferedReader br;
        try {
            DoubleList[][][] localAnc = new DoubleList[genotypeMetaData.genotypeID.length][][];
            for (int i = 0; i < localAnc.length; i++) {
                localAnc[i] = new DoubleList[genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][];
                for (int j = 0; j < localAncestry[i].length; j++) {
                    localAnc[i][j] = new DoubleList[genotypeMetaData.nWayAdmixture[i]];
                    for (int k = 0; k < localAnc[i][j].length; k++) {
                        localAnc[i][j][k] = new DoubleArrayList();
                    }
                }
            }
            String line;
            List<String> temp;
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                br = IOTool.getReader(new File(subDirFile[2], genotypeMetaData.genotypeID[i]+".lai.txt"));
                int haplotypeIndex=0;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    for (int j = 0; j < genotypeMetaData.referencePopList[i].size(); j++) {
                        for (String s : temp) {
                            if (Integer.parseInt(s) == j) {
                                localAnc[i][haplotypeIndex][j].add(1);
                            } else {
                                localAnc[i][haplotypeIndex][j].add(0);
                            }
                        }
                    }
                    haplotypeIndex++;
                    br.readLine();
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
            throw new RuntimeException(e);
        }
        return localAncestry;
    }

    @Override
    public BitSet[][][] extractLocalAncestry_bitset(){
        BitSet[][][] localAncestry = new BitSet[genotypeMetaData.genotypeID.length][][];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new BitSet[genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][];
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
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                br = IOTool.getReader(new File(subDirFile[2], genotypeMetaData.genotypeID[i]+".lai.txt"));
                int haplotypeIndex=0;
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line, " ");
                    for (int j = 0; j < genotypeMetaData.referencePopList[i].size(); j++) {
                        int variantIndex = 0;
                        for (String s : temp) {
                            if (Integer.parseInt(s) == j) {
                                localAncestry[i][haplotypeIndex][j].set(variantIndex, true);
                            } else {
                                localAncestry[i][haplotypeIndex][j].set(variantIndex, false);
                            }
                            variantIndex++;
                        }
                    }
                    haplotypeIndex++;
                    br.readLine();
                }
                br.close();
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        return localAncestry;
    }

}
