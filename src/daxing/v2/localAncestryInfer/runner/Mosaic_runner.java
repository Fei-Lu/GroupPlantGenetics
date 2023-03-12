package daxing.v2.localAncestryInfer.runner;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;

/**
 * usage example
 * Mosaic_runner mosaic_runner = new Mosaic_runner(ParameterFile)
 * double[][][] localAncestry = mosaic_runner.getLocalAncestry()
 */
public class Mosaic_runner extends LocalAncestry {

    /**
     * soft path
     */
    String rscriptPath;

    /**
     * genotypePath file
     */
    GenotypeMetaData genotypeMetaData;

    String fastDirPath; // temp file for mosaic
    String logFile; // all log will be append to this logFile
    String outDir; // outDir
    int coresNumber; // core number, not threads number

    /**
     * working dir for single run
     */
    File[] workingDir;
    String subDirName = "mosaicData";

    public Mosaic_runner(Builder builder){
        this.rscriptPath = builder.rscriptPath;
        this.genotypeMetaData=builder.genotypeMetaData;
        this.fastDirPath=builder.fastDirPath;
        this.logFile=builder.logFile;
        this.outDir=builder.outDir;
        this.coresNumber=builder.coresNumber;
        this.subDirName= builder.subDirName;
        this.makeSubDir();
    }

    public Mosaic_runner(String parameterFile){
        this.initialize(parameterFile);
        this.startRun();
    }

    public void startRun(){
        this.prepareInputFile();
        this.run_mosaic();
        this.transformLocalAncestryRDataToTxt();
    }

    private void initialize(String parameterFile){
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("Rscript")){
                    this.rscriptPath =temp.get(1);
                    continue;
                }
                if (line.startsWith("GenotypeMetaDataPath")){
                    this.genotypeMetaData = new GenotypeMetaData(temp.get(1));
                    continue;
                }
                if (line.startsWith("FastDirPath")){
                    this.fastDirPath=temp.get(1);
                    continue;
                }
                if (line.startsWith("LogFilePath")){
                    this.logFile = temp.get(1);
                    continue;
                }
                if (line.startsWith("OutDir")){
                    this.outDir = temp.get(1);
                    continue;
                }
                if (line.startsWith("coresNumber")){
                    this.coresNumber=Integer.parseInt(temp.get(1));
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
            this.workingDir[i] = new File(outDir, genotypeMetaData.genotypeID[i]);
            this.workingDir[i].mkdir();
        }
    }

    public static class Builder{

        /**
         * Required parameters
         */
        GenotypeMetaData genotypeMetaData;
        String logFile; // all log will be append to this logFile
        String outDir; // outDir

        /**
         * Optional parameters
         *
         */
        String fastDirPath ="/Users/xudaxing/Documents/temp";
        String rscriptPath= "/usr/local/bin/Rscript";
        int coresNumber=2; // core number, not threads number
        String subDirName = "mosaicData";

        public Builder(GenotypeMetaData genotypeMetaData, String logFile, String outDir){
            this.genotypeMetaData=genotypeMetaData;
            this.logFile=logFile;
            this.outDir=outDir;
        }

        public Builder rscriptPath(String rscriptPath){
            this.rscriptPath=rscriptPath;
            return this;
        }

        public Builder fastDirPath(String fastDirPath){
            this.fastDirPath=fastDirPath;
            return this;
        }

        public Builder coresNumber(int coresNumber){
            this.coresNumber=coresNumber;
            return this;
        }

        public Builder subDirName(String subDirName){
            this.subDirName=subDirName;
            return this;
        }

        public Mosaic_runner build(){
            return new Mosaic_runner(this);
        }
    }

    private void prepareInputFile(){
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            this.transformVCF_to_Mosaic_inputFormat(i);
        }
    }

    private void transformVCF_to_Mosaic_inputFormat(int indexOfRun){
        String[] subDir = {this.subDirName};
        File[] subDirFile = new File[subDir.length];
        for (int i = 0; i < subDir.length; i++) {
            subDirFile[i] = new File(workingDir[indexOfRun], subDir[i]);
            subDirFile[i].mkdir();
        }
        IntList posList = new IntArrayList();
        TaxaInfo taxaInfo = this.genotypeMetaData.getTaxaInfo(indexOfRun);
        DoubleList geneticsMapList = new DoubleArrayList();
        try (BufferedReader brGenotype = IOTool.getReader(genotypeMetaData.genotypePath[indexOfRun]);
             BufferedReader brRecombinationMap = IOTool.getReader(genotypeMetaData.recombinationMap[indexOfRun])) {
            String line, popName;
            List<String> temp, headerList;
            StringBuilder sb = new StringBuilder();

            brRecombinationMap.readLine();
            while ((line=brRecombinationMap.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                posList.add(Integer.parseInt(temp.get(0)));
                geneticsMapList.add(Double.parseDouble(temp.get(2)));
            }

            // rates.1
            File rate1 = new File(subDirFile[0], "rates.1");
            BufferedWriter bwRates1 = IOTool.getWriter(rate1);
            sb.setLength(0);
            sb.append(":sites:").append(posList.size());
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            sb.setLength(0);
            for (int i = 0; i < posList.size(); i++) {
                sb.append(posList.getInt(i)).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            sb.setLength(0);
            for (int i = 0; i < geneticsMapList.size(); i++) {
                sb.append(geneticsMapList.getDouble(i)).append(" ");
            }
            sb.deleteCharAt(sb.length()-1);
            bwRates1.write(sb.toString());
            bwRates1.newLine();
            bwRates1.flush();
            bwRates1.close();

            // pop2TaxaMap
            while ((line=brGenotype.readLine()).startsWith("##")){}
            headerList = PStringUtils.fastSplit(line);
            Map<String,IntList> popIndexMap = new HashMap<>();
            List<String> refPop_AdmixedPopList = new ArrayList<>(genotypeMetaData.referencePopList[indexOfRun]);
            refPop_AdmixedPopList.add(genotypeMetaData.admixedPop[indexOfRun]);
            for (String value : refPop_AdmixedPopList) {
                popIndexMap.put(value, new IntArrayList());
            }
            String taxonName;
            for (int i = 9; i < headerList.size(); i++) {
                taxonName = headerList.get(i);
                popName = taxaInfo.getTaxonPop(taxonName);
                popIndexMap.get(popName).add(i);
            }

            // sample.names
            File sampleNamesFile = new File(subDirFile[0], "sample.names");
            BufferedWriter bwSampleName = IOTool.getWriter(sampleNamesFile);
            IntList indexList;
            for (int i = 0; i < refPop_AdmixedPopList.size(); i++) {
                indexList = popIndexMap.get(refPop_AdmixedPopList.get(i));
                for (int j = 0; j < indexList.size(); j++) {
                    sb.setLength(0);
                    sb.append(refPop_AdmixedPopList.get(i)).append(" ");
                    sb.append(headerList.get(indexList.getInt(j))).append(" ");
                    sb.append("0 0 0 2 -9");
                    bwSampleName.write(sb.toString());
                    bwSampleName.newLine();
                }
            }
            bwSampleName.flush();
            bwSampleName.close();

            // snpfile and genofile
            File snpFile1 = new File(subDirFile[0], "snpfile.1");
            BufferedWriter bwSNPFile = IOTool.getWriter(snpFile1);
            Map<String,BufferedWriter> pop2BWMap = new HashMap<>();
            File genoFile;
            BufferedWriter bw;
            for (String s : refPop_AdmixedPopList) {
                genoFile = new File(subDirFile[0], s + "genofile.1");
                bw = IOTool.getWriter(genoFile);
                pop2BWMap.put(s, bw);
            }
            IntList indexesList;
            while ((line=brGenotype.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append("         ").append(temp.get(2)).append(" ");
                sb.append(temp.get(0)).append(" ").append(0).append(" ");
                sb.append(temp.get(1)).append(" ");
                sb.append(temp.get(3)).append(" ");
                sb.append(temp.get(4)).append(" ");
                bwSNPFile.write(sb.toString());
                bwSNPFile.newLine();

                for (Map.Entry<String,IntList> entry: popIndexMap.entrySet()){
                    popName = entry.getKey();
                    indexesList = entry.getValue();
                    bw = pop2BWMap.get(popName);
                    sb.setLength(0);
                    for (int index : indexesList){
                        sb.append(temp.get(index));
                    }
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bwSNPFile.flush();
            bwSNPFile.close();

            for (Map.Entry<String,BufferedWriter> entry : pop2BWMap.entrySet()){
                entry.getValue().flush();
                entry.getValue().close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void run_mosaic(){
        // make sure logFile will be create from first call
        new File(this.logFile).delete();
        String mosaicRpath = Mosaic_runner.class.getResource("mosaic.R").getPath();
        List<Callable<Integer>> callableList = new ArrayList<>();
        long start = System.nanoTime();
        for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append(rscriptPath).append(" ").append(mosaicRpath).append(" ");
            sb.append(genotypeMetaData.admixedPop[i]).append(" ").append(new File(workingDir[i], subDirName)).append("/ ");
            sb.append("-f ").append(fastDirPath).append(" ");
            sb.append("-a ").append(genotypeMetaData.nWayAdmixture[i]).append(" ");
            sb.append("-m ").append(this.coresNumber).append(" ");
            sb.append("-c ").append(genotypeMetaData.chrID[i]).append(":").append(genotypeMetaData.chrID[i]);
            int finalI = i;
            callableList.add(()-> CommandUtils.runOneCommand(sb.toString(), workingDir[finalI].getAbsolutePath(), new File(logFile)));
        }

        List<Integer> results = CommandUtils.run_commands(callableList, coresNumber*2);
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i)!=0){
                System.out.println(results.get(i));
                System.out.println(genotypeMetaData.genotypeID[i]+" run failed");
            }
        }
        System.out.println("mosaic_runner completed in "+ Benchmark.getTimeSpanSeconds(start)+" seconds");
    }

    private void transformLocalAncestryRDataToTxt(){
        String extractMosaicLocalAncestryRpath = Mosaic_runner.class.getResource("extractMosaicLocalAncestry.R").getPath();
        File currentWorkingDir;
        File mosaicInputDataDir;
        try {
            List<Callable<Integer>> callableList = new ArrayList<>();
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                StringBuilder sb = new StringBuilder();
                currentWorkingDir = workingDir[i];
                mosaicInputDataDir = currentWorkingDir.listFiles(dir -> dir.getName().equals(this.subDirName))[0];
                File mosaicResDir = currentWorkingDir.listFiles(dir -> dir.getName().startsWith("MOSAIC_RESULTS"))[0];
                File[] rDataFiles = mosaicResDir.listFiles(dir -> dir.getName().endsWith("RData"));
                String modelParametersRData = Arrays.stream(rDataFiles).filter(file -> !file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
                String localAncestryRData = Arrays.stream(rDataFiles).filter(file -> file.getName().startsWith("localanc")).toArray(File[]::new)[0].getAbsolutePath();
                sb.setLength(0);
                sb.append("Rscript ").append(extractMosaicLocalAncestryRpath).append(" ");
                sb.append(modelParametersRData).append(" ").append(localAncestryRData).append(" ");
                sb.append(mosaicInputDataDir).append("/ ");
                sb.append(mosaicResDir).append("/");
                int finalI = i;
                callableList.add(()->CommandUtils.runOneCommand(sb.toString(), workingDir[finalI].getAbsolutePath(), new File(localAncestryRData)));
            }
            List<Integer> results = CommandUtils.run_commands(callableList, coresNumber*2);
            for (int i = 0; i < results.size(); i++) {
                if (results.get(i)!=0){
                    System.out.println(results.get(i));
                    System.out.println(genotypeMetaData.genotypeID[i]+" run failed");
                }
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public double[][][][] extractLocalAncestry() {
        double[][][][] localAncestry = new double[genotypeMetaData.genotypeID.length][][][];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new double[this.genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new double[genotypeMetaData.nWayAdmixture[i]][];
            }
        }
        DoubleList[][][] localAnc = new DoubleList[genotypeMetaData.genotypeID.length][][];
        File currentWorkingDir;
        File mosaicResDir;
        List<File> ancestryFiles;
        DoubleList[][] haplotype_source_variants_array;
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                currentWorkingDir = workingDir[i];
                mosaicResDir = currentWorkingDir.listFiles(dir -> dir.getName().startsWith("MOSAIC_RESULTS"))[0];
                ancestryFiles = IOTool.getFileListInDirStartsWith(mosaicResDir.getAbsolutePath(), "ancestry");
                haplotype_source_variants_array = new DoubleList[genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][];
                for (int j = 0; j < haplotype_source_variants_array.length; j++) {
                    haplotype_source_variants_array[j] = new DoubleList[genotypeMetaData.nWayAdmixture[i]];
                    for (int k = 0; k < haplotype_source_variants_array[j].length; k++) {
                        haplotype_source_variants_array[j][k] = new DoubleArrayList();
                    }
                }
                for (int j = 0; j < ancestryFiles.size(); j++) {
                    br = IOTool.getReader(ancestryFiles.get(j));
                    br.readLine();
                    while ((line=br.readLine())!=null){
                        temp = PStringUtils.fastSplit(line);
                        for (int k = 0; k < temp.size(); k++) {
                            haplotype_source_variants_array[k][j].add(Double.parseDouble(temp.get(k)));
                        }
                    }
                    br.close();
                }
                localAnc[i] = haplotype_source_variants_array;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        for (int i = 0; i < localAnc.length; i++) {
            for (int j = 0; j < localAnc[i].length; j++) {
                for (int k = 0; k < localAnc[i][j].length; k++) {
                    localAncestry[i][j][k] = localAnc[i][j][k].toDoubleArray();
                }
            }
        }
        return localAncestry;
    }

    @Override
    public BitSet[][][] extractLocalAncestry_bitset() {
        BitSet[][][] localAncestry = new BitSet[genotypeMetaData.genotypeID.length][][];
        for (int i = 0; i < localAncestry.length; i++) {
            localAncestry[i] = new BitSet[this.genotypeMetaData.getTaxaInfo(i).getPopSampleSize(genotypeMetaData.admixedPop[i])][];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new BitSet[genotypeMetaData.nWayAdmixture[i]];
                for (int k = 0; k < localAncestry[i][j].length; k++) {
                    localAncestry[i][j][k] = new BitSet();
                }
            }
        }
        File currentWorkingDir;
        File mosaicResDir;
        List<File> ancestryFiles;
        BufferedReader br;
        try {
            String line;
            List<String> temp;
            for (int i = 0; i < genotypeMetaData.genotypeID.length; i++) {
                currentWorkingDir = workingDir[i];
                mosaicResDir = currentWorkingDir.listFiles(dir -> dir.getName().startsWith("MOSAIC_RESULTS"))[0];
                ancestryFiles = IOTool.getFileListInDirStartsWith(mosaicResDir.getAbsolutePath(), "ancestry");
                for (int j = 0; j < ancestryFiles.size(); j++) {
                    br = IOTool.getReader(ancestryFiles.get(j));
                    br.readLine();
                    int variantIndex = 0 ;
                    boolean ancestryValue;
                    while ((line=br.readLine())!=null){
                        temp = PStringUtils.fastSplit(line);
                        for (int k = 0; k < temp.size(); k++) {
                            ancestryValue = Double.parseDouble(temp.get(k)) > 0.5 ? true : false;
                            localAncestry[i][k][j].set(variantIndex, ancestryValue);
                        }
                        variantIndex++;
                    }
                    br.close();
                }
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        return localAncestry;
    }


}
