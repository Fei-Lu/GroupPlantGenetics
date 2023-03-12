package daxing.v2.localAncestryInfer.evaluation;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.laidp.GenotypeTable;
import daxing.v2.localAncestryInfer.simulation.SimulationMetadata;
import org.apache.commons.math3.stat.correlation.PearsonsCorrelation;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class LocalAncestry {

    static int[] variantNum;

    /**
     * @return local ancestry matrix,
     * dim1 is different run,
     * dim2 is different haplotype
     * dim3 is source population,  first all introgressed populations, then native population
     * dim4 is variants
     */
    abstract protected double[][][][] extractLocalAncestry();

    abstract protected BitSet[][][] extractLocalAncestry_bitset();

    public int[][][] contingencyTable(double[][][][] actual_values){
        double[][][][] inferredValue = this.extractLocalAncestry();
        return LocalAncestry.contingencyTable(inferredValue, actual_values);
    }

    public int[][][] contingencyTable_bitset(BitSet[][][] actual_values){
        BitSet[][][] inferredValue = this.extractLocalAncestry_bitset();
        return LocalAncestry.contingencyTable_bitset(inferredValue, actual_values);
    }

    public int[][][] contingencyTable_2way(double[][][][] actual_values){
        double[][][][] inferredValue = this.extractLocalAncestry();
        return LocalAncestry.contingencyTable_2way(inferredValue, actual_values);
    }

    public double[][] pearsonCorrelation(double[][][][] actual_values){
        double[][][][] inferredValue = this.extractLocalAncestry();
        return LocalAncestry.pearsonCorrelation(inferredValue, actual_values);
    }

    public double[][] meanDeviation(double[][][][] actual_values){
        double[][][][] inferredValue = this.extractLocalAncestry();
        return LocalAncestry.meanDeviation(inferredValue, actual_values);
    }

    /**
     *
     * @param inferredValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @param actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @return Mean deviation, dim1 is different run, dim2 is different haplotype
     */
    public static double[][] meanDeviation(double[][][][] inferredValue, double[][][][] actualValue){
        double[][] meanDeviation = new double[inferredValue.length][inferredValue[0].length];
        double deviationSum, count;
        for (int i = 0; i < inferredValue.length; i++) {
            for (int j = 0; j < inferredValue[i].length; j++) {
                deviationSum = 0;
                count = 0;
                for (int k = 0; k < inferredValue[i][j].length; k++) {
                    for (int l = 0; l < inferredValue[i][j][k].length; l++) {
                        deviationSum+=Math.abs(actualValue[i][j][k][l] - inferredValue[i][j][k][l]);
                        count++;
                    }
                }
                meanDeviation[i][j] = deviationSum/count;
            }
        }
        return meanDeviation;
    }

    /**
     *
     * @param inferredValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @param actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population, dim4 is variants
     * @return Pearson correlation, dim1 is different run, dim2 is different haplotype
     */
    public static double[][] pearsonCorrelation(double[][][][] inferredValue, double[][][][] actualValue){
        List<Callable<Double>> callableList= new ArrayList<>();
        for (int i = 0; i < inferredValue.length; i++) {
            for (int j = 0; j < inferredValue[i].length; j++) {
                int cumCount=0;
                double[] inferredValueArray = new double[inferredValue[i][j].length * inferredValue[i][j][0].length];
                double[] actualValueArray = new double[actualValue[i][j].length * actualValue[i][j][0].length];
                for (int k = 0; k < inferredValue[i][j].length; k++) {
                    System.arraycopy(inferredValue[i][j][k], 0, inferredValueArray, cumCount, inferredValue[i][j][k].length);
                    System.arraycopy(actualValue[i][j][k], 0, actualValueArray, cumCount, actualValue[i][j][k].length);
                    cumCount+=inferredValue[i][j][k].length;
                }
                callableList.add(() -> new PearsonsCorrelation().correlation(inferredValueArray, actualValueArray));
            }
        }
        List<Double> res = CommandUtils.run_commands(callableList, 32);
        double[][] result = new double[inferredValue.length][inferredValue[0].length];
        for (int i = 0; i < result.length; i++) {
            Arrays.fill(result[i], -1);
        }
        for (int i = 0; i < res.size(); i++) {
            int a = i/inferredValue[0].length;
            int b = i%inferredValue[0].length;
            result[a][b] = res.get(i);
        }
        return result;
    }


    /**
     *
     * @param inferredValue
     * @param actualValue
     * @return contingencyTable, dim1 is different run, dim2 is admixed taxon index, dim is one of
     * [count_truePositive, count_falseNegative, count_falsePositive, count_trueNegative]
     */
    public static int[][][] contingencyTable(double[][][][] inferredValue, double[][][][] actualValue){
        int runNum = inferredValue.length;
        int admixedTaxaNum = inferredValue[0].length;
        int contingencyTableNum = 4;
        int sourceNum = inferredValue[0][0].length;
        int[][][] contingencyTable = new int[runNum][admixedTaxaNum][contingencyTableNum];
        double inferred, actual;
        for (int runIndex = 0; runIndex < runNum; runIndex++) {
            for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaNum; admixedTaxonIndex++) {
                int count_truePositive=0;
                int count_falseNegative=0;
                int count_falsePositive=0;
                int count_trueNegative=0;
                for (int sourceIndex = 0; sourceIndex < sourceNum; sourceIndex++) {
                    for (int variantIndex = 0; variantIndex < inferredValue[runIndex][admixedTaxonIndex][sourceIndex].length; variantIndex++) {
                        inferred = inferredValue[runIndex][admixedTaxonIndex][sourceIndex][variantIndex];
                        actual = actualValue[runIndex][admixedTaxonIndex][sourceIndex][variantIndex];
                        if (actual > 0.5 && inferred > 0.5){
                            count_truePositive++;
                        }else if (actual > 0.5 && inferred < 0.5){
                            count_falseNegative++;
                        }else if (actual < 0.5 && inferred > 0.5){
                            count_falsePositive++;
                        }else if (actual < 0.5 && inferred < 0.5){
                            count_trueNegative++;
                        }
                    }
                }
                contingencyTable[runIndex][admixedTaxonIndex][0] = count_truePositive;
                contingencyTable[runIndex][admixedTaxonIndex][1] = count_falseNegative;
                contingencyTable[runIndex][admixedTaxonIndex][2] = count_falsePositive;
                contingencyTable[runIndex][admixedTaxonIndex][3] = count_trueNegative;
            }
        }
        return contingencyTable;
    }

    /**
     *
     * @param inferredValue
     * @param actualValue
     * @return contingencyTable, dim1 is different run, dim2 is admixed taxon index, dim is one of
     * [count_truePositive, count_falseNegative, count_falsePositive, count_trueNegative]
     */
    public static int[][][] contingencyTable_bitset(BitSet[][][] inferredValue, BitSet[][][] actualValue){
        int runNum = inferredValue.length;
        int admixedTaxaNum = inferredValue[0].length;
        int contingencyTableNum = 4;
        int sourceNum = inferredValue[0][0].length;
        int[][][] contingencyTable = new int[runNum][admixedTaxaNum][contingencyTableNum];
        boolean inferred, actual;
        int variantNum;
        for (int runIndex = 0; runIndex < runNum; runIndex++) {
            variantNum = LocalAncestry.variantNum[runIndex];
            for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaNum; admixedTaxonIndex++) {
                int count_truePositive=0;
                int count_falseNegative=0;
                int count_falsePositive=0;
                int count_trueNegative=0;
                for (int sourceIndex = 0; sourceIndex < sourceNum; sourceIndex++) {
                    for (int variantIndex = 0; variantIndex < variantNum; variantIndex++) {
                        inferred = inferredValue[runIndex][admixedTaxonIndex][sourceIndex].get(variantIndex);
                        actual = actualValue[runIndex][admixedTaxonIndex][sourceIndex].get(variantIndex);
                        if (actual && inferred){
                            count_truePositive++;
                        }else if (actual && (!inferred)){
                            count_falseNegative++;
                        }else if ((!actual) && inferred){
                            count_falsePositive++;
                        }else if ((!actual) && (!inferred)){
                            count_trueNegative++;
                        }
                    }
                }
                contingencyTable[runIndex][admixedTaxonIndex][0] = count_truePositive;
                contingencyTable[runIndex][admixedTaxonIndex][1] = count_falseNegative;
                contingencyTable[runIndex][admixedTaxonIndex][2] = count_falsePositive;
                contingencyTable[runIndex][admixedTaxonIndex][3] = count_trueNegative;
            }
        }
        return contingencyTable;
    }

    /**
     *
     * @param inferredValue
     * @param actualValue
     * @return contingencyTable, dim1 is different run, dim2 is admixed taxon index, dim is one of
     * [count_truePositive, count_falseNegative, count_falsePositive, count_trueNegative]
     */
    public static int[][][] contingencyTable_2way(double[][][][] inferredValue, double[][][][] actualValue){
        int runNum = inferredValue.length;
        int admixedTaxaNum = inferredValue[0].length;
        int contingencyTableNum = 4;
        int[][][] contingencyTable = new int[runNum][admixedTaxaNum][contingencyTableNum];
        double inferred0, inferred1, actual0, actual1;
        int inferredSource, actualSource;
        for (int runIndex = 0; runIndex < runNum; runIndex++) {
            for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaNum; admixedTaxonIndex++) {
                int count_truePositive=0;
                int count_falseNegative=0;
                int count_falsePositive=0;
                int count_trueNegative=0;
                for (int variantIndex = 0; variantIndex < inferredValue[runIndex][admixedTaxonIndex][0].length; variantIndex++) {
                    inferred0 = inferredValue[runIndex][admixedTaxonIndex][0][variantIndex]; // introgressed ancestry
                    inferred1 = inferredValue[runIndex][admixedTaxonIndex][1][variantIndex]; // native ancestry
                    actual0 = actualValue[runIndex][admixedTaxonIndex][0][variantIndex]; // introgressed ancestry
                    actual1 = actualValue[runIndex][admixedTaxonIndex][1][variantIndex]; // native ancestry
                    inferredSource = -1;
                    if (inferred0 > 0.5 && inferred1 < 0.5){
                        inferredSource = 1; // here 1 means introgressed ancestry
                    }else if (inferred0 < 0.5 && inferred1 > 0.5){
                        inferredSource = 0; // here 0 means native ancestry
                    }else {
                        inferredSource = 0; // default ancestry is native
                    }
                    actualSource = actual0 > 0.5 ? 1 : 0;

                    if (actualSource == 1 && inferredSource == 1){
                        count_truePositive++;
                    }else if (actualSource == 1 && inferredSource == 0){
                        count_falseNegative++;
                    }else if (actualSource == 0  && inferredSource == 1){
                        count_falsePositive++;
                    }else if (actualSource == 0 && inferredSource == 0){
                        count_trueNegative++;
                    }
                }
                contingencyTable[runIndex][admixedTaxonIndex][0] = count_truePositive;
                contingencyTable[runIndex][admixedTaxonIndex][1] = count_falseNegative;
                contingencyTable[runIndex][admixedTaxonIndex][2] = count_falsePositive;
                contingencyTable[runIndex][admixedTaxonIndex][3] = count_trueNegative;
            }
        }
        return contingencyTable;
    }

    /**
     *
     * @param simulationMetadata simulationMetadata
     * @param simulationDir simulationDir
     * @return actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population (ordered by
     * introgressed, native), dim4 is variants
     */
    public static double[][][][] extractLocalAncestry_actualValue(SimulationMetadata simulationMetadata,
                                                                  String simulationDir){
        String[] demesID = simulationMetadata.getDemesID();
        double[][][][] actualValue = new double[demesID.length][][][];
        int[] admixedSampleSize = simulationMetadata.getAdmixedPopSampleSize();
        List<String>[] referencePopList = simulationMetadata.getReferencePopList();
        File[] tractFiles = new File[demesID.length];
        File[] genotypeFiles = new File[demesID.length];
        for (int i = 0; i < demesID.length; i++) {
            tractFiles[i] = new File(simulationDir, demesID[i]+".tract");
            genotypeFiles[i] = new File(simulationDir, demesID[i]+".vcf");
            actualValue[i] = extractLocalAncestry_actualValue(genotypeFiles[i], tractFiles[i], admixedSampleSize[i],
                    referencePopList[i]);
        }
        return actualValue;
    }

    /**
     *
     * @param simulationMetadata simulationMetadata
     * @param simulationDir simulationDir
     * @return actualValue dim1 is different run, dim2 is different haplotype, dim3 is source population (ordered by
     * introgressed, native), dim4 is variants
     */
    public static BitSet[][][] extractLocalAncestry_actualValue_bitset(SimulationMetadata simulationMetadata,
                                                                       String simulationDir){
        String[] demesID = simulationMetadata.getDemesID();
        BitSet[][][] actualValue = new BitSet[demesID.length][][];
        int[] admixedSampleSize = simulationMetadata.getAdmixedPopSampleSize();
        List<String>[] referencePopList = simulationMetadata.getReferencePopList();
        File[] tractFiles = new File[demesID.length];
        File[] genotypeFiles = new File[demesID.length];
        GenotypeTable[] genotypeTables = new GenotypeTable[demesID.length];
        LocalAncestry.variantNum = new int[demesID.length];
        for (int i = 0; i < demesID.length; i++) {
            tractFiles[i] = new File(simulationDir, demesID[i]+".tract");
            genotypeFiles[i] = new File(simulationDir, demesID[i]+".vcf");
            genotypeTables[i] = new GenotypeTable(genotypeFiles[i].getAbsolutePath());
            LocalAncestry.variantNum[i] = genotypeTables[i].getSiteNumber();
            actualValue[i] = extractLocalAncestry_actualValue_bitset(genotypeTables[i], tractFiles[i], admixedSampleSize[i],
                    referencePopList[i]);
        }
        return actualValue;
    }

    /**
     *
     * @param genotypFile genotype file
     * @param tractFile simulation tract file
     * @param admixedSampleSize sample size of admixed population
     * @param refPopList introgressed population, native population
     * @return actualValue, dim1 is admixed taxon index, dim2 is different source (ordered by introgressed,
     * native), dim3 is variant index
     */
    private static double[][][] extractLocalAncestry_actualValue(File genotypFile, File tractFile,
                                                                int admixedSampleSize,
                                                                List<String> refPopList){
        GenotypeTable genotypeTable = new GenotypeTable(genotypFile.getAbsolutePath());
        double[][][] localAncestry_actualValue = new double[admixedSampleSize][refPopList.size()][genotypeTable.getSiteNumber()];
        Map<String, Integer> pop2Index =
                IntStream.range(0, refPopList.size()).boxed().collect(Collectors.toMap(refPopList::get, i -> i));

        // tractFile only hava introgressed pop info, so here we set default values of native pop to 1
        // the final ele of refPopList is native pop
        String nativePop = refPopList.get(refPopList.size()-1);
        int nativePopIndex = pop2Index.get(nativePop);
        for (int i = 0; i < localAncestry_actualValue.length; i++) {
            Arrays.fill(localAncestry_actualValue[i][nativePopIndex], 1);
        }
        int refPopIndex;
        try (BufferedReader br = IOTool.getReader(tractFile)) {
            String line;
            List<String> temp;
            br.readLine();
            int haplotypeIndex, startHit, endHit, startIndex, endIndex;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                haplotypeIndex = Integer.parseInt(PStringUtils.fastSplit(temp.get(0), "_").get(1));
                refPopIndex = pop2Index.get(temp.get(1));
                startHit = genotypeTable.getSiteIndex("1", Integer.parseInt(temp.get(2)));
                endHit  = genotypeTable.getSiteIndex("1", Integer.parseInt(temp.get(3)));
                startIndex = startHit < 0 ?  -startHit-1 : startHit;
                endIndex = endHit < 0 ? -endHit-1 : endHit;
                Arrays.fill(localAncestry_actualValue[haplotypeIndex][refPopIndex], startIndex, endIndex, 1);
                Arrays.fill(localAncestry_actualValue[haplotypeIndex][nativePopIndex], startIndex, endIndex, 0);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return localAncestry_actualValue;
    }

    /**
     *
     * @param genotypeTable genotype file
     * @param tractFile simulation tract file
     * @param admixedSampleSize sample size of admixed population
     * @param refPopList introgressed population, native population
     * @return actualValue, dim1 is admixed taxon index, dim2 is different source (ordered by introgressed,
     * native), dim3 is variant index
     */
    private static BitSet[][] extractLocalAncestry_actualValue_bitset(GenotypeTable genotypeTable, File tractFile,
                                                                      int admixedSampleSize,
                                                                      List<String> refPopList){
        BitSet[][] localAncestry_actualValue = new BitSet[admixedSampleSize][refPopList.size()];
        for (int i = 0; i < localAncestry_actualValue.length; i++) {
            for (int j = 0; j < localAncestry_actualValue[i].length; j++) {
                localAncestry_actualValue[i][j] = new BitSet();
            }
        }
        Map<String, Integer> pop2Index =
                IntStream.range(0, refPopList.size()).boxed().collect(Collectors.toMap(refPopList::get, i -> i));

        // tractFile only hava introgressed pop info, so here we set default values of native pop to 1
        // the final ele of refPopList is native pop
        String nativePop = refPopList.get(refPopList.size()-1);
        int nativePopIndex = pop2Index.get(nativePop);
        for (int i = 0; i < localAncestry_actualValue.length; i++) {
            localAncestry_actualValue[i][nativePopIndex].set(0, genotypeTable.getSiteNumber(), true);
        }
        int refPopIndex;
        try (BufferedReader br = IOTool.getReader(tractFile)) {
            String line;
            List<String> temp;
            br.readLine();
            int haplotypeIndex, startHit, endHit, startIndex, endIndex;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                haplotypeIndex = Integer.parseInt(PStringUtils.fastSplit(temp.get(0), "_").get(1));
                refPopIndex = pop2Index.get(temp.get(1));
                startHit = genotypeTable.getSiteIndex("1", Integer.parseInt(temp.get(2)));
                endHit  = genotypeTable.getSiteIndex("1", Integer.parseInt(temp.get(3)));
                startIndex = startHit < 0 ?  -startHit-1 : startHit;
                endIndex = endHit < 0 ? -endHit-1 : endHit;
                localAncestry_actualValue[haplotypeIndex][refPopIndex].set(startIndex, endIndex, true);
                localAncestry_actualValue[haplotypeIndex][nativePopIndex].set(startIndex, endIndex, false);
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return localAncestry_actualValue;
    }

    public static void write_RobustnessData(String simulationMetadataFile, String simulationDir,
                                            List<double[][][][]> software_localAncestry,
                                            List<String> softwareList, String outFileRobustnessData){
        SimulationMetadata simulationMetadata = new SimulationMetadata(simulationMetadataFile);
        double[][][][] actual_values = LocalAncestry.extractLocalAncestry_actualValue(simulationMetadata,
                simulationDir);
        double[][][] meanDeviationData = new double[software_localAncestry.size()][][];
        double[][][] pearsonCorrelationData = new double[software_localAncestry.size()][][];
        int[][][][] contingencyTableData = new int[software_localAncestry.size()][][][];
        String[] software = new String[software_localAncestry.size()];
        for (int i = 0; i < software_localAncestry.size(); i++) {
            meanDeviationData[i] = LocalAncestry.meanDeviation(software_localAncestry.get(i), actual_values);
            pearsonCorrelationData[i] = LocalAncestry.pearsonCorrelation(software_localAncestry.get(i), actual_values);
            contingencyTableData[i] = LocalAncestry.contingencyTable(software_localAncestry.get(i), actual_values);
            software[i] = softwareList.get(i);
        }
        LocalAncestry.write_RobustnessData(meanDeviationData, pearsonCorrelationData, contingencyTableData, software,
                simulationMetadata,
                outFileRobustnessData);
    }

    /**
     *
     * @param meanDeviationData dim1 is different software, dim2 is different run, dim3 is different haplotype
     * @param pearsonCorrelationData dim1 is different software, dim2 is different run, dim3 is different haplotype
     * @param software
     * @param simulationMetadata
     * @param outFile_evaluation
     */
    private static void write_RobustnessData(double[][][] meanDeviationData, double[][][] pearsonCorrelationData,
                                             int[][][][] contingencyTableData,
                                             String[] software, SimulationMetadata simulationMetadata,
                                             String outFile_evaluation){

        try (BufferedWriter bw = IOTool.getWriter(outFile_evaluation)) {
            bw.write("DemesID\tSoftware\tAdmixedIndividual\tMeanDeviation\tPearsonCorrelation\tTruePositive" +
                    "\tFalseNegative\tFalsePositive\tTrueNegative");

            bw.newLine();
            StringBuilder sb = new StringBuilder();
            String joined;
            for (int i = 0; i < meanDeviationData.length; i++) {
                for (int j = 0; j < meanDeviationData[i].length; j++) {
                    for (int k = 0; k < meanDeviationData[i][j].length; k++) {
                        sb.setLength(0);
                        sb.append(simulationMetadata.getDemesID()[j]).append("\t");
                        sb.append(software[i]).append("\t");
                        sb.append("tsk_").append(k).append("\t");
                        sb.append(meanDeviationData[i][j][k]).append("\t");
                        sb.append(pearsonCorrelationData[i][j][k]).append("\t");
                        joined = Joiner.on("\t").join(Ints.asList(contingencyTableData[i][j][k]));
                        sb.append(joined);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
