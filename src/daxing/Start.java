package daxing;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.utiles.IOTool;
import daxing.common.utiles.MD5;
import daxing.v2.localAncestryInfer.demography.DemographicModelTools;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
import daxing.v2.localAncestryInfer.laidp.GenotypeTable;
import daxing.v2.localAncestryInfer.runner.*;
import daxing.v2.localAncestryInfer.simulation.Simulation;
import daxing.v2.localAncestryInfer.simulation.SimulationMetadata;
import pgl.infra.utils.Benchmark;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.List;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) throws IOException, InterruptedException {
        String genotypeFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/003_simulation/D014.vcf";
        int windowSize = 200;
        int stepSize = 100;
        String taxaGroupFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.taxaGroup.txt";
        String ancestryAllele = "simulation";
        int conjunctionNum = 2;
        double switchCostScore = 1.5;
        int maxSolutionCount = 32;
        String localAnceOutFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014" +
                ".localAnc_2.txt";
        int threadsNum = 2;
//////
//////
//////        String outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp";

//        String genotypeFile = args[0];
//        int windowSize =Integer.parseInt(args[1]);
//        int stepSize = Integer.parseInt(args[2]);
//        String taxaGroupFile = args[3];
//        String ancestryAllele = args[4];
//        int conjunctionNum = Integer.parseInt(args[5]);
//        double switchCostScore = Double.parseDouble(args[6]);
//        String localAnceOutFile = args[7];
//        int threadsNum = Integer.parseInt(args[8]);
//////////////
        long start = System.nanoTime();
        GenotypeTable.run_LAIDP(genotypeFile, windowSize, stepSize, taxaGroupFile, ancestryAllele, conjunctionNum,
                switchCostScore, maxSolutionCount, localAnceOutFile, threadsNum);

        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s");

        MD5.checkTwoFileMD5("/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.localAnc.txt","/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.localAnc_2.txt");
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s ");

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_2";
//        String outDir2 = "/Users/xudaxing/Desktop/LAIDP_development/threeWay";
//        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.THREE_WAY,outDir2);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.FOUR_WAY, outDir3);

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M";
//        evaluate_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);

    }

    /**
     * pearsonCorrelationData meanDeviationData contingencyTableData
     * @param nWay
     * @param outDir
     */
    public static void evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way nWay, String outDir){

        final String[] DIRS = {"001_parameterFile","002_demes","003_simulation","004_runner","log",
                "005_evaluation"};
        final String[] SOFTWARE = {"loter","elai","mosaic","laidp"};

        File[] dirsFile = new File[DIRS.length];
        for (int i = 0; i < dirsFile.length; i++) {
            dirsFile[i] = new File(outDir, DIRS[i]);
            dirsFile[i].mkdir();
        }

        String[] logFiles = new String[SOFTWARE.length];
        String[] softwareSubDir = new String[SOFTWARE.length];
        for (int i = 0; i < SOFTWARE.length; i++) {
            File file = new File(dirsFile[3], SOFTWARE[i]);
            file.mkdir();
            softwareSubDir[i] = file.getAbsolutePath();
            logFiles[i] = new File(dirsFile[4], SOFTWARE[i]+".log").getAbsolutePath();
        }

        String simulationMetadataOutFile = new File(dirsFile[0], "simulationMetadata.txt").getAbsolutePath();
        String simulationLogFile = new File(dirsFile[4], "simulation.log").getAbsolutePath();
        DemographicModelTools.batchRun(nWay, simulationMetadataOutFile, dirsFile[1].getAbsolutePath());
        Simulation simulation = new Simulation.Builder(simulationMetadataOutFile, simulationLogFile,
                dirsFile[2].getAbsolutePath()).build();
        simulation.run_simulation();

        GenotypeMetaData genotypeMetaData = new GenotypeMetaData(simulationMetadataOutFile, dirsFile[2].getAbsolutePath());

        List<double[][][][]> software_localAncestry = new ArrayList<>();
        List<String> software = new ArrayList<>();

//         使用 Java 8 Streams API
        IntStream.range(0, SOFTWARE.length).forEach(i -> {
            switch(SOFTWARE[i]) {
                case "loter":
                    Loter_runner loterRunner = new Loter_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    loterRunner.startRun();
                    software_localAncestry.add(loterRunner.extractLocalAncestry());
                    software.add(SOFTWARE[i]);
                    break;
                case "elai":
                    ELAI_runner elaiRunner = new ELAI_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    elaiRunner.startRun();
                    software_localAncestry.add(elaiRunner.extractLocalAncestry());
                    software.add(SOFTWARE[i]);
                    break;
                case "mosaic":
                    Mosaic_runner mosaicRunner = new Mosaic_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    mosaicRunner.startRun();
                    software_localAncestry.add(mosaicRunner.extractLocalAncestry());
                    software.add(SOFTWARE[i]);
                    break;
                case "laidp":
                    LAIDP_runner laidpRunner = new LAIDP_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    laidpRunner.startRun();
                    software_localAncestry.add(laidpRunner.extractLocalAncestry());
                    software.add(SOFTWARE[i]);
                default:
                    break;
            }
        });

        LocalAncestry.write_RobustnessData(simulationMetadataOutFile, dirsFile[2].getAbsolutePath(),
                software_localAncestry, software, new File(dirsFile[5], "evaluation.txt").getAbsolutePath());
    }

    public static void evaluate_contingencyTable(DemographicModelTools.N_way nWay, String outDir){

        final String[] DIRS = {"001_parameterFile","002_demes","003_simulation","004_runner","log",
                "005_evaluation"};
        final String[] SOFTWARE = {"loter","elai","laidp"};

        File[] dirsFile = new File[DIRS.length];
        for (int i = 0; i < dirsFile.length; i++) {
            dirsFile[i] = new File(outDir, DIRS[i]);
            dirsFile[i].mkdir();
        }

        String[] logFiles = new String[SOFTWARE.length];
        String[] softwareSubDir = new String[SOFTWARE.length];
        for (int i = 0; i < SOFTWARE.length; i++) {
            File file = new File(dirsFile[3], SOFTWARE[i]);
            file.mkdir();
            softwareSubDir[i] = file.getAbsolutePath();
            logFiles[i] = new File(dirsFile[4], SOFTWARE[i]+".log").getAbsolutePath();
        }

        String simulationMetadataOutFile = new File(dirsFile[0], "simulationMetadata.txt").getAbsolutePath();

        String simulationLogFile = new File(dirsFile[4], "simulation.log").getAbsolutePath();
        DemographicModelTools.batchRun(nWay, simulationMetadataOutFile, dirsFile[1].getAbsolutePath());
        Simulation simulation = new Simulation.Builder(simulationMetadataOutFile, simulationLogFile,
                dirsFile[2].getAbsolutePath()).build();
        simulation.run_simulation();

        SimulationMetadata simulationMetadata = new SimulationMetadata(simulationMetadataOutFile);
//        double[][][][] actual_values = LocalAncestry.extractLocalAncestry_actualValue(simulationMetadata, dirsFile[2].getAbsolutePath());
        BitSet[][][] actual_values = LocalAncestry.extractLocalAncestry_actualValue_bitset(simulationMetadata,
                dirsFile[2].getAbsolutePath());

        GenotypeMetaData genotypeMetaData = new GenotypeMetaData(simulationMetadataOutFile, dirsFile[2].getAbsolutePath());
        int[][][][] software_contingencyTable = new int[SOFTWARE.length][][][];


//         使用 Java 8 Streams API
        IntStream.range(0, SOFTWARE.length).forEach(i -> {
            switch(SOFTWARE[i]) {
                case "loter":
                    Loter_runner loterRunner = new Loter_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    loterRunner.startRun();
                    software_contingencyTable[i]=loterRunner.contingencyTable_2way_bitset(actual_values);
                    break;
                case "elai":
                    ELAI_runner elaiRunner = new ELAI_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    elaiRunner.startRun();
                    software_contingencyTable[i]=elaiRunner.contingencyTable_2way_bitset(actual_values);
                    break;
                case "mosaic":
                    Mosaic_runner mosaicRunner = new Mosaic_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    mosaicRunner.startRun();
                    software_contingencyTable[i] = mosaicRunner.contingencyTable_2way_bitset(actual_values);
                    break;
                case "laidp":
                    LAIDP_runner laidpRunner = new LAIDP_runner.Builder(genotypeMetaData, logFiles[i], softwareSubDir[i]).build();
                    laidpRunner.startRun();
                    software_contingencyTable[i] = laidpRunner.contingencyTable_2way_bitset(actual_values);
                default:
                    break;
            }
        });

        try (BufferedWriter bw = IOTool.getWriter(new File(dirsFile[5], "evaluation.txt").getAbsolutePath())) {
            bw.write("DemesID\tSoftware\tAdmixedIndividual\tTruePositive" +
                    "\tFalseNegative\tFalsePositive\tTrueNegative");

            bw.newLine();
            StringBuilder sb = new StringBuilder();
            String joined;
            for (int i = 0; i < software_contingencyTable.length; i++) {
                for (int j = 0; j < software_contingencyTable[i].length; j++) {
                    for (int k = 0; k < software_contingencyTable[i][j].length; k++) {
                        sb.setLength(0);
                        sb.append(simulationMetadata.getDemesID()[j]).append("\t");
                        sb.append(SOFTWARE[i]).append("\t");
                        sb.append("tsk_").append(k).append("\t");
                        joined = Joiner.on("\t").join(Ints.asList(software_contingencyTable[i][j][k]));
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