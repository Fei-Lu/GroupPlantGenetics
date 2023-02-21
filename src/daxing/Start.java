package daxing;

import daxing.v2.localAncestryInfer.demography.DemographicModelTools;
import daxing.v2.localAncestryInfer.evaluation.Robustness;
import daxing.v2.localAncestryInfer.runner.*;
import daxing.v2.localAncestryInfer.simulation.Simulation;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) throws IOException, InterruptedException {

//        String genotypeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/001_simulatedGenotype/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_ancestral.vcf";
//        String fd_dxyFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/003_fd/003_individualLocalAncestry/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
//        String groupByPop2IndividualFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/simulated_group.txt";
////        String outDir_deme="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
//        String outDir_deme="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp";
//
//        int conjunctionNum=2; // 2
//        double switchCostScore=1.5; // 1.5
//        int maxSolutionCount=20; // 100
//
////        String genotypeFile = args[0];
////        String fd_dxyFile = args[1];
////        String groupByPop2IndividualFile=args[2];
////        String outDir_deme=args[3];
////
////        int conjunctionNum=Integer.parseInt(args[4]); // 2
////        double switchCostScore=Double.parseDouble(args[5]); // 1.5
////        int maxSolutionCount=Integer.parseInt(args[6]); // 100
//////        int threadNum= Integer.parseInt(args[7]);
//
//
//        long start = System.nanoTime();
//        LocalAncestryInferenceStart.inferLocalAncestry("1", new File(genotypeFile),
//                new File(groupByPop2IndividualFile), new File(fd_dxyFile), conjunctionNum, switchCostScore,
//                new File(outDir_deme), maxSolutionCount);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds");

//        String file1 = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/005_twoWay/004_laidpRes/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne/chr1_tsk_61_LAI.txt";
//        String file2 = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/005_twoWay/004_laidpRes/temp/chr1_tsk_61_LAI.txt";
//        MD5.checkTwoFileMD5(file1,file2);


////
//        String genotypeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/001_simulatedGenotype/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_ancestral.vcf";
//        int windowSize = 200;
//        int stepSize = 100;
//        String taxaGroupFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/taxaGroup.txt";
//        String ancestryAllele = "simulation";
//        int conjunctionNum = 2;
//        double switchCostScore = 1.5;
//        String localAnceOutFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp/localAncestry.txt";
//        int threadsNum = 2;
//
//
////        String outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp";
//////
//        String genotypeFile = args[0];
//        int windowSize =Integer.parseInt(args[1]);
//        int stepSize = Integer.parseInt(args[2]);
//        String taxaGroupFile = args[3];
//        String ancestryAllele = args[4];
//        int conjunctionNum = Integer.parseInt(args[5]);
//        double switchCostScore = Double.parseDouble(args[6]);
//        String localAnceOutFile = args[7];
//        int threadsNum = Integer.parseInt(args[8]);
////        long start = System.nanoTime();
//        GenotypeTable.run_LAIDP(genotypeFile, windowSize, stepSize, taxaGroupFile, ancestryAllele, conjunctionNum,
//                switchCostScore, localAnceOutFile, threadsNum);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s");

//        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
//        BitSet[] ancestralAlleleBitSet = genotypeTable.getAncestralAlleleBitSet(ancestryAllele);
//        genotypeTable.calculateLocalAncestry_test(windowSize, stepSize, taxaGroupFile,
//                ancestralAlleleBitSet, conjunctionNum, switchCostScore, threadsNum, outDir);


////
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s ");

//
//        String simulatedTractDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/002_tracts/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
//        String genotypeTableFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/001_simulatedGenotype/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_ancestral.vcf";
//        String laidpResDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
//        String loterFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/005_loterFile/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_LAI.txt";
//        String outFile_accuracy = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/005_twoWay/006_evaluation/contingencyTable_proportion0.1_genetration100_P1P2_divergence0.1Ne.txt";
////        Evaluation.write_accuracy_recall_precision(simulatedTractDir, genotypeTableFile, laidpResDir, outFile_accuracy);
//        Evaluation.write_ContingencyTable(simulatedTractDir, genotypeTableFile, laidpResDir, loterFile,
//                outFile_accuracy);



        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_temp";
////////        String outDir2 = "/Users/xudaxing/Desktop/threeWay";
////////        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
        run(DemographicModelTools.N_way.TWO_WAY, outDir);
//        run(DemographicModelTools.N_way.THREE_WAY,outDir2);
//        run(DemographicModelTools.N_way.FOUR_WAY, outDir3);

    }

    public static void run(DemographicModelTools.N_way nWay, String outDir){

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

        Robustness.write_RobustnessData(simulationMetadataOutFile, dirsFile[2].getAbsolutePath(),
                software_localAncestry, software, new File(dirsFile[5], "evaluation.txt").getAbsolutePath());

    }




}