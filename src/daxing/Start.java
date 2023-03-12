package daxing;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.demography.DemographicModelTools;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
import daxing.v2.localAncestryInfer.runner.*;
import daxing.v2.localAncestryInfer.simulation.Simulation;
import daxing.v2.localAncestryInfer.simulation.SimulationMetadata;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.BitSet;
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
//        int maxSolutionCount=32; // 100
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
////
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
//        String genotypeFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/003_simulation/D014.vcf";
//        int windowSize = 200;
//        int stepSize = 100;
//        String taxaGroupFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/004_runner/laidp/D012/D012.taxaGroup.txt";
//        String ancestryAllele = "simulation";
//        int conjunctionNum = 2;
//        double switchCostScore = 1.5;
//        String localAnceOutFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/004_runner/laidp/D012" +
//                "/D012.localAnc_2.txt";
//        int threadsNum = 2;
//////
//////
//////        String outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp";
//////////
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
////        long start = System.nanoTime();
//        GenotypeTable.run_LAIDP(genotypeFile, windowSize, stepSize, taxaGroupFile, ancestryAllele, conjunctionNum,
//                switchCostScore, localAnceOutFile, threadsNum);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s");

//        MD5.checkTwoFileMD5("/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/005_evaluation/evaluation_contingencyTable.txt","/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/005_evaluation/evaluation.txt");




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

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_2";
//        String outDir2 = "/Users/xudaxing/Desktop/LAIDP_development/threeWay";
//        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.THREE_WAY,outDir2);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.FOUR_WAY, outDir3);

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M";
//////////////////////////////////////////////////////////////////////////////////////        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
//        evaluate_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);
//        evaluate_contingencyTable(DemographicModelTools.N_way.THREE_WAY, outDir);
//        System.out.println();

//        MD5.checkTwoFileMD5("/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/005_evaluation/evaluation_contingencyTable_2way_bitset.txt",
//                "/Users/xudaxing/Desktop/LAIDP_development/twoWay_ancient/005_evaluation/evaluation.txt");

//        int[] gridSource_fragment = {1,1,2,2,2,2,2,1,1};
//        int n_wayAdmixture = 3;
//        int windowSize = 5;
//        double[][]  res = GenotypeTable.get_state_trans_prob2(gridSource_fragment, n_wayAdmixture, windowSize);
//        System.out.println();
//
//        GenotypeTable genotypeTable = new GenotypeTable("/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_generation100/002_simulatedGenotype/simulate_proportion0.1_generation100_addAncestral.vcf");
//        int windowSize = 200;
//        int stepSize = 100;
//        String taxaGroupFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/taxaGroup.txt";
//        BitSet[] ancestralAlleleBitSet = genotypeTable.getAncestralAlleleBitSet("simulation");
//        int threadsNum = 2;
//        String fdOutDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_generation100/004_fd/007_fd_threshold";
//        String individualLocalAncestryDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_generation100/004_fd/008_localAnc";
//        genotypeTable.calculateLocalAncestry_check_fd_gridSource(windowSize, stepSize, taxaGroupFile,
//                ancestralAlleleBitSet, threadsNum, fdOutDir, individualLocalAncestryDir);

//        String file1 = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/003_twoWay_proportion0.1_generation100/004_fd/004_individualLocalAncestry/IndividualFd_tsk_61simulated.txt";
//        String file2 = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/003_twoWay_proportion0.1_generation100/004_fd/temp/SimulatedLocalAncestry_tsk_61simulated.txt";
//        MD5.checkTwoFileMD5(file1,file2);


//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype = {1,1,0,0,0,1,0,0,1,1,1,1};

//        double[][] alts = {{-1, -2, -2, -3, -4, -1, -2, -5},
//                           {-3, -2, -1, -1, -3, -4, -6, -1},
//                           {-2,-3, -2, -4, -4, -4, -5, -6}};
//        double[][] state_tran = {{-0.1, -3, -3},
//                                 {-3,-0.1,-3},
//                                 {-3,-3,-0.1}};
//
//        double[] start_prob = {-0.2, -2, -2};
//
////        double[][] alts = {{0.1, 0.01, 0.4, 0.4, 0.2, 0.4, 0.5, 0.2},
////                {0.1, 0.2, 0.3, 0.1, 0.2, 0.4, 0.2, 0.4},
////                {0.2,0.1, 0.3, 0.1, 0.3, 0.1, 0.02, 0.01}};
////        double[][] state_tran = {{0.5, 0.25, 0.25},
////                {0.25,0.5,0.25},
////                {0.25,0.25,0.5}};
////
////        double[] start_prob = {0.3, 0.3, 0.4};
//        int[] obs = {1,1,0,0,1,1,0,1};
//
//        int[][] path = HMM.hmm(alts, state_tran, start_prob, obs);
//        int[] path1 = HMM.hmm2(alts, state_tran, start_prob, obs);
//        System.out.println();

//        List<String> stcList = new ArrayList<>();
//        stcList.add("1");
//        stcList.add("2");
//        stcList.add("3");
//        stcList.add("4");
//        stcList.add("5");
//        stcList.add("6");
//        Map<String, Source> map = new HashMap<>();
//        map.put("1", Source.NATIVE);
//        map.put("2", Source.NATIVE);
//        map.put("3", Source.INTROGRESSED_1);
//        map.put("4", Source.INTROGRESSED_1);
//        map.put("5", Source.INTROGRESSED_2);
//        map.put("6", Source.INTROGRESSED_2);
//
//
////        long start0 = System.nanoTime();
////        for (int i = 0; i < 10_000; i++) {
////            Solution.getCandidateSolution(srcGenotype, queryGenotype, 1.5, stcList, map);
////        }
////        long start1 = System.nanoTime();
////        for (int i = 0; i < 10_000; i++) {
////            Solution.getCandidateSolution2(srcGenotype, queryGenotype, 1.5);
////        }
////        long start2 = System.nanoTime();
////        System.out.println((double) (start1-start0)/(start2-start1));
////
////        System.out.println(start1-start0);
////        System.out.println(start2-start1);
//
//        IntSet[] res1 = Solution.getCandidateSolution2(srcGenotype, queryGenotype, 1.5);
//        System.out.println();


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