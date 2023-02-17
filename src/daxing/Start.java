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
////
//////////
//////////////
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



//        String genotypeTableFile = "/Users/xudaxing/Desktop/threeWay/003_simulation/D003.vcf";
//        int windowSize = 200;
//        int stepSize = 100;
//        int threadsNum = 3;
//        GenotypeTable genotypeTable = new GenotypeTable(genotypeTableFile);
//        int variantsNum = genotypeTable.getSiteNumber();
//
//        int[] pop_taxonIndex_admixed = IntStream.range(0, 30).toArray();
//        int[] pop_taxonIndex_native = IntStream.range(30,60).toArray();
//        int[][] pop_taxonIndex_introgressed = new int[2][];
//        pop_taxonIndex_introgressed[0] = IntStream.range(60,90).toArray();
//        pop_taxonIndex_introgressed[1] = IntStream.range(90,120).toArray();
//
//        int[] windowStartIndexArray = GenotypeTable.getWindowStartIndex(stepSize, genotypeTable.getSiteNumber());
//        double[][][] dxy_windows_admixed = genotypeTable.calculateAdmixedDxy(pop_taxonIndex_admixed, pop_taxonIndex_native,
//                pop_taxonIndex_introgressed, windowStartIndexArray, windowSize);
//        double[][] dxy_pairwise_nativeIntrogressed = genotypeTable.calculatePairwiseDxy(pop_taxonIndex_admixed,
//                pop_taxonIndex_native,pop_taxonIndex_introgressed,
//                windowStartIndexArray, windowSize);
//
//
//        BitSet[] ancestralAlleleBitSet = new BitSet[2];
//        ancestralAlleleBitSet[0] = new BitSet();
//        ancestralAlleleBitSet[1] = new BitSet();
//
//        int[][] popTaxaIndices = new int[1+30+2][];
//        popTaxaIndices[0] = IntStream.range(30,60).toArray();
//        int k=0;
//        for (int i = 1; i < popTaxaIndices.length-2; i++) {
//            popTaxaIndices[i] = new int[1];
//            popTaxaIndices[i][0] = k;
//            k++;
//        }
//        popTaxaIndices[31] = IntStream.range(60,90).toArray();
//        popTaxaIndices[32] = IntStream.range(90,120).toArray();
//
//        double[][] dafs = genotypeTable.calculateDaf(threadsNum, popTaxaIndices, ancestralAlleleBitSet);
//
//        double[] dafs_native = dafs[0];
//        double[][] dafs_admixed = new double[30][variantsNum];
//        double[][] dafs_introgressed = new double[2][variantsNum];
//
//        for (int i = 0; i < dafs.length; i++) {
//            if (i == 0) continue;
//            if (i > 30){
//                dafs_introgressed[i-31] = dafs[i];
//            }else {
//                dafs_admixed[i-1] = dafs[i];
//            }
//
//        }
//
//        double[][][] fd = GenotypeTable.calculate_fd(threadsNum, dafs_native, dafs_admixed, dafs_introgressed,
//                windowStartIndexArray,
//                windowSize, variantsNum);
//        int[][] gridSource = GenotypeTable.calculateSource(fd, dxy_pairwise_nativeIntrogressed, dxy_windows_admixed);
//
//        int conjunctionNum = 2;
//        int[] sourceTaxIndices = IntStream.range(30,120).toArray();
//        String taxaGroupFile = "/Users/xudaxing/Desktop/taxaGroup.txt";
//        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);
//        long start = System.nanoTime();
//        double switchCostScore = 1.5;
////        BitSet[][] localAnc = genotypeTable.calculateLocalAncestry(windowSize, stepSize, taxaGroupFile, ancestralAlleleBitSet,
////                conjunctionNum, switchCostScore, threadsNum);
//        BitSet[][] localAnc = GenotypeTable.calculateLocalAncestry(windowStartIndexArray, windowSize, variantsNum,
//                gridSource, conjunctionNum, genotypeTable, pop_taxonIndex_admixed, sourceTaxIndices, 3, 1.5,
//                taxaGroup, threadsNum);
//        System.out.println("completed in "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
//        String outFile = "/Users/xudaxing/Desktop/localAncestry.txt";
//        GenotypeTable.write_localAncestry(localAnc, outFile, variantsNum);
//        System.out.println();


//        String genotypeFile = "/Users/xudaxing/Desktop/threeWay/003_simulation/D003.vcf";
//        int windowSize = 200;
//        int stepSize = 100;
//        String taxaGroupFile = "/Users/xudaxing/Desktop/taxaGroup.txt";
//        String ancestryAllele = "simulation";
//        int conjunctionNum = 2;
//        double switchCostScore = 1.5;
//        String localAnceOutFile = "/Users/xudaxing/Desktop/localAncestry.txt";
//        int threadsNum = 2;

//        String genotypeFile = args[0];
//        int windowSize =Integer.parseInt(args[1]);
//        int stepSize = Integer.parseInt(args[2]);
//        String taxaGroupFile = args[3];
//        String ancestryAllele = args[4];
//        int conjunctionNum = Integer.parseInt(args[5]);
//        double switchCostScore = Double.parseDouble(args[6]);
//        String localAnceOutFile = args[7];
//        int threadsNum = Integer.parseInt(args[8]);
//        GenotypeTable.run_LAIDP(genotypeFile, windowSize, stepSize, taxaGroupFile, ancestryAllele, conjunctionNum,
//                switchCostScore, localAnceOutFile, threadsNum);


//        double switchCostScore= 1.5;
//        int[][] srcGenotype = {{0,1,0,1,0,1,0,0,0,0,1,1},
//                            {0,0,0,1,0,1,1,0,0,0,1,1},
//                            {0,0,1,0,1,0,0,0,1,0,1,1},
//                            {0,0,0,0,1,0,1,0,1,1,1,1},
//                            {1,1,0,0,0,0,1,1,1,1,0,0},
//                            {1,0,0,1,0,0,1,1,1,1,0,0}};
//        int[] queryGenotype =       {1,1,0,0,0,1,0,0,1,1,1,1};


//        String src= "0,1,2,3,4,5";
//        List<String> srcIndiList = PStringUtils.fastSplit(src,",");
//
//        Map<String, Source> taxaSourceMap = new HashMap<>();
//        taxaSourceMap.put("0", Source.NATIVE_SOURCE_0);
//        taxaSourceMap.put("1", Source.NATIVE_SOURCE_0);
//        taxaSourceMap.put("2", Source.INTROGRESSED_SOURCE_1);
//        taxaSourceMap.put("3", Source.INTROGRESSED_SOURCE_1);
//        taxaSourceMap.put("4", Source.INTROGRESSED_SOURCE_2);
//        taxaSourceMap.put("5", Source.INTROGRESSED_SOURCE_2);
//
//        double[][] miniCost = Solution.getMiniCostScore(srcGenotype, queryGenotype, switchCostScore);
//        IntList[] solution = Solution.getCandidateSolution(srcGenotype, queryGenotype, switchCostScore, srcIndiList,
//                taxaSourceMap);
//        System.out.println();
//        long start = System.nanoTime();
//        for (int i = 0; i < 1_000_000; i++) {
//            IntList[] solution = Solution.getCandidateSolution(srcGenotype, queryGenotype, switchCostScore, srcIndiList,
//                    taxaSourceMap);
//        }
//        long start1 = System.nanoTime();
//        for (int i = 0; i < 1_000_000; i++) {
//            IntList[] solution2 = Solution.getCandidateSolution2(srcGenotype, queryGenotype, switchCostScore, srcIndiList,
//                    taxaSourceMap);
//        }
//        long start2 = System.nanoTime();
//        long duration1 = start1-start;
//        long duration2 = start2-start1;
//        System.out.println(duration1);
//        System.out.println(duration2);
//        System.out.println((double) duration2/duration1);
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


//        String modelOutPaht = "/Users/xudaxing/Desktop/test.yaml";
//        int[] splitEventTime = {10000, 9000, 1000};
//        double[] admixtureProportion = {0.1, 0.2};
//        int[] admixtureTime = {100, 50};
//        DemographicModel demographicModel = DemographicModelTools.threeWayBuilder_equilibriumPopulation(splitEventTime,
//                admixtureProportion, admixtureTime);
////        DemographicModel demographicModel = DemographicModelTools.threeWayBuilder_equilibriumPopulation(splitEventTime,
////                admixtureProportion, admixtureTime);
//        DemographicModelTools.writeModel(demographicModel, new File(modelOutPaht));

        String outDir = "/Users/xudaxing/Desktop/twoWay";
//        String outDir2 = "/Users/xudaxing/Desktop/threeWay";
//        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
        run(DemographicModelTools.N_way.TWO_WAY, outDir);
//        run(DemographicModelTools.N_way.THREE_WAY,outDir2);
//        run(DemographicModelTools.N_way.FOUR_WAY, outDir3);

//        String a = "/Users/xudaxing/Desktop/multipleWay/simulationMetadata.txt";
//        String b = "/Users/xudaxing/Desktop/multipleWay";
//        DemographicModelTools.batchRun_multipleWay(a, b);

//        int[] splitEventTime = {10000, 9000, 500};
//        double[] admixtureProportion = {0.1, 0.2};
//        int[] admixtureTime = {500, 100};
//        DemographicModel demographicModel =DemographicModelTools.threeWayBuilder_equilibriumPopulation(splitEventTime, admixtureProportion, admixtureTime);
//        DemographicModelTools.writeModel(demographicModel, new File("/Users/xudaxing/Desktop/multipleWay/test.yaml"));

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
                default:
                    break;
            }
        });

        Robustness.write_RobustnessData(simulationMetadataOutFile, dirsFile[2].getAbsolutePath(),
                software_localAncestry, software, new File(dirsFile[5], "evaluation.txt").getAbsolutePath());

    }




}