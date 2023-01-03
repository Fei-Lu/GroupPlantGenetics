package daxing;

import daxing.v2.mosaic.Mosaic_utils;

import java.io.IOException;

public class Start {

    public static void main(String[] args) throws IOException, InterruptedException {

//        String genotypeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/001_simulatedGenotype/simulate_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne_ancestral.vcf";
//        String fd_dxyFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/003_fd/003_individualLocalAncestry/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
//        String groupByPop2IndividualFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/simulated_group.txt";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/06_two_way_admixture_proportion0.1_genetration100_P1P2_divergence0.1Ne";
////////
////////////
//        int conjunctionNum=2; // 2
//        double switchCostScore=1.5; // 1.5
//        int maxSolutionCount=20; // 100
////
////        String genotypeFile = args[0];
////        String fd_dxyFile = args[1];
////        String groupByPop2IndividualFile=args[2];
////        String outDir=args[3];
////
////        int conjunctionNum=Integer.parseInt(args[4]); // 2
////        double switchCostScore=Double.parseDouble(args[5]); // 1.5
////        int maxSolutionCount=Integer.parseInt(args[6]); // 100
//////        int threadNum= Integer.parseInt(args[7]);
////
////
//        long start = System.nanoTime();
//        LocalAncestryInferenceStart.inferLocalAncestry("1", new File(genotypeFile),
//                new File(groupByPop2IndividualFile), new File(fd_dxyFile), conjunctionNum, switchCostScore,
//                new File(outDir), maxSolutionCount);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds");

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

//        String simulatedTractDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/002_twoWay_proportion0.1_generation500/003_simulatedTract";
//        String genotypeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/002_twoWay_proportion0.1_generation500/002_simulatedGenotype/simulate_proportion0.1_generation500.vcf";
//        String laidpResDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/002_twoWay_proportion0.1_generation500/005_laidpRes";
//        String loterResFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/002_twoWay_proportion0.1_generation500/006_loterFile/simulate_proportion0.1_generation500_LAI.txt";
//        String contingencyTableOutFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/002_twoWay_proportion0.1_generation500/007_evaluation/contingencyTable.txt";
//        Evaluation.write_ContingencyTable(simulatedTractDir, genotypeFile, laidpResDir, loterResFile, contingencyTableOutFile);
//        System.out.println();
//
//        Utils.addAncestral("/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/004_twoWay_proportion0.05_generation100/002_simulatedGenotype/simulate.vcf",
//                "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/004_twoWay_proportion0.05_generation100/002_simulatedGenotype/simulate_ancestral.vcf");

//        String pythonInterpreterPath = "/Users/xudaxing/anaconda3/envs/ABBA_BABA/bin/python";
//        String parseVCFPy = "/Users/xudaxing/Software/ABBA_BABA/genomics_general/VCF_processing/parseVCF.py";
//        String simulatedGenotypeWithAncestral = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_generation100/002_simulatedGenotype/simulate_proportion0.1_generation100_addAncestral.vcf";
//        String ploidyFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_generation100/004_fd/001_group/ploidy.txt";
//        String simulatedGenotypeWithAncestralGeno_outFile = "/Users/xudaxing/Documents/deleteriousMutation" +
//                "/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0" +
//                ".1_generation100/temp/simulate_proportion0.1_generation100_addAncestral.geno";
//        Simulation.transform_vcf2Geno(pythonInterpreterPath, parseVCFPy, simulatedGenotypeWithAncestral, ploidyFile,
//                simulatedGenotypeWithAncestralGeno_outFile);

//        String pythonInterpreterPath = "/Users/xudaxing/anaconda3/envs/Msprime/bin/python";
//        String laidp_simulatePy = "/Users/xudaxing/IdeaProjects/PlantGenetics/GroupPlantGenetics/src/daxing/v2/localAncestryInfer/LAIDP.py";
//        String graphFileDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/000_graphs";
//        int sequence_length = 100_000_000;
//        int random_seed = 1;
//        double mutation_rate = 7.1e-9;
//        String workingDirectory = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay";
//        String logDir= "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/log";
//        int threadsNum = 10;
//        String pythonInterpreterPath_fd = "/Users/xudaxing/anaconda3/envs/ABBA_BABA/bin/python";
//        String parseVCFPy = "/Users/xudaxing/Software/ABBA_BABA/genomics_general/VCF_processing/parseVCF.py";
//        String ploidyFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt" +
//                "/021_Simulation/ploidy.txt";
//        String abbababaWindowPy = "/Users/xudaxing/Software/ABBA_BABA/genomics_general/ABBABABAwindows.py";
//        String popsFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/simulated_group.txt";
//        String rscriptInterpreterPath = "Rscript";
//        String transformFdResToIndividualLocalAncestryR = "/Users/xudaxing/IdeaProjects/PlantGenetics/GroupPlantGenetics/src/daxing/v2/localAncestryInfer/transformFdResToIndividualLocalAncestry.R";
//        Simulation.simulate_multiple_genotype_tract(pythonInterpreterPath, laidp_simulatePy, graphFileDir,
//                sequence_length, random_seed, mutation_rate, workingDirectory, logDir, threadsNum,
//                pythonInterpreterPath_fd, parseVCFPy, ploidyFile, abbababaWindowPy, popsFile, rscriptInterpreterPath,transformFdResToIndividualLocalAncestryR);

//        String pythonInterpreterPath = "/Users/xudaxing/anaconda3/envs/Msprime/bin/python";
//        String demeFilesDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/000_graphs";
//        String demesPy = "/Users/xudaxing/IdeaProjects/PlantGenetics/GroupPlantGenetics/src/daxing/v2/localAncestryInfer/Demes.py";
//        String graph_outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/000_graphs";
//        String workingDirectory = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation";
//        String logDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/log";
//        int threadsNum = 6;
//        Simulation.plot_deme(pythonInterpreterPath, demeFilesDir, demesPy, graph_outDir, workingDirectory, logDir, threadsNum);


//        String genotypeFileDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/001_simulatedGenotype";
//        String taxaInfoFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/taxaInfo.txt";
//        String recombinationMapDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/000_simulatedGenotype_recombinationMap";
//        String outDir_Mosaic = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/007_mosaic";
//        Mosaic_utils.transformVCF_to_Mosaic_inputFormat(genotypeFileDir,taxaInfoFile, recombinationMapDir, outDir_Mosaic);

        String fastFilePath = "/Users/xudaxing/Desktop/fastFiles/";
        String admixedPop = "French";
        String mosaicRunDir = "/Users/xudaxing/Desktop/mosaic/";
        String logFile = "/Users/xudaxing/Desktop/log.txt";
        int admixedIndividualNumber = 56;
        int variantsNum = 14863;
        double[][][] res= Mosaic_utils.build_mosaic_sh(fastFilePath, admixedPop, mosaicRunDir, logFile,
                admixedIndividualNumber,
                variantsNum, 2, 6, 18);
        System.out.println();

    }



}