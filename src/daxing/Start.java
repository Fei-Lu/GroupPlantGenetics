package daxing;

import daxing.v2.localAncestryInfer.Evaluation;
import daxing.v2.localAncestryInfer.LocalAncestryInferenceStart;
import daxing.v2.simulate.Utils;
import pgl.infra.utils.Benchmark;

import java.io.File;

public class Start {

    public static void main(String[] args) {

//        String genotypeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/002_simulatedGenotype/simulate_ancestral.vcf";
//        String fd_dxyFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/004_fd/004_individualLocalAncestry";
//        String groupByPop2IndividualFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/004_fd/001_group/simulated_group.txt";
//        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/005_laidpRes";
////////
////////////
//        int conjunctionNum=2; // 2
//        double switchCostScore=1.5; // 1.5
//        int maxSolutionCount=20; // 100
//        int threadNum= 10;
//////
//////        String genotypeFile = args[0];
//////        String fd_dxyFile = args[1];
//////        String groupByPop2IndividualFile=args[2];
//////        String outDir=args[3];
//////
//////        int conjunctionNum=Integer.parseInt(args[4]); // 2
//////        double switchCostScore=Double.parseDouble(args[5]); // 1.5
//////        int maxSolutionCount=Integer.parseInt(args[6]); // 100
////////        int threadNum= Integer.parseInt(args[7]);
//////
//////
//        long start = System.nanoTime();
//        LocalAncestryInferenceStart.inferLocalAncestry("1", new File(genotypeFile),
//                new File(groupByPop2IndividualFile), new File(fd_dxyFile), conjunctionNum, switchCostScore,
//                new File(outDir), maxSolutionCount);
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" seconds");
//
//
//        String simulatedTractDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/003_simulatedTract";
//        String genotypeTableFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/002_simulatedGenotype/simulate_proportion0.1_generation100_addAncestral.vcf";
//        String laidpResDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/005_laidpRes";
//        String outFile_accuracy = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/006_evaluation/contingencyTable.txt";
////        Evaluation.write_accuracy_recall_precision(simulatedTractDir, genotypeTableFile, laidpResDir, outFile_accuracy);
//        Evaluation.write_ContingencyTable(simulatedTractDir, genotypeTableFile, laidpResDir, outFile_accuracy);
//
//        Utils.addAncestral("/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/002_simulatedGenotype/simulate.vcf",
//                "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/002_simulatedGenotype/simulate_ancestral.vcf");
    }



}