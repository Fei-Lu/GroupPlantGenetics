package daxing;

import com.google.common.base.Joiner;
import com.google.common.primitives.Ints;
import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.demography.DemographicModelTools;
import daxing.v2.localAncestryInfer.evaluation.LocalAncestry;
import daxing.v2.localAncestryInfer.runner.*;
import daxing.v2.localAncestryInfer.simulation.Simulation;
import daxing.v2.localAncestryInfer.simulation.SimulationMetadata;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.BitSet;
import java.util.Date;
import java.util.List;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) throws IOException, InterruptedException {
//        String genotypeFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/003_simulation/D014.vcf";
//        int windowSize = 200;
//        int stepSize = 100;
//        String taxaGroupFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.taxaGroup.txt";
//        String ancestryAllele = "simulation";
//        int conjunctionNum = 2;
//        double switchCostScore = 1.5;
//        int maxSolutionCount = 32;
//        String localAnceOutFile = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014" +
//                ".localAnc_2.txt";
//        int threadsNum = 2;
////////
////////
////////        String outDir = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/005_twoWay/004_laidpRes/temp";
//
////        String genotypeFile = args[0];
////        int windowSize =Integer.parseInt(args[1]);
////        int stepSize = Integer.parseInt(args[2]);
////        String taxaGroupFile = args[3];
////        String ancestryAllele = args[4];
////        int conjunctionNum = Integer.parseInt(args[5]);
////        double switchCostScore = Double.parseDouble(args[6]);
////        String localAnceOutFile = args[7];
////        int threadsNum = Integer.parseInt(args[8]);
////////////////
//        long start = System.nanoTime();
//        GenotypeTable.run_LAIDP(genotypeFile, windowSize, stepSize, taxaGroupFile, ancestryAllele, conjunctionNum,
//                switchCostScore, maxSolutionCount, localAnceOutFile, threadsNum);
//
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s");
//
//        MD5.checkTwoFileMD5("/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.localAnc.txt","/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M/004_runner/laidp/D014/D014.localAnc_2.txt");
//        System.out.println(Benchmark.getTimeSpanSeconds(start)+" s ");

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_2";
//        String outDir2 = "/Users/xudaxing/Desktop/LAIDP_development/threeWay";
//        String outDir3 = "/Users/xudaxing/Desktop/fourWay";
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.THREE_WAY,outDir2);
//        evaluate_pearsonCorrelation_meanDeviation_contingencyTable(DemographicModelTools.N_way.FOUR_WAY, outDir3);

//        String outDir = "/Users/xudaxing/Desktop/LAIDP_development/twoWay_100M";
//        evaluate_contingencyTable(DemographicModelTools.N_way.TWO_WAY, outDir);

//        LAIDP_CLI.startFromCLI(args);
//        GenotypeTable genotypeTable = new GenotypeTable("/Users/xudaxing/Desktop/LAIDP_development/twoWay_recent_test/003_simulation/D001.vcf");
//        String title = args[0];
//        String shFile = args[1];
//        String workingDirectory = args[2];
//        String logDir = args[3];
//        int threadsNum = Integer.parseInt(args[4]);
//        CommandUtils.runSH_multipleCommands(title, shFile, workingDirectory, logDir, threadsNum);
//        String genoMatrixFile = args[0];
//        String outVCF = args[1];
//        transform2VCF(genoMatrixFile,outVCF);

//        rename(args[0]);
//        rename("/Users/xudaxing/Desktop/temp");

//        String title = args[0];
//        String shFile = args[1];
//        String workingDirectory = args[2];
//        String logDir = args[3];
//        boolean ifRedictStandoutToLog = args[4].equals("T") ? true : false;
//        int threadsNum = Integer.parseInt(args[5]);
//        CommandUtils.runSH_multipleCommands(title, shFile, workingDirectory, logDir, ifRedictStandoutToLog, threadsNum);

//        String a = "/Users/xudaxing/Desktop/btau_buffalo.txt";
//        String b = "/Users/xudaxing/Desktop/btau_buffalo2.txt";
//        temp(a, b);

//        String genomeFaFile = args[0];
////        String prefix = args[1];
//        String outFile = args[1];
//
//        FastaBit fastaBit= new FastaBit(genomeFaFile);
//        int seqNum = fastaBit.getSeqNumber();
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < seqNum; i++) {
//            sb.setLength(0);
////            sb.append(prefix).append("_");
//            sb.append(i+1);
//            fastaBit.setName(sb.toString(), i);
//        }
//        fastaBit.writeFasta(outFile, IOFileFormat.Text);
//        String title = args[0];
//        String shFile = args[1];
//        String workingDirectory = args[2];
//        String logDir = args[3];
////        boolean ifRedictStandoutToLog = Boolean.getBoolean(args[4]);
//        int threadsNum = Integer.parseInt(args[4]);
////        CommandUtils.runSH_multipleCommands(title, shFile, workingDirectory, logDir, ifRedictStandoutToLog, threadsNum);
//        CommandUtils.runSH_OneCommand(title, shFile, workingDirectory, logDir, threadsNum);

//        List<String> ab = WheatLineage.D.getChr();
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < ab.size(); i++) {
//            sb.append("chr").append(ab.get(i)).append(" ");
//        }
//        sb.deleteCharAt(sb.length()-1);
//        System.out.println(sb);
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

    public static void transform2VCF(String genoMatrixFile, String outVCF){
        try (BufferedReader br = IOTool.getReader(genoMatrixFile);
             BufferedWriter bw = IOTool.getWriter(outVCF)) {
            String line;
            List<String> temp, headerList;
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append("##fileformat=VCFv4.2\n");
            SimpleDateFormat f = new SimpleDateFormat("yyyyMMdd");
            sb.append("##fileDate=").append(f.format(new Date())).append("\n");
            sb.append("##FILTER=<ID=PASS,Description=\"All filters passed\">").append("\n");
            sb.append("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">").append("\n");
            sb.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t");
            headerList = PStringUtils.fastSplit(br.readLine());
            sb.append(String.join("\t", headerList.subList(8, headerList.size())));
            bw.write(sb.toString());
            bw.newLine();
            String refAllele, altAllele;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                refAllele = temp.get(2);
                sb.setLength(0);
                sb.append(temp.get(0)).append("\t");
                sb.append(temp.get(1)).append("\t");
                sb.append(temp.get(0)).append("_").append(temp.get(1)).append("\t");
                sb.append(temp.get(2)).append("\t");
                sb.append(temp.get(3)).append("\t");
                sb.append(".").append("\t");
                sb.append("PASS").append("\t");
                sb.append(".").append("\t");
//                sb.append("NS=").append(temp.get(4)).append(";");
//                sb.append("MF=").append(temp.get(5)).append("\t");
                sb.append("GT").append("\t");
                for (int i = 8; i < temp.size(); i++) {
                    if (temp.get(i).equals("-")){
                        sb.append(".|.").append("\t");
                    }else {
                        sb.append(temp.get(i).equals(refAllele) ? 0 : 1).append("|").append(temp.get(i).equals(refAllele) ? 0 : 1).append("\t");
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void rename(String faDir){
        List<File> files = IOTool.getFileListInDirEndsWith(faDir, ".fa.gz");
        int[] chrID = files.stream().map(File::getName).mapToInt(s -> Integer.parseInt(s.substring(3,6))).toArray();
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll("fa.gz","fa")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->{
            int chr = chrID[e];
            File file = files.get(e);
            try (BufferedReader br = IOTool.getReader(file);
                 BufferedWriter bw = IOTool.getWriter(new File(faDir, outNames[e]))) {
                boolean ifFirst = true;
                String line;
                StringBuilder sb = new StringBuilder();
                while ((line= br.readLine())!=null){
                    if (ifFirst){
                        sb.setLength(0);
                        sb.append(">chr").append(chr);
                        bw.write(sb.toString());
                        bw.newLine();
                        ifFirst= false;
                    }else {
                        bw.write(line);
                        bw.newLine();
                    }
                }
                bw.flush();
            }catch (Exception exception){
                exception.printStackTrace();
            }
        });
    }


    public static void temp(String inputFile, String outFile){
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            StringBuilder sb = new StringBuilder();
            int btauChr, buffaloChr;
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                btauChr = Integer.parseInt(temp.get(1));
                buffaloChr = Integer.parseInt(temp.get(0));
                sb.setLength(0);
                sb.append(PStringUtils.getNDigitNumber(3, btauChr)).append("\t");
                sb.append(PStringUtils.getNDigitNumber(3, buffaloChr));
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void splitFasta(String fastaFile, String outDir){
        try (BufferedReader br = IOTool.getReader(fastaFile)) {
            String line;
            List<String> temp;
            StringBuilder sb = new StringBuilder();
            while ((line=br.readLine())!=null){
                sb.setLength(0);
                if (line.startsWith(">")){
                    sb.append(">").append(line.substring(4));
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }



}