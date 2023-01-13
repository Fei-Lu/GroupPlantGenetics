package daxing.v2.localAncestryInfer.simulation;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Simulation {

    /**
     * Python packages need to be installed
     * demes, demesdraw, matplotlib, msprime, numpy
     */
    String pythonInterpreterPath;

    String demesPathFile;

    public static void runCommand(String command, String commandInfo){
        String[] commands = StringUtils.split(command, " ");
        ProcessBuilder processBuilder = new ProcessBuilder(commands);
        int exitCode=-1;
        try {
            Process process = processBuilder.start();
            exitCode = process.waitFor();
        } catch (IOException | InterruptedException e) {
            e.printStackTrace();
        }
        assert exitCode == 0: commandInfo+" run failed";
        System.out.println(commandInfo+ " finished");
    }

    public static String plot_deme(File pythonInterpreterPath, File demeFilesDir, File demesPy,
                                   File graph_outFile){
        StringBuilder sb = new StringBuilder();
        sb.append(pythonInterpreterPath.getAbsolutePath()).append(" ").append(demesPy.getAbsolutePath()).append(" ");
        sb.append(demeFilesDir.getAbsolutePath()).append(" ").append(graph_outFile.getAbsolutePath());
//        Simulation.runCommand(sb.toString(), "plot_deme");
        return sb.toString();
    }

    public static void plot_deme(String pythonInterpreterPath, String demeFilesDir, String demesPy,
                                 String graph_outDir, String workingDirectory,
                                 String logDir, int threadsNum){
        List<File> demeFiles = IOTool.getFileListInDirEndsWith(demeFilesDir, ".yaml");
        String[] outFileNames = demeFiles.stream().map(File::getName).map(s -> s.replaceAll("yaml","pdf")).toArray(String[]::new);
        List<String> commands = new ArrayList<>();
        for (int i = 0; i < demeFiles.size(); i++) {
            commands.add(Simulation.plot_deme(new File(pythonInterpreterPath), demeFiles.get(i), new File(demesPy),
                    new File(graph_outDir, outFileNames[i])));
        }
        CommandUtils.runSH_OneCommand("plot_demes", commands, workingDirectory, logDir, threadsNum);
    }

    public static String simulate_genotype_tract(String pythonInterpreterPath, String laidp_simulatePy,
                                                 String graphFile,
                                String outFile_simulateGenotype,
                                String outDir_tracts,
                                int sequence_length, int random_seed, double mutation_rate){
        StringBuilder sb = new StringBuilder();
        sb.append(pythonInterpreterPath).append(" ").append(laidp_simulatePy).append(" ");
        sb.append(graphFile).append(" ").append(outFile_simulateGenotype).append(" ");
        sb.append(outDir_tracts).append(" ").append(sequence_length).append(" ");
        sb.append(random_seed).append(" ").append(mutation_rate);
        return sb.toString();
//        Simulation.runCommand(sb.toString(), "simulate genotype and get tract");
    }

    public static void simulate_multiple_genotype_tract(String pythonInterpreterPath, String laidp_simulatePy,
                                                        String graphFileDir, int sequence_length, int random_seed,
                                                        double mutation_rate, String workingDirectory,
                                                        String logDir, int threadsNum,
                                                        String pythonInterpreterPath_fd, String parseVCFPy,
                                                        String ploidyFile, String abbababaWindowPy, String popsFile,
                                                        String rscriptInterpreterPath, String transformFdResToIndividualLocalAncestryR){
        List<File> graphs = IOTool.getFileListInDirEndsWith(graphFileDir, "yaml");
        double[] admixture_proportion = new double[graphs.size()];
        int[] admixture_generation = new int[graphs.size()];
        double[] divergenceTime = new double[graphs.size()];
        List<String> temp;
        for (int i = 0; i < graphs.size(); i++) {
            temp = PStringUtils.fastSplit(graphs.get(i).getName(), "_");
            admixture_proportion[i] = Double.parseDouble(temp.get(4).replaceAll("proportion",""));
            admixture_generation[i] = Integer.parseInt(temp.get(5).replaceAll("genetration",""));
            divergenceTime[i] = Double.parseDouble(temp.get(7).replaceAll("divergence", "").replaceAll("Ne.yaml",""));
        }
        List<String> commandsList = new ArrayList<>();
        String[] subDir = {"001_simulatedGenotype", "002_tracts", "003_fd","004_laidpRes"};
        String[] subDir_2 = new String[graphs.size()];
        File[] subDirFile  = new File[subDir.length];
        File[] subDir_2_File = new File[graphs.size()];
        for (int i = 0; i < graphs.size(); i++) {
            subDir_2[i] = graphs.get(i).getName().replaceAll(".yaml", "");
        }
        for (int i = 0; i < subDir.length; i++) {
            subDirFile[i] = new File(workingDirectory, subDir[i]);
            subDirFile[i].mkdir();
        }
        for (int i = 0; i < graphs.size(); i++) {
            subDir_2_File[i] = new File(subDirFile[1], subDir_2[i]);
            subDir_2_File[i].mkdir();
        }

//        for (int i = 0; i < graphs.size(); i++) {
//            commandsList.add(simulate_genotype_tract(pythonInterpreterPath, laidp_simulatePy, graphs.get(i).getAbsolutePath(),
//                    new File(subDirFile[0],"simulate_two_way_admixture_proportion"+admixture_proportion[i]+"_genetration"+admixture_generation[i]+
//                            "_P1P2_divergence"+divergenceTime[i]+"Ne.vcf").getAbsolutePath(),
//                    subDir_2_File[i].getAbsolutePath(), sequence_length,
//                    random_seed, mutation_rate));
//        }
//        CommandUtils.runSH("simulate_multiple", commandsList, workingDirectory, logDir, threadsNum);


//        List<File> genotypeFiles = IOTool.getFileListInDirEndsWith(subDirFile[0].getAbsolutePath(), "vcf");
//        String[] outFileName_GenotypeWithAncestral = genotypeFiles.stream().map(File::getName).map(s -> s.replaceAll(
//                ".vcf","_ancestral.vcf")).toArray(String[]::new);
//        for (int i = 0; i < graphs.size(); i++) {
//            Simulation.addAncestralSample(genotypeFiles.get(i).getAbsolutePath(), new File(subDirFile[0],
//                    outFileName_GenotypeWithAncestral[i]).getAbsolutePath());
//        }

        String[] subDir_fd = {"001_geno", "002_fdRes", "003_individualLocalAncestry"};
        File[] subDirFile_fd = new File[subDir_fd.length];
        for (int i = 0; i < subDir_fd.length; i++) {
            subDirFile_fd[i] = new File(subDirFile[2], subDir_fd[i]);
            subDirFile_fd[i].mkdir();
        }
//        List<File> genotypeWithAncestralFiles = IOTool.getFileListInDirEndsWith(subDirFile[0].getAbsolutePath(),
//                "_ancestral.vcf");
//        String[] outNames_GenotypeWithAncestralFiles =
//                genotypeWithAncestralFiles.stream().map(File::getName).map(s -> s.replaceAll("vcf", "geno")).toArray(String[]::new);
//        commandsList.clear();
//        for (int i = 0; i < genotypeWithAncestralFiles.size(); i++) {
//            commandsList.add(Simulation.transform_vcf2Geno(pythonInterpreterPath_fd, parseVCFPy,
//                    genotypeWithAncestralFiles.get(i).getAbsolutePath(), ploidyFile, new File(subDirFile_fd[0],
//                            outNames_GenotypeWithAncestralFiles[i]).getAbsolutePath()));
//        }
//        CommandUtils.runSH("transform_vcf2Geno", commandsList, workingDirectory, logDir, threadsNum);


//        List<File> genoFiles = IOTool.getFileListInDirEndsWith(subDirFile_fd[0].getAbsolutePath(), ".geno");
//        commandsList.clear();
//        int[] p2TaxaID= IntStream.range(60, 90).toArray();
//        int[] haploidID = IntStream.range(0, 90).toArray();
//        for (int i = 0; i < genoFiles.size(); i++) {
//            commandsList.addAll(Simulation.run_fd(pythonInterpreterPath_fd, abbababaWindowPy,
//                    genoFiles.get(i).getAbsolutePath(), subDirFile_fd[1].getAbsolutePath(), p2TaxaID, 4, 0,
//                    "ancestral", popsFile, haploidID, admixture_proportion[i], admixture_generation[i],
//                    divergenceTime[i]));
//        }
//        CommandUtils.runSH("fd", commandsList, workingDirectory, logDir, threadsNum);

//        Simulation.transformFdResToIndividualLocalAncestry(rscriptInterpreterPath,
//                transformFdResToIndividualLocalAncestryR, subDirFile_fd[1].getAbsolutePath()
//                , subDirFile_fd[2].getAbsolutePath());

    }

    public static void addAncestralSample(String simulatedGenotype, String simulatedGenotypeWithAncestral_outFile){
        try (BufferedReader br = IOTool.getReader(simulatedGenotype);
             BufferedWriter bw = IOTool.getWriter(simulatedGenotypeWithAncestral_outFile)) {
            String line;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            StringBuilder sb = new StringBuilder();
            sb.append(line).append("\t").append("ancestral");
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                sb.setLength(0);
                sb.append(line).append("\t").append(0);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static String transform_vcf2Geno(String pythonInterpreterPath, String parseVCFPy,
                                          String simulatedGenotypeWithAncestral,
                                          String ploidyFile,
                                          String simulatedGenotypeWithAncestralGeno_outFile){
        StringBuilder sb = new StringBuilder();
        sb.append(pythonInterpreterPath).append(" ");
        sb.append(parseVCFPy).append(" ");
        sb.append("-i").append(" ").append(simulatedGenotypeWithAncestral).append(" ");
        sb.append("--skipIndels --skipMono").append(" ").append("--ploidyFile ").append(ploidyFile).append(" ");
        sb.append("-o ").append(simulatedGenotypeWithAncestralGeno_outFile);
        return sb.toString();
    }

    public static List<String> run_fd(String pythonInterpreterPath, String abbababaWindowPy, String genoFile,
                              String outDir, int[] p2TaxaID, int p1PopID, int p3PopID, String outGroup,
                              String popsFile, int[] haploidID, double admixture_proportion, int admixture_generation
            , double divergenceTime){
        List<String> commandsList = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < p2TaxaID.length; i++) {
            sb.setLength(0);
            sb.append(pythonInterpreterPath).append(" ").append(abbababaWindowPy).append(" ");
            sb.append("--windType sites").append(" ").append("-g ").append(genoFile).append(" ");
            sb.append("-f haplo").append(" ");
            sb.append("-o ").append(outDir).append("/simulate_proportion").append(admixture_proportion).append("_");
            sb.append("genetration").append(admixture_generation).append("_divergence").append(divergenceTime).append("Ne");
            sb.append("_tsk_");
            sb.append(p2TaxaID[i]).append(".csv ");
            sb.append("-w 200 -m 10 --overlap 100 -P1 ").append(p1PopID).append(" ");
            sb.append("-P2 tsk_").append(p2TaxaID[i]).append(" ");
            sb.append("-P3 ").append(p3PopID).append(" ").append("-O ").append(outGroup).append(" ");
            sb.append("-T 1 ").append("--popsFile ").append(popsFile).append(" ");
            sb.append("--writeFailedWindows --haploid ");
            for (int hapID : haploidID){
                sb.append("tsk_").append(hapID).append(",");
            }
            sb.append(outGroup);
            commandsList.add(sb.toString());
        }
        return commandsList;
    }

    public static void transformFdResToIndividualLocalAncestry(String rscriptInterpreterPath,
                                                               String transformFdResToIndividualLocalAncestryR,
                                                               String fdResDir,
                                                               String individualLocalAncestryOutDir){
        StringBuilder sb = new StringBuilder();
        sb.append(rscriptInterpreterPath).append(" ");
        sb.append(transformFdResToIndividualLocalAncestryR).append(" ");
        sb.append(fdResDir).append(" ");
        sb.append(individualLocalAncestryOutDir);
        Simulation.runCommand(sb.toString(), "transformFdResToIndividualLocalAncestry");
    }


}
