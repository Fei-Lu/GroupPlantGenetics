package daxing.v2.localAncestryInfer.simulation;

import daxing.common.sh.CommandUtils;
import daxing.common.utiles.IOTool;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Callable;

public class Simulation {

    /**
     * Python packages need to be installed
     * demes, demesdraw, matplotlib, msprime, numpy
     */
    String pythonInterpreterPath;

    SimulationMetadata simulationMetadata;

    String logFile;
    String outDir;

    double mutation_rate;
    int random_seed;
    int threadsNum;

    public Simulation(String pythonInterpreterPath, String simulationMetadata, String logFile,
                      String outDir, double mutation_rate, int random_seed, int threadsNum){

        this.pythonInterpreterPath=pythonInterpreterPath;
        this.simulationMetadata=new SimulationMetadata(simulationMetadata);
        this.logFile=logFile;
        this.outDir=outDir;
        this.mutation_rate=mutation_rate;
        this.random_seed=random_seed;
        this.threadsNum=threadsNum;
    }

    public Simulation(Builder builder){
        this.pythonInterpreterPath=builder.pythonInterpreterPath;
        this.simulationMetadata=builder.simulationMetadata;
        this.logFile=builder.logFile;
        this.outDir=builder.outDir;
        this.mutation_rate=builder.mutation_rate;
        this.random_seed=builder.random_seed;
        this.threadsNum=builder.threadsNum;
    }

    public Simulation(String parameterFile){
        this.initialize(parameterFile);
        this.run_simulation();
    }

    private void initialize(String parameterFile){
        try (BufferedReader br = IOTool.getReader(parameterFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                if (line.startsWith("##")) continue;
                temp = PStringUtils.fastSplit(line, ":");
                if (line.startsWith("pythonInterpreterPath")){
                    this.pythonInterpreterPath=temp.get(1);
                    continue;
                }
                if (line.startsWith("SimulationMetadata")){
                    this.simulationMetadata = new SimulationMetadata(temp.get(1));
                    continue;
                }
                if (line.startsWith("LogFilePath")){
                    this.logFile=temp.get(1);
                    continue;
                }
                if (line.startsWith("OutDir")){
                    this.outDir=temp.get(1);
                    continue;
                }
                if (line.startsWith("mutation_rate")){
                    this.mutation_rate=Double.parseDouble(temp.get(1));
                    continue;
                }
                if (line.startsWith("random_seed")){
                    this.random_seed=Integer.parseInt(temp.get(1));
                    continue;
                }
                if (line.startsWith("threadsNum")){
                    this.threadsNum=Integer.parseInt(temp.get(1));
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static class Builder {

        /**
         * Required parameters
         */
        private SimulationMetadata simulationMetadata;
        private String logFile;
        private String outDir;

        /**
         * Optional parameters
         */
        private String pythonInterpreterPath = "/Users/xudaxing/anaconda3/envs/Msprime/bin/python";
        private double mutation_rate = 7.1e-9;
        private int random_seed = 1;
        private int threadsNum = 2;

        public Builder(String simulationMetadata, String logFile, String outDir){
            this.simulationMetadata=new SimulationMetadata(simulationMetadata);
            this.logFile=logFile;
            this.outDir=outDir;
        }

        public Builder pythonInterpreterPath(String pythonInterpreterPath){
            this.pythonInterpreterPath=pythonInterpreterPath;
            return this;
        }

        public Builder mutation_rate(double mutation_rate){
            this.mutation_rate=mutation_rate;
            return this;
        }

        public Builder random_seed(int random_seed){
            this.random_seed=random_seed;
            return this;
        }

        public Builder threadsNum(int threadsNum){
            this.threadsNum=threadsNum;
            return this;
        }

        public Simulation build(){
            return new Simulation(this);
        }

    }

    public int getTotalNumberOfSimulation(){
        return this.simulationMetadata.demesID.length;
    }

    public String getDemesID(int simulationIndex){
        return this.simulationMetadata.demesID[simulationIndex];
    }

    public String getDemesPath(int simulationIndex){
        return this.simulationMetadata.demesPath[simulationIndex];
    }

    public String getAdmixedPop(int simulationIndex){
        return this.simulationMetadata.admixedPop[simulationIndex];
    }

    public String getNativePop(int simulationIndex){
        return this.simulationMetadata.nativePop[simulationIndex];
    }

    public List<String> getIntrogressedPop(int simulationIndex){
        return this.simulationMetadata.introgressedPop[simulationIndex];
    }

    public int getSampleSizeOfAdmixedPop(int simulationIndex){
        return this.simulationMetadata.admixedPopSampleSize[simulationIndex];
    }

    public int getSampleSizeOfNativePop(int simulationIndex){
        return this.simulationMetadata.nativePopSampleSize[simulationIndex];
    }

    public IntList getSampleSizeOfIntrogressedPop(int simulationIndex){
        return this.simulationMetadata.introgressedPopSampleSize[simulationIndex];
    }

    public int getSeqLen(int simulationIndex){
        return this.simulationMetadata.sequenceLen[simulationIndex];
    }

    public double getRecombinationRate(int simulationIndex){
        return this.simulationMetadata.recombinationRate[simulationIndex];
    }

    public int getRandom_seed() {
        return random_seed;
    }

    public double getMutation_rate() {
        return mutation_rate;
    }

    public String getOutDir() {
        return outDir;
    }

    public String getLogFile() {
        return logFile;
    }

    public int getThreadsNum() {
        return threadsNum;
    }

    public void run_simulation(){
        new File(this.logFile).delete();
        long start = System.nanoTime();
        String simulationPyPath = Simulation.class.getResource("Simulation.py").getPath();
        List<Callable<Integer>> callableList = new ArrayList<>();
        for (int i = 0; i < this.getTotalNumberOfSimulation(); i++) {
            StringBuilder sb = new StringBuilder();
            sb.setLength(0);
            sb.append(this.pythonInterpreterPath).append(" ");
            sb.append(simulationPyPath).append(" ");
            sb.append(this.getDemesPath(i)).append(" ");
            sb.append(new File(outDir, this.getDemesID(i)+".graph.pdf").getAbsolutePath()).append(" ");
            sb.append(new File(outDir, this.getDemesID(i)+".vcf")).append(" ");
            sb.append(new File(outDir, this.getDemesID(i)+".taxaInfo")).append(" ");
            sb.append(new File(outDir, this.getDemesID(i)+".recombinationMap")).append(" ");
            sb.append(new File(outDir, this.getDemesID(i)+".tract")).append(" ");
            sb.append(this.getSeqLen(i)).append(" ");
            sb.append(this.getRandom_seed()).append(" ");
            sb.append(this.getMutation_rate()).append(" ");
            sb.append(this.getRecombinationRate(i)).append(" ");
            sb.append(this.getAdmixedPop(i)).append(",");
            sb.append(this.getNativePop(i)).append(",");
            sb.append(String.join(",", this.getIntrogressedPop(i))).append(" ");
            sb.append(this.getSampleSizeOfAdmixedPop(i)).append(",");
            sb.append(this.getSampleSizeOfNativePop(i)).append(",");
            sb.append(StringUtils.join(this.getSampleSizeOfIntrogressedPop(i), ","));
            callableList.add(()->CommandUtils.runOneCommand(sb.toString(), this.getOutDir(),
                    new File(this.getLogFile())));
        }
        List<Integer> results = CommandUtils.run_commands(callableList, this.getThreadsNum());
        int failedCommandCount=0;
        for (int i = 0; i < results.size(); i++) {
            if (results.get(i)!=0){
                System.out.println(this.getDemesID(i)+" simulation run failed");
                failedCommandCount++;
            }
        }
        if (failedCommandCount > 0){
            System.out.println("Total "+failedCommandCount+" simulations run failed");
        }else {
            System.out.println("All simulation run successful");
        }
        System.out.println("simulation runner completed in "+ Benchmark.getTimeSpanSeconds(start)+" seconds");
    }

//    public static void addAncestralSample(String simulatedGenotype, String simulatedGenotypeWithAncestral_outFile){
//        try (BufferedReader br = IOTool.getReader(simulatedGenotype);
//             BufferedWriter bw = IOTool.getWriter(simulatedGenotypeWithAncestral_outFile)) {
//            String line;
//            while ((line=br.readLine()).startsWith("##")){
//                bw.write(line);
//                bw.newLine();
//            }
//            StringBuilder sb = new StringBuilder();
//            sb.append(line).append("\t").append("ancestral");
//            bw.write(sb.toString());
//            bw.newLine();
//            while ((line=br.readLine())!=null){
//                sb.setLength(0);
//                sb.append(line).append("\t").append(0);
//                bw.write(sb.toString());
//                bw.newLine();
//            }
//            bw.flush();
//        } catch (IOException e) {
//            e.printStackTrace();
//        }
//    }
//
//    public static String transform_vcf2Geno(String pythonInterpreterPath, String parseVCFPy,
//                                          String simulatedGenotypeWithAncestral,
//                                          String ploidyFile,
//                                          String simulatedGenotypeWithAncestralGeno_outFile){
//        StringBuilder sb = new StringBuilder();
//        sb.append(pythonInterpreterPath).append(" ");
//        sb.append(parseVCFPy).append(" ");
//        sb.append("-i").append(" ").append(simulatedGenotypeWithAncestral).append(" ");
//        sb.append("--skipIndels --skipMono").append(" ").append("--ploidyFile ").append(ploidyFile).append(" ");
//        sb.append("-o ").append(simulatedGenotypeWithAncestralGeno_outFile);
//        return sb.toString();
//    }
//
//    public static List<String> run_fd(String pythonInterpreterPath, String abbababaWindowPy, String genoFile,
//                              String outDir, int[] p2TaxaID, int p1PopID, int p3PopID, String outGroup,
//                              String popsFile, int[] haploidID, double admixture_proportion, int admixture_generation
//            , double divergenceTime){
//        List<String> commandsList = new ArrayList<>();
//        StringBuilder sb = new StringBuilder();
//        for (int i = 0; i < p2TaxaID.length; i++) {
//            sb.setLength(0);
//            sb.append(pythonInterpreterPath).append(" ").append(abbababaWindowPy).append(" ");
//            sb.append("--windType sites").append(" ").append("-g ").append(genoFile).append(" ");
//            sb.append("-f haplo").append(" ");
//            sb.append("-o ").append(outDir).append("/simulate_proportion").append(admixture_proportion).append("_");
//            sb.append("genetration").append(admixture_generation).append("_divergence").append(divergenceTime).append("Ne");
//            sb.append("_tsk_");
//            sb.append(p2TaxaID[i]).append(".csv ");
//            sb.append("-w 200 -m 10 --overlap 100 -P1 ").append(p1PopID).append(" ");
//            sb.append("-P2 tsk_").append(p2TaxaID[i]).append(" ");
//            sb.append("-P3 ").append(p3PopID).append(" ").append("-O ").append(outGroup).append(" ");
//            sb.append("-T 1 ").append("--popsFile ").append(popsFile).append(" ");
//            sb.append("--writeFailedWindows --haploid ");
//            for (int hapID : haploidID){
//                sb.append("tsk_").append(hapID).append(",");
//            }
//            sb.append(outGroup);
//            commandsList.add(sb.toString());
//        }
//        return commandsList;
//    }
//
//    public static void transformFdResToIndividualLocalAncestry(String rscriptInterpreterPath,
//                                                               String transformFdResToIndividualLocalAncestryR,
//                                                               String fdResDir,
//                                                               String individualLocalAncestryOutDir){
//        StringBuilder sb = new StringBuilder();
//        sb.append(rscriptInterpreterPath).append(" ");
//        sb.append(transformFdResToIndividualLocalAncestryR).append(" ");
//        sb.append(fdResDir).append(" ");
//        sb.append(individualLocalAncestryOutDir);
//        Simulation.runCommand(sb.toString(), "transformFdResToIndividualLocalAncestry");
//    }


}
