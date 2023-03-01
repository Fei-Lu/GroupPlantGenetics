package daxing.v2.localAncestryInfer.demography;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;
import com.google.common.primitives.Ints;
import daxing.common.utiles.IOTool;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class DemographicModelTools {

    public enum N_way{
        TWO_WAY,
        THREE_WAY,
        FOUR_WAY
    }

    public static DemographicModel readModel(File modelPath){
        DemographicModel demographicModel;
        ObjectMapper mapper = new ObjectMapper(new YAMLFactory());
        mapper.findAndRegisterModules();
        try {
            demographicModel = mapper.readValue(modelPath, DemographicModel.class);
            demographicModel.trim_default();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return demographicModel;
    }

    public static void writeModel(DemographicModel demographicModel, File outFile_model){
        YAMLFactory yamlFactory = new YAMLFactory();
        yamlFactory.disable(YAMLGenerator.Feature.WRITE_DOC_START_MARKER);
        yamlFactory.enable(YAMLGenerator.Feature.MINIMIZE_QUOTES);
        yamlFactory.enable(YAMLGenerator.Feature.INDENT_ARRAYS_WITH_INDICATOR);
        ObjectMapper mapper = new ObjectMapper(yamlFactory);
        mapper.findAndRegisterModules();
        mapper.disable(SerializationFeature.WRITE_DATES_AS_TIMESTAMPS);
        mapper.setSerializationInclusion(JsonInclude.Include.NON_DEFAULT);
        try {
            mapper.writeValue(outFile_model, demographicModel);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    /**
     *
     * When multiple pulses are specified with the same time, the migration pulses occur in the order in which they are written.
     * Consider the following two pulses into deme A at time 100.
     * pulses:
     * - sources: [B]
     *   dest: A
     *   time: 100
     *   proportions: [0.1]
     * - sources: [C]
     *   dest: A
     *   time: 100
     *   proportions: [0.2]
     * The second pulse replaces 20% of A’s ancestry, including 20% of the ancestry that was inherited from B in the first pulse.
     * So immediately after time 100, A has 20% ancestry from C but only 8% ancestry from B.
     * @param admixtureProportion
     * @return transformedAdmixtureProportion
     */
    public static double[] transform_proportion(double[] admixtureProportion){
        assert  admixtureProportion.length > 0: "error, check admixtureProportion";

        if (admixtureProportion.length  > 1 ){
            double[] transformedAdmixtureProportion = new double[admixtureProportion.length];
            for (int i = 0; i < admixtureProportion.length; i++){
                double proportion=1;
                for (int j = i+1; j < admixtureProportion.length; j++) {
                    proportion*=1-admixtureProportion[j];
                }
                transformedAdmixtureProportion[i] = admixtureProportion[i]/proportion;
            }
            return transformedAdmixtureProportion;
        }else return admixtureProportion;
    }

    private static List<Deme> getEquilibriumPopulationDemes(String[] demeNames, int[] epochEndTime,
                                                            String[] ancestors, int populationSize){
        List<Deme> demes = new ArrayList<>();
        List<Epoch> epoches;
        Epoch epoch;
        Deme deme;
        List<String> ancestorsList;
        for (int i = 0; i < demeNames.length; i++) {
            epoches = new ArrayList<>();
            epoch = new Epoch(populationSize, epochEndTime[i], populationSize);
            epoches.add(epoch);
            if (i > 0){
                ancestorsList = new ArrayList<>();
                ancestorsList.add(ancestors[i]);
                deme = new Deme.Builder(demeNames[i], epoches).ancestors(ancestorsList).build();
            }else {
                deme = new Deme.Builder(demeNames[i], epoches).build();
            }
            demes.add(deme);
        }
        return demes;
    }

//    private static List<Deme> getTwoWayNonEquilibriumPopulationDemes(int[] splitEventTime, int populationSize){
//        String[] demeNames = {"A","B","C","D","E"};
//        String[] ancestors = {null,"A","A","B","B"};
//        int[] epochEndTime = new int[splitEventTime.length * 2 +1];
//        System.arraycopy(splitEventTime, 0, epochEndTime, 0, splitEventTime.length);
//
//        List<Deme> demes = new ArrayList<>();
//        List<Epoch> epoches;
//        Epoch epoch;
//        Deme deme;
//        List<String> ancestorsList;
//        for (int i = 0; i < demeNames.length; i++) {
//            epoches = new ArrayList<>();
//            epoch = new Epoch(populationSize, epochEndTime[i], populationSize);
//            epoches.add(epoch);
//            if (i > 0){
//                ancestorsList = new ArrayList<>();
//                ancestorsList.add(ancestors[i]);
//                deme = new Deme.Builder(demeNames[i], epoches).ancestors(ancestorsList).build();
//                if (demeNames[i].equals("C")){
//                    epoch = new Epoch();
//                }
//            }else {
//                deme = new Deme.Builder(demeNames[i], epoches).build();
//            }
//            demes.add(deme);
//        }
//
//        return demes;
//    }

    private static List<Pulse> getPulses(String[] sourcePop, String destPop,
                                         double[] admixtureProportion, int[] admixtureTime){
        List<String> sourcePopList = new ArrayList<>();
        sourcePopList.addAll(Arrays.asList(sourcePop));
        List<Double> admixtureProportionList = new ArrayList<>();
        for (double proportion : DemographicModelTools.transform_proportion(admixtureProportion)){
            admixtureProportionList.add(proportion);
        }
        List<Pulse> pulses = new ArrayList<>();
        Pulse pulse;
        for (int i = 0; i < sourcePopList.size(); i++) {
            pulse = new Pulse(sourcePopList.subList(i,i+1), destPop, admixtureProportionList.subList(i,i+1), admixtureTime[i]);
            pulses.add(pulse);
        }
        return pulses;
    }


    /**
     * support two-way and even multiple-way admixture
     * @param demeNames
     * @param splitEventTime split event time, from past to present, such ad 10000, 9000, 1000
     * @param ancestors only one ancestor per population
     * @param populationSize fixed
     * @param sourcePop source population per admixed population
     * @param destPop desPop
     * @param admixtureProportion admixtureProportion
     * @param admixtureTime admixtureTime
     * @return DemographicModel
     */
    public static DemographicModel equilibriumPopulationModelBuilder(String[] demeNames, int[] splitEventTime,
                                                                      String[] ancestors, int populationSize,
                                                                      String[] sourcePop, String destPop,
                                                                      double[] admixtureProportion,
                                                                      int[] admixtureTime){
        int[] epochEndTime = new int[splitEventTime.length * 2 +1];
        System.arraycopy(splitEventTime, 0, epochEndTime, 0, splitEventTime.length);
        List<Deme> demes = DemographicModelTools.getEquilibriumPopulationDemes(demeNames, epochEndTime, ancestors, populationSize);
        List<Pulse> pulses = DemographicModelTools.getPulses(sourcePop, destPop,
                admixtureProportion, admixtureTime);
        return new DemographicModel.Builder("generations", demes).pulses(pulses).build();
    }

    /**
     *
     * @param splitEventTime 2 dimensional array for two-way admixture
     * @param admixtureProportion
     * @param admixtureTime
     * @return DemographicModel
     */
    public static DemographicModel twoWayBuilder_equilibriumPopulation(int[] splitEventTime,
                                                                     double admixtureProportion,
                                                                     int admixtureTime){
        String[] demeNames = {"A","B","C","D","E"};
        String[] ancestors = {null,"A","A","B","B"};
        int populationSize = 10000;
        String[] sourcePop = new String[1];
        sourcePop[0] = "C";
        String destPop = "E";
        double[] admixtureProportionArray = new double[1];
        admixtureProportionArray[0] = admixtureProportion;
        int[] admixtureTimeArray = new int[1];
        admixtureTimeArray[0] = admixtureTime;

        return DemographicModelTools.equilibriumPopulationModelBuilder(demeNames, splitEventTime, ancestors,
                populationSize, sourcePop, destPop, admixtureProportionArray, admixtureTimeArray);
    }

    /**
     *
     * @param splitEventTime 3 dimensional array for three-way admixture
     * @param admixtureProportion 2 dimensional array for three-way admixture
     * @param admixtureTime 2 dimensional array for three-way admixture
     * @return DemographicModel
     */
    public static DemographicModel threeWayBuilder_equilibriumPopulation(int[] splitEventTime,
                                                                       double[] admixtureProportion,
                                                                       int[] admixtureTime){
        String[] demeNames = {"A","B","C","D","E","F","G"};
        String[] ancestors = {null,"A","B","A","B","C","C"};
        int populationSize = 10000;
        String[] sourcePop = new String[2];
        sourcePop[0] = "D";
        sourcePop[1] = "E";
        String destPop = "F";
        return DemographicModelTools.equilibriumPopulationModelBuilder(demeNames, splitEventTime, ancestors,
                populationSize, sourcePop, destPop, admixtureProportion, admixtureTime);
    }

    public static DemographicModel fourWayBuilder_equilibriumPopulation(int[] splitEventTime,
                                                                         double[] admixtureProportion,
                                                                         int[] admixtureTime){
        String[] demeNames = {"A","B","C","D","E","F","G","H","I"};
        String[] ancestors = {null,"A","B","C","A","B","C","D","D"};
        int populationSize = 10000;
        String[] sourcePop = new String[3];
        sourcePop[0] = "E";
        sourcePop[1] = "F";
        sourcePop[2] = "G";
        String destPop = "H";
        return DemographicModelTools.equilibriumPopulationModelBuilder(demeNames, splitEventTime, ancestors,
                populationSize, sourcePop, destPop, admixtureProportion, admixtureTime);
    }

    public static void batchRun(N_way nWay, String simulationMetadataOutFile, String outDir){
        switch (nWay){
            case TWO_WAY:
                batchRun_twoWay(simulationMetadataOutFile, outDir);
                break;
            case THREE_WAY:
                batchRun_threeWay(simulationMetadataOutFile, outDir);
                break;
            case FOUR_WAY:
                batchRun_fourWay(simulationMetadataOutFile, outDir);
                break;
        }
    }

    private static void batchRun_twoWay(String simulationMetadataOutFile, String outDir){
//        String[] admixed_native_introgressed_pop = {"E","D","C"};
//        int[] sampleSize = {30,30,30};
//        int equilibriumPopulationSize = 10000;
//        int[] splitEventTime_0 = {10000};
//        int[] splitEventTime_1 = {500, 1000, 2000, 4000, 8000};
//        double[] ratio_admixture_divergence = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64};
//        double[] admixtureProportion = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32};
//        int sequenceLength = 100_000_000;
//        double recombination_rate = 1e-8;
//        double mutation_rate = 1e-8; // 7.1e-9

        String[] admixed_native_introgressed_pop = {"E","D","C"};
        int[] sampleSize = {30,30,30};
        int equilibriumPopulationSize = 10000;
        int[] splitEventTime_0 = {10000};
        int[] splitEventTime_1 = {8000};
        double[] ratio_admixture_divergence = {0.1, 0.2, 0.4, 0.8};
        double[] admixtureProportion = {0.1};
        int sequenceLength = 1_000_000;
        double recombination_rate = 1e-8;
        double mutation_rate = 1e-8; // 7.1e-9


        int[] splitEventTime;
        List<DemographicModel> demographicModels = new ArrayList<>();
        List<String> modelName = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        DemographicModel demographicModel;
        double admixtureTime, divergenceTime0, divergenceTime1;
        StringBuilder sb_Metadata = new StringBuilder();
        String header = "DemesID\tDemePath\tAdmixedPop\tNativePop\tIntrogressedPopulation\tSequenceLength" +
                "\tRecombination_rate\tMutation_rate";
        try (BufferedWriter bw = IOTool.getWriter(simulationMetadataOutFile)) {
            bw.write(header);
            bw.newLine();
            int count= 1;
            for (int j : splitEventTime_0) {
                for (int k : splitEventTime_1) {
                    for (double value : ratio_admixture_divergence) {
                        for (double v : admixtureProportion) {
                            splitEventTime = new int[2];
                            splitEventTime[0] = j;
                            splitEventTime[1] = k;
                            admixtureTime = k * value;
                            divergenceTime0 = j;
                            divergenceTime1 = k;
                            demographicModel = DemographicModelTools.twoWayBuilder_equilibriumPopulation(splitEventTime,
                                    v, (int) admixtureTime);
                            demographicModel.trim_default();
                            demographicModels.add(demographicModel);

                            // model name
                            sb.setLength(0);
                            sb.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("_");
                            sb.append("way_").append(2).append("_").append("adp_").append(admixed_native_introgressed_pop[0]).append("_");
                            sb.append("nap_").append(admixed_native_introgressed_pop[1]).append("_");
                            sb.append("inp_").append(admixed_native_introgressed_pop[2]).append("_");
                            sb.append("ep_").append(equilibriumPopulationSize).append("_");
                            sb.append("pr_").append(v).append("_");
                            sb.append("at_").append((int) admixtureTime).append("_");
                            sb.append("dta_").append((int) divergenceTime0).append("_");
                            sb.append("dtb_").append((int) divergenceTime1).append(".yaml");
                            modelName.add(sb.toString());

                            // simulation metadata
                            sb_Metadata.setLength(0);
                            sb_Metadata.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("\t");
                            sb_Metadata.append(new File(outDir, sb.toString()).getAbsolutePath()).append("\t");
                            sb_Metadata.append(admixed_native_introgressed_pop[0]).append(":").append(sampleSize[0]).append("\t");
                            sb_Metadata.append(admixed_native_introgressed_pop[1]).append(":").append(sampleSize[1]).append("\t");
                            sb_Metadata.append(admixed_native_introgressed_pop[2]).append(":").append(sampleSize[2]).append("\t");
                            sb_Metadata.append(sequenceLength).append("\t").append(recombination_rate).append("\t");
                            sb_Metadata.append(mutation_rate);
                            bw.write(sb_Metadata.toString());
                            bw.newLine();
                            count++;
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        for (int i = 0; i < demographicModels.size(); i++) {
            DemographicModelTools.writeModel(demographicModels.get(i), new File(outDir, modelName.get(i)));
        }


    }

    private static void batchRun_threeWay(String simulationMetadataOutFile, String outDir){

        String[] admixed_native_introgressed_pop = {"F","G","D","E"};
        int[] sampleSize = {30,30,30,30};
        int equilibriumPopulationSize = 10000;
        int[] splitEventTime_0 = {10000}; // D
        int[] splitEventTime_1 = {8000}; // E
        int[] splitEventTime_2 = {500};  // F and G
        double[] ratio_admixture_divergence_0 = {0.4, 0.3}; // D
        double[] ratio_admixture_divergence_1 = {0.2, 0.1}; // E
        double[] admixtureProportion_0 = {0.1, 0.2}; // D
        double[] admixtureProportion_1 = {0.2, 0.3}; // E
        int sequenceLength = 1_000_000;
        double recombination_rate = 1e-8;
        double mutation_rate = 1e-8; // 7.1e-9

//        String[] admixed_native_introgressed_pop = {"F","G","D","E"};
//        int[] sampleSize = {30,30,30,30};
//        int equilibriumPopulationSize = 10000;
//        int[] splitEventTime_0 = {10000};
//        int[] splitEventTime_1 = {8000};
//        int[] splitEventTime_2 = {500};
//        double[] ratio_admixture_divergence_0 = {0.3};
//        double[] ratio_admixture_divergence_1 = {0.1};
//        double[] admixtureProportion_0 = {0.1};
//        double[] admixtureProportion_1 = {0.2};
//        int sequenceLength = 1_000_000;
//        double recombination_rate = 1e-8;
//        double mutation_rate = 1e-8; // 7.1e-9


        int[] splitEventTime;
        List<DemographicModel> demographicModels = new ArrayList<>();
        List<String> modelName = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        DemographicModel demographicModel;
        int[] admixtureTime;
        double[] admixtureProportion;
        int divergenceTime0, divergenceTime1, divergenceTime2;
        StringBuilder sb_Metadata = new StringBuilder();
        String header = "DemesID\tDemePath\tAdmixedPop\tNativePop\tIntrogressedPopulation\tSequenceLength" +
                "\tRecombination_rate\tMutation_rate";
        try (BufferedWriter bw = IOTool.getWriter(simulationMetadataOutFile)) {
            bw.write(header);
            bw.newLine();
            int count= 1;
            for (int value : splitEventTime_0) {
                for (int i : splitEventTime_1) {
                    for (int j : splitEventTime_2) {
                        for (int k = 0; k < ratio_admixture_divergence_0.length; k++) {
                            for (int l = 0; l < admixtureProportion_0.length; l++) {
                                splitEventTime = new int[3];
                                splitEventTime[0] = value;
                                splitEventTime[1] = i;
                                splitEventTime[2] = j;
                                admixtureTime = new int[2];
                                admixtureTime[0] = (int) (j * ratio_admixture_divergence_0[k]);
                                admixtureTime[1] = (int) (j * ratio_admixture_divergence_1[k]);
                                admixtureProportion = new double[2];
                                admixtureProportion[0] = admixtureProportion_0[l];
                                admixtureProportion[1] = admixtureProportion_1[l];
                                divergenceTime0 = value;
                                divergenceTime1 = i;
                                divergenceTime2 = j;
                                demographicModel =
                                        DemographicModelTools.threeWayBuilder_equilibriumPopulation(splitEventTime,
                                                admixtureProportion, admixtureTime);
                                demographicModel.trim_default();
                                demographicModels.add(demographicModel);

                                // model name
                                sb.setLength(0);
                                sb.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("_");
                                sb.append("way_").append(3).append("_").append("adp_").append(admixed_native_introgressed_pop[0]).append("_");
                                sb.append("nap_").append(admixed_native_introgressed_pop[1]).append("_");
                                sb.append("in0_").append(admixed_native_introgressed_pop[2]).append("_");
                                sb.append("in1_").append(admixed_native_introgressed_pop[3]).append("_");
                                sb.append("ep_").append(equilibriumPopulationSize).append("_");
                                sb.append("pr0_").append(admixtureProportion_0[l]).append("_");
                                sb.append("pr1_").append(admixtureProportion_1[l]).append("_");
                                sb.append("at0_").append(admixtureTime[0]).append("_");
                                sb.append("at1_").append(admixtureTime[1]).append("_");
                                sb.append("dt0_").append(divergenceTime0).append("_");
                                sb.append("dt1_").append(divergenceTime1).append("_");
                                sb.append("dt2_").append(divergenceTime2);
                                sb.append(".yaml");
                                modelName.add(sb.toString());

                                // simulation metadata
                                String[] introgressedPop = new String[2];
                                int[] introgressedPop_sampleSize = new int[2];
                                System.arraycopy(admixed_native_introgressed_pop, 2, introgressedPop, 0, 2);
                                System.arraycopy(sampleSize, 2, introgressedPop_sampleSize, 0, 2);
                                sb_Metadata.setLength(0);
                                sb_Metadata.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("\t");
                                sb_Metadata.append(new File(outDir, sb.toString()).getAbsolutePath()).append("\t");
                                sb_Metadata.append(admixed_native_introgressed_pop[0]).append(":").append(sampleSize[0]).append("\t");
                                sb_Metadata.append(admixed_native_introgressed_pop[1]).append(":").append(sampleSize[1]).append("\t");
                                sb_Metadata.append(String.join(",", introgressedPop)).append(":");
                                sb_Metadata.append(Ints.join(",", introgressedPop_sampleSize)).append("\t");
                                sb_Metadata.append(sequenceLength).append("\t").append(recombination_rate).append("\t");
                                sb_Metadata.append(mutation_rate);
                                bw.write(sb_Metadata.toString());
                                bw.newLine();
                                count++;
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        for (int i = 0; i < demographicModels.size(); i++) {
            DemographicModelTools.writeModel(demographicModels.get(i), new File(outDir, modelName.get(i)));
        }
    }

    private static void batchRun_fourWay(String simulationMetadataOutFile, String outDir){

        String[] admixed_native_introgressed_pop = {"H","I","E","F","G"};
        int[] sampleSize = {30,30,30,30,30};
        int equilibriumPopulationSize = 10000;
        int[] splitEventTime_0 = {10000};
        int[] splitEventTime_1 = {9000};
        int[] splitEventTime_2 = {8000};
        int[] splitEventTime_3 = {100, 1000, 4000};
        double[] ratio_admixture_divergence_0 = {0.3};
        double[] ratio_admixture_divergence_1 = {0.2};
        double[] ratio_admixture_divergence_2 = {0.1};
        double[] admixtureProportion_0 = {0.1};
        double[] admixtureProportion_1 = {0.2};
        double[] admixtureProportion_2 = {0.3};
        int sequenceLength = 1_000_000;
        double recombination_rate = 1e-8;
        double mutation_rate = 1e-8; // 7.1e-9

        int n_way = 4;

        int[] splitEventTime;
        List<DemographicModel> demographicModels = new ArrayList<>();
        List<String> modelName = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        DemographicModel demographicModel;
        int[] admixtureTime;
        double[] admixtureProportion;
        int divergenceTime0, divergenceTime1, divergenceTime2, divergenceTime3;
        StringBuilder sb_Metadata = new StringBuilder();
        String header = "DemesID\tDemePath\tAdmixedPop\tNativePop\tIntrogressedPopulation\tSequenceLength" +
                "\tRecombination_rate\tMutation_rate";
        try (BufferedWriter bw = IOTool.getWriter(simulationMetadataOutFile)) {
            bw.write(header);
            bw.newLine();
            int count= 1;
            for (int j : splitEventTime_0) {
                for (int element : splitEventTime_1) {
                    for (int item : splitEventTime_2) {
                        for (int value : splitEventTime_3) {
                            for (int k = 0; k < ratio_admixture_divergence_0.length; k++) {
                                for (int l = 0; l < admixtureProportion_0.length; l++) {
                                    splitEventTime = new int[n_way];
                                    splitEventTime[0] = j;
                                    splitEventTime[1] = element;
                                    splitEventTime[2] = item;
                                    splitEventTime[3] = value;
                                    admixtureTime = new int[n_way - 1];
                                    admixtureTime[0] = (int) (value * ratio_admixture_divergence_0[k]);
                                    admixtureTime[1] = (int) (value * ratio_admixture_divergence_1[k]);
                                    admixtureTime[2] = (int) (value * ratio_admixture_divergence_2[k]);
                                    admixtureProportion = new double[n_way - 1];
                                    admixtureProportion[0] = admixtureProportion_0[l];
                                    admixtureProportion[1] = admixtureProportion_1[l];
                                    admixtureProportion[2] = admixtureProportion_2[l];
                                    divergenceTime0 = j;
                                    divergenceTime1 = element;
                                    divergenceTime2 = item;
                                    divergenceTime3 = value;
                                    demographicModel =
                                            DemographicModelTools.fourWayBuilder_equilibriumPopulation(splitEventTime,
                                                    admixtureProportion, admixtureTime);
                                    demographicModel.trim_default();
                                    demographicModels.add(demographicModel);

                                    // model name
                                    sb.setLength(0);
                                    sb.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("_");
                                    sb.append("way_").append(n_way).append("_").append("adp_").append(admixed_native_introgressed_pop[0]).append("_");
                                    sb.append("nap_").append(admixed_native_introgressed_pop[1]).append("_");
                                    sb.append("in0_").append(admixed_native_introgressed_pop[2]).append("_");
                                    sb.append("in1_").append(admixed_native_introgressed_pop[3]).append("_");
                                    sb.append("in2_").append(admixed_native_introgressed_pop[4]).append("_");
                                    sb.append("ep_").append(equilibriumPopulationSize).append("_");
                                    sb.append("pr0_").append(admixtureProportion_0[l]).append("_");
                                    sb.append("pr1_").append(admixtureProportion_1[l]).append("_");
                                    sb.append("pr2_").append(admixtureProportion_2[l]).append("_");
                                    sb.append("at0_").append(admixtureTime[0]).append("_");
                                    sb.append("at1_").append(admixtureTime[1]).append("_");
                                    sb.append("at2_").append(admixtureTime[2]).append("_");
                                    sb.append("dt0_").append(divergenceTime0).append("_");
                                    sb.append("dt1_").append(divergenceTime1).append("_");
                                    sb.append("dt2_").append(divergenceTime2).append("_");
                                    sb.append("dt3_").append(divergenceTime3);
                                    sb.append(".yaml");
                                    modelName.add(sb.toString());

                                    // simulation metadata
                                    String[] introgressedPop = new String[n_way - 1];
                                    int[] introgressedPop_sampleSize = new int[n_way - 1];
                                    System.arraycopy(admixed_native_introgressed_pop, 2, introgressedPop, 0, n_way - 1);
                                    System.arraycopy(sampleSize, 2, introgressedPop_sampleSize, 0, n_way - 1);
                                    sb_Metadata.setLength(0);
                                    sb_Metadata.append("D").append(PStringUtils.getNDigitNumber(3, count)).append("\t");
                                    sb_Metadata.append(new File(outDir, sb.toString()).getAbsolutePath()).append("\t");
                                    sb_Metadata.append(admixed_native_introgressed_pop[0]).append(":").append(sampleSize[0]).append("\t");
                                    sb_Metadata.append(admixed_native_introgressed_pop[1]).append(":").append(sampleSize[1]).append("\t");
                                    sb_Metadata.append(String.join(",", introgressedPop)).append(":");
                                    sb_Metadata.append(Ints.join(",", introgressedPop_sampleSize)).append("\t");
                                    sb_Metadata.append(sequenceLength).append("\t").append(recombination_rate).append("\t");
                                    sb_Metadata.append(mutation_rate);
                                    bw.write(sb_Metadata.toString());
                                    bw.newLine();
                                    count++;
                                }
                            }
                        }
                    }
                }
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        for (int i = 0; i < demographicModels.size(); i++) {
            DemographicModelTools.writeModel(demographicModels.get(i), new File(outDir, modelName.get(i)));
        }
    }



}
