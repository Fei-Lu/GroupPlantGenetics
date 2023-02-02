package daxing.v2.localAncestryInfer.demography;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;
import daxing.common.utiles.IOTool;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class DemographicModelTools {

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
    private static double[] transform_twoPulse_threeWay_proportion(double[] admixtureProportion){
        assert  admixtureProportion.length ==2 ||  admixtureProportion.length ==1: "error, check admixtureProportion";
        if (admixtureProportion.length ==2 ){
            double[] transformedAdmixtureProportion = new double[2];
            transformedAdmixtureProportion[0] = admixtureProportion[0]/(1-admixtureProportion[1]);
            transformedAdmixtureProportion[1] = admixtureProportion[1];
            return transformedAdmixtureProportion;
        }else if (admixtureProportion.length==1)return admixtureProportion;
        else return new double[0];
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
            epoch = new Epoch(populationSize, epochEndTime[i]);
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

    private static List<Pulse> getEquilibriumPopulationPulses(String[] sourcePop, String destPop,
                                                              double[] admixtureProportion, int[] admixtureTime){
        List<String> sourcePopList = new ArrayList<>();
        for (String pop : sourcePop){
            sourcePopList.add(pop);
        }
        List<Double> admixtureProportionList = new ArrayList<>();
        for (double proportion : DemographicModelTools.transform_twoPulse_threeWay_proportion(admixtureProportion)){
            admixtureProportionList.add(proportion);
        }
        List<Pulse> pulses = new ArrayList<>();
        Pulse pulse;
        for (int i = 0; i < sourcePopList.size(); i++) {
            pulse = new Pulse(sourcePopList, destPop, admixtureProportionList, admixtureTime[i]);
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
     * @return
     */
    public static DemographicModel equilibriumPopulationModelBuilder(String[] demeNames, int[] splitEventTime,
                                                                      String[] ancestors, int populationSize,
                                                                      String[] sourcePop, String destPop,
                                                                      double[] admixtureProportion,
                                                                      int[] admixtureTime){
        int[] epochEndTime = new int[splitEventTime.length * 2 +1];
        System.arraycopy(splitEventTime, 0, epochEndTime, 0, splitEventTime.length);
        List<Deme> demes = DemographicModelTools.getEquilibriumPopulationDemes(demeNames, epochEndTime, ancestors, populationSize);
        List<Pulse> pulses = DemographicModelTools.getEquilibriumPopulationPulses(sourcePop, destPop,
                admixtureProportion, admixtureTime);
        return new DemographicModel.Builder("generations", demes).pulses(pulses).build();
    }

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

    public static void batchRun_twoWay(String simulationMetadataOutFile, String outDir){
//        String[] admixed_native_introgressed_pop = {"E","D","C"};
//        int[] sampleSize = {30,30,30};
//        int equilibriumPopulationSize = 10000;
//        int[] splitEventTime_0 = {10000};
//        int[] splitEventTime_1 = {500, 1000, 2000, 4000, 8000};
//        double[] ratio_admixture_divergence = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64};
//        double[] admixtureProportion = {0.01, 0.02, 0.04, 0.08, 0.16, 0.32};
//        int sequenceLength = 100_000_000;
//        double recombination_rate = 1e-8;


        String[] admixed_native_introgressed_pop = {"E","D","C"};
        int[] sampleSize = {30,30,30};
        int equilibriumPopulationSize = 10000;
        int[] splitEventTime_0 = {10000};
        int[] splitEventTime_1 = {500};
        double[] ratio_admixture_divergence = {0.01};
        double[] admixtureProportion = {0.1, 0.2};
        int sequenceLength = 10_000_000;
        double recombination_rate = 1e-8;


        int[] splitEventTime;
        List<DemographicModel> demographicModels = new ArrayList<>();
        List<String> modelName = new ArrayList<>();
        StringBuilder sb = new StringBuilder();
        DemographicModel demographicModel;
        double admixtureTime, divergenceTime0, divergenceTime1;
        StringBuilder sb_Metadata = new StringBuilder();
        String header = "DemesID\tDemePath\tAdmixedPop\tNativePop\tIntrogressedPopulation\tSequenceLength" +
                "\tRecombination_rate";
        try (BufferedWriter bw = IOTool.getWriter(simulationMetadataOutFile)) {
            bw.write(header);
            bw.newLine();
            int count= 1;
            for (int i = 0; i < splitEventTime_0.length; i++) {
                for (int j = 0; j < splitEventTime_1.length; j++) {
                    for (int k = 0; k < ratio_admixture_divergence.length; k++) {
                        for (int l = 0; l < admixtureProportion.length; l++) {
                            splitEventTime = new int[2];
                            splitEventTime[0] = splitEventTime_0[i];
                            splitEventTime[1] = splitEventTime_1[j];
                            admixtureTime = splitEventTime_1[j] * ratio_admixture_divergence[k];
                            divergenceTime0 = splitEventTime_0[i];
                            divergenceTime1 = splitEventTime_1[j];
                            demographicModel = DemographicModelTools.twoWayBuilder_equilibriumPopulation(splitEventTime,
                                    admixtureProportion[l], (int) admixtureTime);
                            demographicModel.trim_default();
                            demographicModels.add(demographicModel);

                            // model name
                            sb.setLength(0);
                            sb.append("way_").append(2).append("_").append("adp_").append(admixed_native_introgressed_pop[0]).append("_");
                            sb.append("nap_").append(admixed_native_introgressed_pop[1]).append("_");
                            sb.append("inp_").append(admixed_native_introgressed_pop[2]).append("_");
                            sb.append("ep_").append(equilibriumPopulationSize).append("_");
                            sb.append("pr_").append(admixtureProportion[l]).append("_");
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
                            sb_Metadata.append(sequenceLength).append("\t").append(recombination_rate);
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



}