package daxing.v2.localAncestryInfer.simulation;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;
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
     * @param demeNames String[5]
     * @param divergenceTime int[2], the first dim is divergence time between agent population of introgressed ancestry and
     *                       admixed population, the second dim is divergence time between agent population of
     *                       native ancestry and admixed population; such as 10000, 1000
     * @param ancestors only one ancestor per population
     * @param populationSize fixed
     * @param sourcePop only one source population per admixed population
     * @param desPop desPop
     * @param admixtureProportion admixtureProportion
     * @param admixtureTime admixtureTime
     * @return
     */
    private static DemographicModel equilibriumPopulation_onePulse_twoWay_Builder(String[] demeNames, int[] divergenceTime,
                                                                String[] ancestors, int populationSize,
                                                                String sourcePop, String desPop,
                                                                double admixtureProportion,
                                                                int admixtureTime){
        int[] epochEndTime = {divergenceTime[0], divergenceTime[1], 0, 0, 0};
        List<String> sourcePopList = new ArrayList<>();
        sourcePopList.add(sourcePop);
        List<Double> admixtureProportionList = new ArrayList<>();
        admixtureProportionList.add(admixtureProportion);
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
        List<Pulse> pulses = new ArrayList<>();
        Pulse pulse;
        for (int i = 0; i < sourcePopList.size(); i++) {
            pulse = new Pulse(sourcePopList, desPop, admixtureProportionList, admixtureTime);
            pulses.add(pulse);
        }
        return new DemographicModel.Builder("generations", demes).pulses(pulses).build();
    }

    public static DemographicModel equilibriumPopulation_onePulse_twoWay_Builder(int[] divergenceTime,
                                                                                 double admixtureProportion,
                                                                                 int admixtureTime){
        String[] demeNames = {"A","B","C","D","E"};
        String[] ancestors = {null,"A","A","B","B"};
        int populationSize = 10000;
        return DemographicModelTools.equilibriumPopulation_onePulse_twoWay_Builder(demeNames, divergenceTime, ancestors,
                populationSize, "C", "E", admixtureProportion, admixtureTime);
    }

}
