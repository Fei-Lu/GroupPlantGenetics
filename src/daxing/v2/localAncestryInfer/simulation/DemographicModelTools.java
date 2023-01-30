package daxing.v2.localAncestryInfer.simulation;

import com.fasterxml.jackson.annotation.JsonInclude;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.dataformat.yaml.YAMLFactory;
import com.fasterxml.jackson.dataformat.yaml.YAMLGenerator;
import java.io.File;
import java.io.IOException;

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
}
