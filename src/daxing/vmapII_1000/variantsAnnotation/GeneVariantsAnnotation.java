package daxing.vmapII_1000.variantsAnnotation;

import daxing.common.IOTool;

import java.io.File;
import java.util.List;

public class GeneVariantsAnnotation {

    public static void extractVariantsInfo(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");

    }
}
