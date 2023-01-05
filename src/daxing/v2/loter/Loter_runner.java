package daxing.v2.loter;

import com.google.common.collect.Table;

import java.io.File;
import java.util.List;

public class Loter_runner {

    List<File> genotypeList; // simulated genotype list
    String taxaInfoFile;

    /**
     * Taxon\tPopulation\tPopulationName
     * tsk_0\t2\tC
     * tsk_1\t2\tC
     * tsk_2\t2\tC
     */
    Table<String, String, String> taxaInfoMap;

    String logFile; // all log will be append to this logFile

    String outDir; // outDir

    String[] subDir = {"loterGroup", "loterGroupVCF", "LAI"};

    public Loter_runner(String parameterFile){

    }
}
