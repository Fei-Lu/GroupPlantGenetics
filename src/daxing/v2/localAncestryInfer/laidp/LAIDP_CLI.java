package daxing.v2.localAncestryInfer.laidp;

import org.apache.commons.cli.*;
import pgl.infra.utils.Benchmark;
import java.util.BitSet;

public class LAIDP_CLI {


    /**
     * Required parameters
     */
    String genotypeFile;
    String taxaGroupFile;
    String ancestralAllele;
    String localAnceOutFile;


    /**
     * Optional parameters
     */
    int windowSize;
    int stepSize;
    int conjunctionNum;
    double switchCostScore;

    int maxSolutionCount;

    int threadsNum;

    public LAIDP_CLI(CommandLine line){
        this.genotypeFile=line.getOptionValue("g");
        if (line.hasOption("w")){
            this.windowSize=Integer.parseInt(line.getOptionValue("w"));
        }else {
            this.windowSize=200;
        }
        if (line.hasOption("s")){
            this.stepSize = Integer.parseInt(line.getOptionValue("s"));
        }else {
            this.stepSize = 100;
        }
        this.taxaGroupFile = line.getOptionValue("taxaGroup");
        this.ancestralAllele = line.getOptionValue("ancestral");
        if (line.hasOption("conjunctionNum")){
            this.conjunctionNum = Integer.parseInt(line.getOptionValue("conjunctionNum"));
        }else {
            this.conjunctionNum= 2;
        }
        if (line.hasOption("switchCost")){
            this.switchCostScore=Double.parseDouble(line.getOptionValue("switchCost"));
        }else {
            this.switchCostScore = 1.5;
        }
        if (line.hasOption("maxSolution")){
            this.maxSolutionCount = Integer.parseInt(line.getOptionValue("maxSolution"));
        }else {
            this.maxSolutionCount = 32;
        }
        this.localAnceOutFile = line.getOptionValue("out");
        if (line.hasOption("t")){
            this.threadsNum = Integer.parseInt(line.getOptionValue("t"));
        }else {
            this.threadsNum = 1;
        }
    }

    public void startRun(){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
        BitSet[] ancestralAlleleBitSet = genotypeTable.getAncestralAlleleBitSet(ancestralAllele);
        BitSet[][] localAnc = genotypeTable.calculateLocalAncestry(windowSize, stepSize, taxaGroupFile,
                ancestralAlleleBitSet, conjunctionNum, switchCostScore, maxSolutionCount, threadsNum);
        genotypeTable.write_localAncestry(localAnc, localAnceOutFile, taxaGroupFile);
    }

    private static Options SetOptions(){
        Option genotypeFile = new Option("g", "genotypeFile",true,"<fileName>, path of genotype file in vcf format");
        Option windowSize = new Option("w", "windowSize", true, "<integer> window size in variants");
        Option stepSize = new Option("s", "stepSize", true, "<integer> step size in variants");
        Option taxaGroupFile = new Option("taxaGroup","taxaGroupFile", true, "<fileName>, path of taxaGroup" +
                " file");
        Option ancestalyAllele = new Option("ancestral", "ancestralAllele", true, "<integer|fileName> ancestral " +
                "allele, " +
                "integer represent outGroup sample index, separated by comma when multiple outGoups exists; when it " +
                "is string, ancestral allele file is required");
        Option conjunctionNum = new Option("conjunctionNum", "conjunctionNum", true, "the number of grid when define local " +
                "introgressed interval");
        Option switchCostScore = new Option("switchCost", "switchCostScore",true, "switch cost score when haplotype" +
                " " +
                "switching occurs");
        Option maxSolutionCount = new Option("maxSolution","maxSolutionCount",true, "upper bound of the candidate solution");
        Option localAncestyOutFile = new Option("out", "outFile",true, "prefix of outfile");
        Option threadsNum = new Option("t", "threadsNum",true, "thread number");
        Options options = new Options();
        options.addOption(genotypeFile);
        options.addOption(windowSize);
        options.addOption(stepSize);
        options.addOption(taxaGroupFile);
        options.addOption(ancestalyAllele);
        options.addOption(conjunctionNum);
        options.addOption(switchCostScore);
        options.addOption(maxSolutionCount);
        options.addOption(localAncestyOutFile);
        options.addOption(threadsNum);
        return options;
    }

    private static String getCmdLineSyntax(){
        return "java -jar LAIDP.jar -g genotypeFile" + " " +
                "-taxaGroup taxaGroupFile" + " " +
                "-ancestral ancestral" + " " +
                "-out outFile";
    }

    private static String getHeader(){
        return "LAIDP v1.0";
    }

    private static String getFooter(){
        return "Author: Daxing Xu, email: dxxu@genetics.ac.cn" + "\n" +
                "See below for detailed documentation" + "\n" +
                "https://github.com/";
    }

    public static void startFromCLI(String[] args){
        Options options = LAIDP_CLI.SetOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine line;
        HelpFormatter helpFormatter = new HelpFormatter();
        if (args.length < 1) {
            helpFormatter.printHelp(getCmdLineSyntax(), getHeader(), options,getFooter());
            System.exit(1);
        }
        long start = System.nanoTime();
        try {
            line=parser.parse(options, args);
            new LAIDP_CLI(line).startRun();
        } catch (ParseException e) {
            System.err.println("Parsing failed.  Reason: " + e.getMessage());
            helpFormatter.printHelp(getCmdLineSyntax(), getHeader(), options, getFooter());
            System.exit(1);
        }
        System.out.println("Completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        System.out.println("Done");
    }
}
