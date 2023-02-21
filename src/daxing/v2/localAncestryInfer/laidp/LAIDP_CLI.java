package daxing.v2.localAncestryInfer.laidp;

import daxing.common.utiles.IOTool;
import org.apache.commons.cli.*;
import java.io.BufferedWriter;
import java.io.IOException;
import java.util.BitSet;
import java.util.List;

public class LAIDP_CLI {

    String genotypeFile;
    int windowSize;
    int stepSize;
    String admixedPopulation;
    String nativePopulation;
    String[] introgressedPopulation;
    String ancestralAllele;
    int conjunctionNum;
    double switchCostScore;
    String localAnceOutFile;
    int threadsNum;

    public LAIDP_CLI(CommandLine line){
        this.genotypeFile=line.getOptionValue("g");
        this.windowSize = Integer.parseInt(line.getOptionValue("w"));
        this.stepSize = Integer.parseInt(line.getOptionValue("s"));
        this.admixedPopulation = line.getOptionValue("a");
        this.nativePopulation = line.getOptionValue("n");
        this.introgressedPopulation = line.getOptionValues("i");
        this.ancestralAllele = line.getOptionValue("ancestral");
        this.conjunctionNum = Integer.parseInt(line.getOptionValue("conjunctionNum"));
        if (line.hasOption("switchCost")){
            this.switchCostScore = Double.parseDouble(line.getOptionValue("switchCost"));
        }else {
            this.switchCostScore = 1.5;
        }
        this.localAnceOutFile = line.getOptionValue("out");
        this.threadsNum = Integer.parseInt(line.getOptionValue("t"));
    }

//    public void startRun(){
//        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
//        BitSet[] ancestralAlleleBitSet = genotypeTable.getAncestralAlleleBitSet(ancestralAllele);
//        BitSet[][] localAnc = genotypeTable.calculateLocalAncestry(windowSize, stepSize, taxaGroupFile,
//                ancestralAlleleBitSet, conjunctionNum, switchCostScore, threadsNum);
//        int variantsNum = genotypeTable.getSiteNumber();
//        LAIDP_CLI.write_localAncestry(localAnc, localAnceOutFile, variantsNum, taxaGroupFile);
//    }

    private static CommandLine buildFromCLI(String[] args){
        Option genotypeFile = new Option("g", "genotypeFile",true,"<fileName>, path of genotype file in vcf format");
        Option windowSize = new Option("w", "windowSize", true, "<integer> window size in variants");
        Option stepSize = new Option("s", "stepSize", true, "<integer> step size in variants");
        Option admixedPop = new Option("a", "admixed", true, "<string> admixed population");
        Option nativePopulation = new Option("n", "native", true, "<string> native population");
        Option introgressedPopulation = new Option("i", "introgressed", true, "<string> introgressed populations, " +
                "when " +
                "multiple introgressed populations exists, use a comma separator");
        Option ancestalyAllele = new Option("ancestral", "ancestralAllele", true, "<integer|fileName> ancestral " +
                "allele, " +
                "integer represent outGroup sample, separated by comma when multiple outGoups exists; when it " +
                "is string, ancestral allele file is required");
        Option conjunctionNum = new Option("conjunctionNum", "conjunctionNum", true, "the number of grid when define local " +
                "introgressed interval");
        Option switchCostScore = new Option("switchCost", "switchCostScore",true, "switch cost score when haplotype" +
                " " +
                "switching occurs");
        Option localAncestyOutFile = new Option("out", "outFile",true, "prefix of outfile");
        Option threadsNum = new Option("t", "threadsNum",true, "thread number");
        Options options = new Options();
        options.addOption(genotypeFile);
        options.addOption(windowSize);
        options.addOption(stepSize);
        options.addOption(admixedPop);
        options.addOption(nativePopulation);
        options.addOption(introgressedPopulation);
        options.addOption(ancestalyAllele);
        options.addOption(conjunctionNum);
        options.addOption(switchCostScore);
        options.addOption(localAncestyOutFile);
        options.addOption(threadsNum);
        CommandLineParser parser = new DefaultParser();
        CommandLine line = null;
        HelpFormatter helpFormatter = new HelpFormatter();
        try {
            line=parser.parse(options, args);
        } catch (ParseException e) {
            System.err.println("Parsing failed.  Reason: " + e.getMessage());
            helpFormatter.printHelp("LAIDP", options);
            System.exit(1);
        }
        LAIDP_CLI laidpCli = new LAIDP_CLI(line);

        return line;
    }

    public static void write_localAncestry(BitSet[][] localAncestry, String localAncestryOutFile, int variantsNum,
                                           String taxaGroupFile){
        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);
        try (BufferedWriter bw = IOTool.getWriter(localAncestryOutFile)) {
            List<String> admixedTaxaList = taxaGroup.getTaxaOf(Source.ADMIXED);
            StringBuilder sb = new StringBuilder();
            sb.append(String.join("\t", admixedTaxaList));
            bw.write(sb.toString());
            bw.newLine();
            int admixedTaxonNum = localAncestry.length;
            int sourceNum = localAncestry[0].length;
            int ancestry;
            for (int variantIndex = 0; variantIndex < variantsNum; variantIndex++) {
                sb.setLength(0);
                for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxonNum; admixedTaxonIndex++) {
                    for (int sourceIndex = 0; sourceIndex < sourceNum; sourceIndex++) {
                        ancestry = localAncestry[admixedTaxonIndex][sourceIndex].get(variantIndex) ? 1 : 0;
                        sb.append(ancestry).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    sb.append("\t");
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
}
