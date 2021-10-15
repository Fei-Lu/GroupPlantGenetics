package daxing.load.esg;

import daxing.common.utiles.IOTool;
import pgl.infra.utils.Benchmark;

import java.io.BufferedWriter;
import java.io.IOException;
import java.util.List;

public class Go {

    public static void start(){

        String allRatioOEFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/006_vmap2_1062/007_complementary_IndiviExpectation/004_OE/002_RatioOE_triasdsBlock_merge/all.RatioOE.txt.gz";
        String allTaxaESG="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/006_vmap2_1062/007_complementary_IndiviExpectation/004_OE/004_ESG_perTaxon_merge/All.ESG.txt.gz";
        String outFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/006_vmap2_1062" +
                "/007_complementary_IndiviExpectation/004_OE/res.txt";
        TaxaOERatio taxaOERatio = TaxaOERatio.getInstance(allRatioOEFile);
        AllTaxaESG taxaESG= AllTaxaESG.getInstance(allTaxaESG);
        List<String> nonNATaxaList;
        long start=System.nanoTime();
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            StringBuilder sb =new StringBuilder();
            bw.write("TriadsBlockID\tIfPseudo\tNonNATaxaNum\tESG\tExpected\tObserved");
            bw.newLine();
            String triadsBlockID;
            int ifPseudoValue;
            int nonNATaxaNum;
            double[] expectedCount;
            int[] observedCount;
            for (int i = 0; i < taxaESG.getTriadsNum(); i++) {
                for (IfPseudoHexaploid ifPseudoHexaploid : IfPseudoHexaploid.values()){
                    nonNATaxaList = taxaESG.getNonNATaxaList(i, ifPseudoHexaploid);
                    expectedCount = taxaOERatio.getESGExpectedCount(nonNATaxaList);
                    observedCount = taxaESG.getESGObservedCount(i, ifPseudoHexaploid);
                    for (ESG esg : ESG.values()){
                        sb.setLength(0);
                        triadsBlockID = taxaESG.getTriadsBlockID(i);
                        ifPseudoValue = ifPseudoHexaploid.getValue();
                        nonNATaxaNum = taxaESG.getNonNATaxaNum(i, ifPseudoHexaploid);
                        sb.append(triadsBlockID).append("\t");
                        sb.append(ifPseudoValue).append("\t");
                        if (nonNATaxaNum == 0){
                            sb.append(nonNATaxaNum).append("\t");
                            sb.append(esg.name()).append("\t");
                            sb.append("NA").append("\t");
                            sb.append("NA");
                            bw.write(sb.toString());
                            bw.newLine();
                        }else {
                            sb.append(nonNATaxaNum).append("\t");
                            sb.append(esg.name()).append("\t");
                            sb.append(expectedCount[esg.getValue()]).append("\t");
                            sb.append(observedCount[esg.getValue()]);
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                }
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Writing spend "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
    }
}
