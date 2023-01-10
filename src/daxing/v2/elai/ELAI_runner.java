package daxing.v2.elai;

import daxing.common.utiles.IOTool;
import daxing.v2.localAncestryInfer.GenotypeTable;
import daxing.v2.localAncestryInfer.TaxaInfo;
import it.unimi.dsi.fastutil.ints.IntList;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.List;

public class ELAI_runner {

    String eLaiPath;
    IntList nWayAdmixture;
    List<String> referencePopList;
    String admixedPop;

    String posFile;

    int expectationMaximizationSteps;

    IntList lowerLayerClusters;

    String preﬁx_outFile;

    public static void prepareInputFile(String genotypeFile, String taxaInfoFile, String admixedPop,
                                        List<String> referencePopList, String outDirBimBam,
                                        String posFile){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
        TaxaInfo taxaInfo = new TaxaInfo(taxaInfoFile);
        String[] allPopOutName = new String[referencePopList.size()+1];
        allPopOutName[0] = admixedPop;
        for (int i = 0; i < referencePopList.size(); i++) {
            allPopOutName[i+1] = referencePopList.get(i);
        }
        int[][] taxaIndex_allPop = new int[allPopOutName.length][];
        List<String> pop_sampleList;
        for (int i = 0; i < taxaIndex_allPop.length; i++) {
            pop_sampleList = taxaInfo.getTaxaListOf(allPopOutName[i]);
            taxaIndex_allPop[i] = new int[pop_sampleList.size()];
            for (int j = 0; j < taxaIndex_allPop[i].length; j++) {
                taxaIndex_allPop[i][j] = genotypeTable.getTaxonIndex(pop_sampleList.get(j));
            }
        }
        GenotypeTable[] genotypeTables = new GenotypeTable[allPopOutName.length];
        for (int i = 0; i < genotypeTables.length; i++) {
            genotypeTables[i] = genotypeTable.getSubGenotypeTableByTaxa(taxaIndex_allPop[i]);
        }
        BufferedWriter[] bwBimBams = new BufferedWriter[allPopOutName.length];
        for (int i = 0; i < bwBimBams.length; i++) {
            bwBimBams[i] = IOTool.getWriter(new File(outDirBimBam, allPopOutName[i]+".inp"));
        }
        try {
            StringBuilder sb = new StringBuilder();
            BufferedWriter bw = IOTool.getWriter(posFile);
            for (int i = 0; i < genotypeTable.getSiteNumber(); i++) {
                sb.setLength(0);
                sb.append(genotypeTable.getSnps()[i].getSnpID()).append(" ");
                sb.append(genotypeTable.getSnps()[i].getPos()).append(" ");
                sb.append(genotypeTable.getSnps()[i].getChr());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
            int taxonIndex;
            char alleleBase;
            String[] subGenotypeTaxa;
            for (int i = 0; i < genotypeTables.length; i++) {
                pop_sampleList = taxaInfo.getTaxaListOf(allPopOutName[i]);
                sb.setLength(0);
                sb.append(pop_sampleList.size()).append(" =  \n");
                sb.append(genotypeTables[i].getSiteNumber()).append("\n");
                subGenotypeTaxa = genotypeTables[i].getTaxa();
                sb.append("IND,");
                sb.append(String.join(",", subGenotypeTaxa));
                bwBimBams[i].write(sb.toString());
                bwBimBams[i].newLine();
                for (int j = 0; j < genotypeTables[i].getSiteNumber(); j++) {
                    sb.setLength(0);
                    sb.append(genotypeTables[i].getSnps()[j].getSnpID()).append(",");
                    for (int k = 0; k < pop_sampleList.size(); k++) {
                        taxonIndex = genotypeTables[i].getTaxonIndex(pop_sampleList.get(k));
                        alleleBase = genotypeTables[i].getAlleleBase(j, taxonIndex);
                        sb.append(alleleBase).append(alleleBase).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bwBimBams[i].write(sb.toString());
                    bwBimBams[i].newLine();
                }
                bwBimBams[i].flush();
                bwBimBams[i].close();
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
