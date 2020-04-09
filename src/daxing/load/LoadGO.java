package daxing.load;

import daxing.common.IOTool;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class LoadGO {

    public static void go(String pgfFile, String exonSNPAnnoDir, String exonVCFDir, String taxaListFile,
                          String outDir){
        TranscriptDB transcriptDB=new TranscriptDB(pgfFile, exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxaListFile, outDir);
        for (int i = 0; i < exonVCFFiles.size(); i++) {
            go(transcriptDB, exonVCFFiles.get(i), taxonOutDirMap, i+1);
        }
    }

    private static void go(TranscriptDB transcriptDB, File exonVCFFile,
                           Map<String, File> taxonOutDirMap, int chr){
        try (BufferedReader br = IOTool.getReader(exonVCFFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            String[] geneNames=transcriptDB.getGeneName(chr);
            IndividualChrLoad[] taxonLoads=new IndividualChrLoad[taxonNames.size()];
            for (int i = 0; i < taxonLoads.length; i++) {
                taxonLoads[i]=new IndividualChrLoad(taxonNames.get(i), geneNames, chr);
            }
            int pos, geneNameIndex;
            String geneName, genotype;
            List<String> genotypeList;
            boolean isRefAlleleAncestral=true;
            boolean isSyn, isNonsyn, isDeleterious;
            GeneLoad geneLoad;
            byte genotypeByte;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                pos=Integer.parseInt(temp.get(1));
                if (!transcriptDB.contain(chr, pos)) continue;
                if (!transcriptDB.hasAncestral(chr, pos)) continue;
                isSyn=transcriptDB.isSyn(chr, pos);
                isNonsyn=transcriptDB.isNonsyn(chr, pos);
                if (isSyn==false && isNonsyn==false) continue;
                isDeleterious=transcriptDB.isDeleterious(chr, pos);
                isRefAlleleAncestral=transcriptDB.isRefAlleleAncestral(chr, pos);
                geneName=transcriptDB.getGeneName(chr, pos);
                geneNameIndex=Arrays.binarySearch(geneNames, geneName);
                for (int i = 0; i < taxonNames.size(); i++) {
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
                    if (genotype.equals("0/1")) continue;
                    geneLoad=new GeneLoad(geneName);
                    genotypeByte=GeneLoad.caculateGenotype(genotype, isRefAlleleAncestral);
                    if (isSyn){
                        geneLoad.addSyn(genotypeByte);
                    }else if (isDeleterious){
                        geneLoad.addHGDeleterious(genotypeByte);
                        geneLoad.addNonsyn(genotypeByte);
                    } else {
                        geneLoad.addNonsyn(genotypeByte);
                    }
                    taxonLoads[i].addGeneLoad(geneNameIndex, geneLoad);
                }
            }
            File outDir;
            for (int i = 0; i < taxonLoads.length; i++) {
                outDir=taxonOutDirMap.get(taxonNames.get(i));
                taxonLoads[i].write(outDir.getAbsolutePath());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Map<String, File> getTaxonOutDirMap(String taxaListFile, String outDir){
        Map<String,File> taxonOutDirMap=new HashMap<>();
        RowTable<String> taxonTable=new RowTable<>(taxaListFile);
        List<String> taxonNames=taxonTable.getColumn(0);
        File taxonDir;
        for (int i = 0; i < taxonNames.size(); i++) {
            taxonDir=new File(outDir, taxonNames.get(i));
            taxonDir.mkdir();
            taxonOutDirMap.put(taxonNames.get(i), taxonDir);
        }
        return taxonOutDirMap;
    }
}
