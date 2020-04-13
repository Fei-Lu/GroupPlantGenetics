package daxing.load;

import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class LoadGO {

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String vmapIIGroupFile, String triadFile,
                          String outDir){
        DateTime.getDateTimeOfNow();
        String[] subdir={"001_count","002_countMerge","003_retainTriad","004_model","005_modelMerge"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles=IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(vmapIIGroupFile, new File(outDir, subdir[0]).getAbsolutePath());
        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxonOutDirMap, e+1));
        EightModelUtils.merge(new File(outDir, subdir[0]).getAbsolutePath(), vmapIIGroupFile, new File(outDir,
                subdir[1]).getAbsolutePath());
        EightModelUtils.retainTriad(triadFile, new File(outDir, subdir[1]).getAbsolutePath(), new File(outDir,
                subdir[2]).getAbsolutePath());
        countGeneNum(new File(outDir, subdir[2]), new File(outDir, subdir[4]));
//        EightModelUtils.countEightModel(new File(outDir, subdir[2]).getAbsolutePath(), vmapIIGroupFile,
//                new File(outDir, subdir[3]).getAbsolutePath());
//        EightModelUtils.mergeModel(new File(outDir, subdir[3]).getAbsolutePath(), new File(outDir, subdir[4]).getAbsolutePath());
//        addGeneNumToModelMergedFile(new File(outDir, subdir[4]+"/geneNum.txt.gz"),
//                new File(outDir, subdir[4]+"/modelMerged.txt.gz"),
//                new File(outDir, subdir[4]+"/modelMergedGeneNum.txt.gz"));
        DateTime.getDateTimeOfNow();
    }

    public static void go(File exonSNPAnnoFile, File exonVCFFile,
                           Map<String, File> taxonOutDirMap, int chr){
        try (BufferedReader br = IOTool.getReader(exonVCFFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            TranscriptDB transcriptDB=new TranscriptDB(exonSNPAnnoFile.getAbsolutePath());
            String[] geneNames=transcriptDB.getGeneName();
            IndividualChrLoad[] taxonLoads=new IndividualChrLoad[taxonNames.size()];
            for (int i = 0; i < taxonLoads.length; i++) {
                taxonLoads[i]=new IndividualChrLoad(taxonNames.get(i), geneNames, chr);
            }
            int pos;
            String geneName, genotype;
            List<String> genotypeList;
            boolean isRefAlleleAncestral=true;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            byte[] indexGenotype; // syn nonsyn del
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
                for (int i = 0; i < taxonNames.size(); i++) {
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
                    if (genotype.equals("0/1")) continue;
                    genotypeByte=GeneLoad.caculateGenotype(genotype, isRefAlleleAncestral);
                    indexGenotype=new byte[2];
                    if (isSyn){
                        indexGenotype[0]=0;
                        indexGenotype[1]=genotypeByte;
                    }else if (isDeleterious){
                        indexGenotype[0]=2;
                        indexGenotype[1]=genotypeByte;
                    } else {
                        indexGenotype[0]=1;
                        indexGenotype[1]=genotypeByte;
                    }
                    taxonLoads[i].addGenotype(geneName, indexGenotype);
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

    private static Map<String, File> getTaxonOutDirMap(String vmapIIGroupFile, String outDir){
        Map<String,File> taxonOutDirMap=new HashMap<>();
        RowTable<String> taxonTable=new RowTable<>(vmapIIGroupFile);
        List<String> taxonNames=taxonTable.getColumn(0);
        File taxonDir;
        for (int i = 0; i < taxonNames.size(); i++) {
            taxonDir=new File(outDir, taxonNames.get(i));
            taxonDir.mkdir();
            taxonOutDirMap.put(taxonNames.get(i), taxonDir);
        }
        return taxonOutDirMap;
    }

    private static void countGeneNum(File inputDir,File outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir.getAbsolutePath());
        RowTable<String> rowTable;
        try {
            BufferedWriter bw=IOTool.getTextGzipWriter(new File(outDir, "geneNum.txt.gz"));
            bw.write("taxon\tgeneNum");
            bw.newLine();
            String taxondName;
            StringBuilder sb;
            for (int i = 0; i < files.size(); i++) {
                rowTable=new RowTableTool<>(files.get(i).getAbsolutePath());
                taxondName=PStringUtils.fastSplit(files.get(i).getName(), ".").get(0);
                sb=new StringBuilder();
                sb.append(taxondName).append("\t").append(rowTable.getRowNumber());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            bw.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void addGeneNumToModelMergedFile(File geneNumFile, File modelMergedFile, File outFile){
        RowTableTool<String> geneNumTable=new RowTableTool<>(geneNumFile.getAbsolutePath());
        Map<String, String> taxonGeneNumMap=geneNumTable.getHashMap(0,1);
        RowTable<String> modelTable=new RowTableTool<>(modelMergedFile.getAbsolutePath());
        List<String> taxonNameList=modelTable.getColumn(4);
        List<String> geneNumList=new ArrayList<>();
        int geneNum;
        for (int i = 0; i < taxonNameList.size(); i++) {
            geneNum=Integer.parseInt(taxonGeneNumMap.get(taxonNameList.get(i)));
            geneNumList.add(String.valueOf(geneNum/3));
        }
        modelTable.addColumn("triadNum", geneNumList);
        modelTable.writeTextTable(outFile.getAbsolutePath(), IOFileFormat.TextGzip);
    }
}