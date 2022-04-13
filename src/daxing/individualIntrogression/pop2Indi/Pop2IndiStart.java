package daxing.individualIntrogression.pop2Indi;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.DateTime;
import daxing.common.utiles.IOTool;
import daxing.individualIntrogression.individual.Donor;
import daxing.individualIntrogression.individual.IndividualBurden_individualFd;
import daxing.individualIntrogression.window.PopulationIndividualFdStart;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.GeneLoad;
import daxing.load.ancestralSite.SNPAnnotation;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Pop2IndiStart {

    public static void start(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile,
                             SNPAnnotation.MethodCallDeleterious methodCallDeleterious, String popFdFile,
                             String individualFdDir,
                             String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subDirs = {"001_window","002_individual"};
        for (String subDir : subDirs){
            new File(outDir, subDir).mkdir();
        }
        calculateWindowBasedBurden(exonSNPAnnoDir, exonVCFDir, taxa_InfoDBFile, methodCallDeleterious, popFdFile,
                individualFdDir, new File(outDir, subDirs[0]).getAbsolutePath());
        calculateDonorBurdenPerIndividual(individualFdDir, exonSNPAnnoDir, exonVCFDir, new File(outDir, subDirs[1]).getAbsolutePath(),
                taxa_InfoDBFile, methodCallDeleterious);
        System.out.println(DateTime.getDateTimeOfNow());
    }

    /**
     * window based burden
     * @param exonSNPAnnoDir
     * @param exonVCFDir
     * @param taxa_InfoDBFile
     * @param methodCallDeleterious
     * @param popFdFile
     * @param individualFdDir
     * @param outDir
     */
    public static void calculateWindowBasedBurden(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile,
                             SNPAnnotation.MethodCallDeleterious methodCallDeleterious, String popFdFile,
                             String individualFdDir,
                             String outDir){
        long start = System.nanoTime();
        PopulationIndividualFdStart.start(exonSNPAnnoDir, exonVCFDir, taxa_InfoDBFile, methodCallDeleterious,
                popFdFile, individualFdDir, outDir);
        System.out.println("Window based burden completed in "+Benchmark.getTimeSpanSeconds(start)+ " s");
    }

    /**
     * individual based burden
     * @param individualFdDir
     * @param exonSNPAnnoDir
     * @param exonVCFDir
     * @param outDir
     * @param taxaInfoFile
     * @param methodCallDeleterious
     */
    public static void calculateDonorBurdenPerIndividual(String individualFdDir, String exonSNPAnnoDir,
                                                         String exonVCFDir, String outDir, String taxaInfoFile,
                                                         SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        System.out.println(DateTime.getDateTimeOfNow());
        long start = System.nanoTime();
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        String[] outNames= exonVCFFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",
                "_IndividualFdDonor.txt.gz")).toArray(String[]::new);
        IndividualsFdDxyDonor individualsFdDxyDonor = new IndividualsFdDxyDonor(individualFdDir);
        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e-> calculateDonorBurdenPerIndividual(individualsFdDxyDonor,
                exonAnnoFiles.get(e), exonVCFFiles.get(e), new File(outDir, outNames[e]), e+1, taxaInfoFile,
                methodCallDeleterious));
        System.out.println("Individual based burden completed in "+ Benchmark.getTimeSpanSeconds(start)+ " s");
    }

    private static void calculateDonorBurdenPerIndividual(IndividualsFdDxyDonor individualsFdDxyDonor, File exonSNPAnnoFile,
                                                          File exonVCFFile,
                                                          File outFile, int chr, String taxaInfoFile,
                                                          SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        long start=System.nanoTime();
        Map<String, String> introgressionID2vcfIDMap = RowTableTool.getMap(taxaInfoFile, 35, 0);
        Map<String, String> vcfID2IntrogressionIDMap = RowTableTool.getMap(taxaInfoFile, 0, 35);
        List<String> fdVCFIDList = individualsFdDxyDonor.getTaxa().stream().map(s -> introgressionID2vcfIDMap.get(s)).collect(Collectors.toList());
        try (BufferedReader br = IOTool.getReader(exonVCFFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            ChrSNPAnnoDB transcriptDB=new ChrSNPAnnoDB(exonSNPAnnoFile);
            IndividualBurden_individualFd[] taxonBurden=new IndividualBurden_individualFd[taxonNames.size()];
            for (int i = 0; i < taxonBurden.length; i++) {
                taxonBurden[i]=new IndividualBurden_individualFd(taxonNames.get(i), chr);
            }
            int pos, depth, refPos;
            String genotype, subLine, refChr, introgressionID;
            List<String> genotypeList, depthList;
            boolean isRefAlleleAncestral;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            byte[] indexGenotype; // syn nonsyn del
            EnumSet<Donor> donorEnumSet;
            while ((line=br.readLine())!=null){
                subLine= line.substring(0,20);
                temp=PStringUtils.fastSplit(subLine);
                pos=Integer.parseInt(temp.get(1));
                if (!transcriptDB.contain(chr, pos)) continue;
                if (!transcriptDB.hasAncestral(chr, pos)) continue;
                isSyn=transcriptDB.isSyn(chr, pos);
                isNonsyn=transcriptDB.isNonsyn(chr, pos);
                isDeleterious=transcriptDB.isDeleterious(chr, pos, methodCallDeleterious);
                if (!(isSyn || isNonsyn || isDeleterious)) continue;
                isRefAlleleAncestral=transcriptDB.isRefAlleleAncestral(chr, pos);
                temp = PStringUtils.fastSplit(line);
                for (int i = 0; i < taxonNames.size(); i++) {
                    if (!fdVCFIDList.contains(taxonNames.get(i))) continue;
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
//                    if (genotype.equals("0/1")) continue;
                    depthList=PStringUtils.fastSplit(genotypeList.get(1),",");
                    depth=Integer.parseInt(depthList.get(0))+Integer.parseInt(depthList.get(1));
                    if (depth < 2) continue;
                    genotypeByte= GeneLoad.caculateGenotype(genotype, isRefAlleleAncestral);
                    indexGenotype=new byte[2];
                    indexGenotype[1]=genotypeByte;
                    if (isSyn){
                        indexGenotype[0]=0;
                    }else if (isDeleterious){
                        indexGenotype[0]=2;
                    } else {
                        indexGenotype[0]=1;
                    }
                    introgressionID = vcfID2IntrogressionIDMap.get(taxonNames.get(i));
                    refChr = RefV1Utils.getChromosome(chr, pos);
                    refPos = RefV1Utils.getPosOnChromosome(chr, pos);
                    donorEnumSet = individualsFdDxyDonor.getDonorFrom(introgressionID,refChr, refPos);
                    for (Donor donor: donorEnumSet){
                        taxonBurden[i].addLoadCount(donor, indexGenotype);
                    }
                }
            }
            bw.write("Taxa\tChrID\tDonor\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn" +
                    "\tnumDerivedInNonsyn\tnumHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious" +
                    "\tnumHeterInHGDeleterious");
            bw.newLine();
            for (IndividualBurden_individualFd individualBurden:taxonBurden){
                if (!fdVCFIDList.contains(individualBurden.getTaxon())) continue;
                bw.write(individualBurden.toString());
                bw.newLine();
            }
            bw.flush();
            System.out.println(outFile.getName()+ " completed in "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
