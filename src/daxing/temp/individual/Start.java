package daxing.temp.individual;

import daxing.common.table.RowTableTool;
import daxing.common.utiles.DateTime;
import daxing.common.utiles.IOTool;
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
import java.util.stream.IntStream;

public class Start {

    public static void calculateDonorBurdenPerIndividual(String individualFdDir, String exonSNPAnnoDir,
                                                         String exonVCFDir, String outDir, String taxaInfoFile,
                                                         SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        String[] outNames= exonVCFFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",
                        "_IndividualFdDonor.txt.gz")).toArray(String[]::new);
        IndividualFd[] individualFds=getIndividualFd(individualFdDir, taxaInfoFile);
        IntStream.range(0, exonVCFFiles.size()).forEach(e-> calculateDonorBurdenPerIndividual(individualFds,
                exonAnnoFiles.get(e), exonVCFFiles.get(e), new File(outDir, outNames[e]), e+1 ,methodCallDeleterious));
    }

    private static void calculateDonorBurdenPerIndividual(IndividualFd[] individualFds, File exonSNPAnnoFile,
                                                          File exonVCFFile,
                                                          File outFile, int chr,
                                                          SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        long start=System.nanoTime();
        List<String> fdTaxaList = new ArrayList<>();
        for (IndividualFd individualFd : individualFds){
            fdTaxaList.add(individualFd.taxon);
        }
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
            String genotype, subLine, refChr;
            List<String> genotypeList, depthList;
            boolean isRefAlleleAncestral;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            byte[] indexGenotype; // syn nonsyn del
            EnumSet<Donor> donorEnumSet;
            IndividualFd individualFd;
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
                    if (!fdTaxaList.contains(taxonNames.get(i))) continue;
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
                    individualFd=new IndividualFd(taxonNames.get(i));
                    int index = Arrays.binarySearch(individualFds,individualFd);
                    refChr = RefV1Utils.getChromosome(chr, pos);
                    refPos = RefV1Utils.getPosOnChromosome(chr, pos);
                    donorEnumSet=individualFds[index].getDonorFrom(refChr, refPos);
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
                if (!fdTaxaList.contains(individualBurden.taxon)) continue;
                bw.write(individualBurden.toString());
                bw.newLine();
            }
            bw.flush();
            System.out.println(outFile.getName()+ " completed in "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static IndividualFd[] getIndividualFd(String individualFdDir, String taxaInfo){
        List<File> individualFdFiles = IOTool.getFileListInDirEndsWith(individualFdDir, ".gz");
        Map<String, String> introgressionIDTaxaMap= RowTableTool.getMap(taxaInfo, 35, 0);
        IndividualFd[] individualFds = new IndividualFd[individualFdFiles.size()];
        File file;
        String taxon;
        for (int i = 0; i < individualFds.length; i++) {
            file = individualFdFiles.get(i);
            taxon = individualFdFiles.get(i).getName().substring(13,20);
            individualFds[i] = new IndividualFd(file.getAbsolutePath(), introgressionIDTaxaMap.get(taxon));
        }
        return individualFds;
    }
}
