package daxing.temp;

import com.google.common.base.Enums;
import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.GeneLoad;
import daxing.load.ancestralSite.SNPAnnotation;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

public class Go {

    public static void calculateIntrogressionRegionBurden(String introgressionRegionFile, String exonSNPAnnoDir,
                                                          String exonVCFDir,
                                                           String taxa_InfoDBFile, String outDir,
                                                          SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        String[] outNames=
                exonVCFFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",
                        "_introgressionRegion.txt.gz")).toArray(String[]::new);
        IntStream.range(0, exonVCFFiles.size()).forEach(e->calculateIntrogressionRegionBurden(introgressionRegionFile, exonAnnoFiles.get(e),
                exonVCFFiles.get(e), new File(outDir, outNames[e]), e+1, taxa_InfoDBFile, methodCallDeleterious));
    }

    private static void calculateIntrogressionRegionBurden(String introgressionRegionFile, File exonSNPAnnoFile,
                                                          File exonVCFFile,
                                                          File outFile, int chr,
                                                          String taxaInfoFile,
                                                          SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
        Map<String, String> taxaMap= RowTableTool.getMap(taxaInfoFile, 0, 36);
        Map<String, String> taxaP2Map=new HashMap<>();
        for (Map.Entry<String, String> entry: taxaMap.entrySet()){
            if(Enums.getIfPresent(IntrogressionRegion.P2.class, entry.getValue()).isPresent()){
                taxaP2Map.put(entry.getKey(), entry.getValue());
            }
        }
        IntrogressionRegion introgressionRegion = new IntrogressionRegion(introgressionRegionFile);
        try (BufferedReader br = IOTool.getReader(exonVCFFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            ChrSNPAnnoDB transcriptDB=new ChrSNPAnnoDB(exonSNPAnnoFile);
            IndividualBurden[] taxonBurden=new IndividualBurden[taxonNames.size()];
            for (int i = 0; i < taxonBurden.length; i++) {
                taxonBurden[i]=new IndividualBurden(taxonNames.get(i), 10);
            }
            int pos, depth;
            String genotype, subLine;
            List<String> genotypeList, depthList;
            boolean isRefAlleleAncestral;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            byte[] indexGenotype; // syn nonsyn del
            String sub;
            IntrogressionRegion.P2 p2;
            WheatLineage wheatLineage;
            double[] p3Bins;
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
                sub=transcriptDB.getGeneName(chr, pos).substring(8,9);
                wheatLineage=WheatLineage.valueOf(sub);
                temp = PStringUtils.fastSplit(line);
                for (int i = 0; i < taxonNames.size(); i++) {
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
//                    double corrRatio=chrPosCorrRatioMap.get(chr, pos);
                    if (isSyn){
                        indexGenotype[0]=0;
                    }else if (isDeleterious){
                        indexGenotype[0]=2;
                    } else {
                        indexGenotype[0]=1;
                    }
                    //                    taxonLoads[i].addGenotype(geneName, indexGenotype, corrRatio);
//                    taxonLoads[i].addGenotype(geneName, indexGenotype, methodCallDeleterious);
                    p2=IntrogressionRegion.P2.valueOf(taxaP2Map.get(taxonNames.get(i)));
                    p3Bins = introgressionRegion.getP3fdModifyBins(chr,pos,p2);
                    for (IntrogressionRegion.P3 p3: IntrogressionRegion.P3.values()){
                        taxonBurden[i].addLoadCount(wheatLineage,p3Bins[p3.index],p3, indexGenotype);
                    }
                }
            }
            // dim1: Sub, dim2: fdBins, dim3: p3, dim4: burden 9 column
            bw.write("Taxa\tSub\tFdModifyBins\tP3\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn" +
                    "\tnumDerivedInNonsyn\tnumHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious" +
                    "\tnumHeterInHGDeleterious");
            bw.newLine();
            for (IndividualBurden individualBurden:taxonBurden){
                bw.write(individualBurden.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
