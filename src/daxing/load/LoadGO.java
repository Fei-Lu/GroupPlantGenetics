package daxing.load;

import daxing.common.*;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class LoadGO {

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/003_exonSNPAnno/001_byChrID";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/002_exonSNPVCF";
        String vmapIIGroupFile="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/vmapGroup.txt";
        String triadFile="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/triadGenes1.1_cdsLen_geneHC.txt";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/002_analysis/014_deleterious/004_analysis/003_geneLoadIndividual_GeneHC";
        go(exonSNPAnnoDir, exonVCFDir, vmapIIGroupFile, triadFile, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String vmapIIGroupFile, String triadFile, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_retainTriad","004_standardization"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles=IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(vmapIIGroupFile, new File(outDir, subdir[0]).getAbsolutePath());
        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxonOutDirMap, e+1));
        merge(new File(outDir, subdir[0]).getAbsolutePath(), vmapIIGroupFile, new File(outDir, subdir[1]).getAbsolutePath());
        retainTriad(triadFile, new File(outDir, subdir[1]).getAbsolutePath(), new File(outDir, subdir[2]).getAbsolutePath());
        normalizedTriadByAncestralNum(new File(outDir, subdir[2]).getAbsolutePath(), new File(outDir, subdir[3]).getAbsolutePath());
    }

    private static void go(File exonSNPAnnoFile, File exonVCFFile,
                           Map<String, File> taxonOutDirMap, int chr){
        try (BufferedReader br = IOTool.getReader(exonVCFFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){}
            temp= PStringUtils.fastSplit(line);
            List<String> taxonNames=temp.subList(9, temp.size());
            ChrSNPAnnoDB transcriptDB=new ChrSNPAnnoDB(exonSNPAnnoFile);
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

    public static void merge(String inputDir, String vmapIIGroupFile, String outDir){
        List<File> dirs=IOTool.getVisibleDir(inputDir);
        RowTableTool<String> table=new RowTableTool<>(vmapIIGroupFile);
        Map<String,String> taxonGroupMap=table.getHashMap(0, 15);
        Predicate<File> p= f->taxonGroupMap.get(f.getName()).equals("Landrace_Europe") || taxonGroupMap.get(f.getName()).equals("Cultivar");
        List<File> ldCl=dirs.stream().filter(p).collect(Collectors.toList());
        List<File> temp;
        RowTableTool<String> tableTool;
        String outName;
        for (int i = 0; i < ldCl.size(); i++) {
            temp= IOUtils.getVisibleFileListInDir(ldCl.get(i).getAbsolutePath());
            tableTool=new RowTableTool<>(temp.get(0).getAbsolutePath());
            outName= PStringUtils.fastSplit(temp.get(0).getName(), ".").get(1);
            for (int j = 1; j < temp.size(); j++) {
                tableTool.add(new RowTableTool<>(temp.get(j).getAbsolutePath()));
            }
            tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
        }
    }

    public static void retainTriad(String triadFile, String inputDir, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        IntStream.range(0, files.size()).parallel().forEach(e->retainTriad(triadFile, files.get(e), outDir));
    }

    private static void retainTriad(String triadFile, File inputFile, String outDir){
        Triad triad=new Triad(triadFile);
        RowTableTool<String> table=new RowTableTool<>(inputFile.getAbsolutePath());
        List<String> geneNames=table.getColumn(0);
        Collections.sort(geneNames);
        String[] threeName;
        int indexABD[]=new int[3];
        String triadID, taxon, temp;
        taxon=PStringUtils.fastSplit(inputFile.getName(), ".").get(0);
        StringBuilder sb;
        BufferedWriter bw;
        int[] cdsLen;
        try {
            bw=IOTool.getTextGzipWriter(new File(outDir, taxon + ".triad.txt.gz"));
            bw.write("TriadID\tcdsLen\tGeneName\tnumSyn\tnumDerivedInSyn\tnumNonsyn\tnumDerivedInNonsyn" +
                    "\tnumHGDeleterious\tnumDerivedInHGDeleterious\tsubgenome");
            bw.newLine();
            String subgenome, geneName;
            for (int i = 0; i < triad.getRowNum(); i++) {
                threeName=triad.getTriad(i);
                for (int j = 0; j < indexABD.length; j++) {
                    indexABD[j]= Collections.binarySearch(geneNames, threeName[j]);
                }
                if (indexABD[0] < 0) continue;
                if (indexABD[1] < 0) continue;
                if (indexABD[2] < 0) continue;
                triadID=triad.getTraidID(i);
                cdsLen=triad.getCDSLen(i);
                for (int j = 0; j < indexABD.length; j++) {
                    temp=String.join("\t", table.getRow(indexABD[j]));
                    sb=new StringBuilder();
                    sb.append(triadID).append("\t").append(cdsLen[j]).append("\t");
                    sb.append(temp).append("\t");
                    geneName=table.getRow(indexABD[j]).get(0);
                    subgenome=geneName.substring(8,9);
                    sb.append(subgenome);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            bw.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void normalizedTriadByAncestralNum(String inputDir,String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        IntStream.range(0, files.size()).forEach(f->normalizedTriadByAncestralNum(files.get(f), outDir));
    }

    public static void normalizedTriadByAncestralNum(File inputFile, String outDir){
        String fileName=inputFile.getName().replaceAll("triad", "triad.normalized");
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getTextGzipWriter(new File(outDir, fileName))) {
            bw.write("TriadID\tnormalizedNumDerivedInSynA" +
                    "\tnormalizedNumDerivedInSynB\tnormalizedNumDerivedInSynD\tsynRegion" +
                    "\tnormalizedNumDerivedInNonsynA" +
                    "\tnormalizedNumDerivedInNonsynB\tnormalizedNumDerivedInNonsynD\tnonsynRegion" +
                    "\tnormalizedNumDerivedInHGDeleteriousA\tnormalizedNumDerivedInHGDeleteriousB" +
                    "\tnormalizedNumDerivedInHGDeleteriousD\tdelRegion");
            bw.newLine();
            br.readLine();
            String line, triadID, region;
            List<String>[] temp;
            int[] cdsLen;
            int[] snpNum;
            double[] derivedSyn;
            double[] derivedNonsyn;
            double[] deleterious;
            StringBuilder sb;
            double[] normalizedDerivedSyn;
            double[] normalizedDerivedNonSyn;
            double[] normalizedDerivedDel;
            int index;
            while ((line=br.readLine())!=null){
                temp=new List[3];
                temp[0]=PStringUtils.fastSplit(line);
                temp[1]=PStringUtils.fastSplit(br.readLine());
                temp[2]=PStringUtils.fastSplit(br.readLine());
                cdsLen=new int[3];
                snpNum=new int[3];
                derivedSyn=new double[3];
                derivedNonsyn=new double[3];
                deleterious=new double[3];
                normalizedDerivedSyn=new double[3];
                normalizedDerivedNonSyn=new double[3];
                normalizedDerivedDel=new double[3];
                triadID=temp[0].get(0);
                for (int i = 0; i < cdsLen.length; i++) {
                    cdsLen[i]=Integer.parseInt(temp[i].get(1));
                    snpNum[i]=Integer.parseInt(temp[i].get(3))+Integer.parseInt(temp[i].get(5));
                    derivedSyn[i]=Double.parseDouble(temp[i].get(4));
                    derivedNonsyn[i]=Double.parseDouble(temp[i].get(6));
                    deleterious[i]=Double.parseDouble(temp[i].get(8));
                }
                for (int i = 0; i < 3; i++) {
                    normalizedDerivedSyn[i]=1000*10*derivedSyn[i]/(cdsLen[i]*snpNum[i]);
                    normalizedDerivedNonSyn[i]=1000*10*derivedNonsyn[i]/(cdsLen[i]*snpNum[i]);
                    normalizedDerivedDel[i]=1000*10*deleterious[i]/(cdsLen[i]*snpNum[i]);
                }
                sb=new StringBuilder();
                sb.append(triadID).append("\t");
                for (int i = 0; i < 3; i++) {
                    if (Double.isNaN(normalizedDerivedSyn[i])){
                        sb.append("NA").append("\t");
                    }else {
                        sb.append(NumberTool.format(normalizedDerivedSyn[i], 5)).append("\t");
                    }
                }
                if (Double.isNaN(normalizedDerivedSyn[0]) || Double.isNaN(normalizedDerivedSyn[1]) || Double.isNaN(normalizedDerivedSyn[2])){
                    region="NA";
                }else {
                    region=Standardization.getNearestPointIndex(normalizedDerivedSyn).getRegion();
                }
                sb.append(region).append("\t");
                for (int i = 0; i < 3; i++) {
                    if (Double.isNaN(normalizedDerivedNonSyn[i])){
                        sb.append("NA").append("\t");
                    }else {
                        sb.append(NumberTool.format(normalizedDerivedNonSyn[i], 5)).append("\t");
                    }
                }
                if (Double.isNaN(normalizedDerivedNonSyn[0]) || Double.isNaN(normalizedDerivedNonSyn[1]) || Double.isNaN(normalizedDerivedNonSyn[2])){
                    region="NA";
                }else {
                    region=Standardization.getNearestPointIndex(normalizedDerivedNonSyn).getRegion();
                }
                sb.append(region).append("\t");
                for (int i = 0; i < 3; i++) {
                    if (Double.isNaN(normalizedDerivedDel[i])){
                        sb.append("NA").append("\t");
                    }else {
                        sb.append(NumberTool.format(normalizedDerivedDel[i], 5)).append("\t");
                    }
                }
                if (Double.isNaN(normalizedDerivedDel[0]) || Double.isNaN(normalizedDerivedDel[1]) || Double.isNaN(normalizedDerivedDel[2])){
                    region="NA";
                }else {
                    region=Standardization.getNearestPointIndex(normalizedDerivedDel).getRegion();
                }
                sb.append(region);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }



}