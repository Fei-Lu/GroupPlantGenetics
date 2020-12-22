package daxing.load.ancestralSite;

import com.google.common.collect.Table;
import daxing.common.*;
import org.apache.commons.math3.distribution.BinomialDistribution;
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

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/002_exonVCFAnno/002_exonAnnotationByDerivedSift/001_exonSNPAnnotationByChrID";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/001_exonVCF";
        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String triadFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/triadGenes1.1_cdsLen_geneHC.txt";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftLoadComplementary";
        go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, triadFile, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile, String triadFile, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_retainTriad","004_standardization"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles=IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
        IntStream.range(0, exonVCFFiles.size()).forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxonOutDirMap, e+1));
        merge(new File(outDir, subdir[0]).getAbsolutePath(), new File(outDir, subdir[1]).getAbsolutePath());
//        retainTriad(triadFile, taxa_InfoDBFile, new File(outDir, subdir[1]).getAbsolutePath(),
//                new File(outDir, subdir[2]).getAbsolutePath());
//        normalizedTriadByAncestralNum(new File(outDir, subdir[2]).getAbsolutePath(), new File(outDir, subdir[3]).getAbsolutePath());
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
            int pos, depth;
            String geneName, genotype;
            List<String> genotypeList, depthList;
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
                if (!(isSyn || isNonsyn)) continue;
                isDeleterious=transcriptDB.isDeleterious(chr, pos);
                isRefAlleleAncestral=transcriptDB.isRefAlleleAncestral(chr, pos);
                geneName=transcriptDB.getGeneName(chr, pos);
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
                    if (isSyn){
                        indexGenotype[0]=0;
                    }else if (isDeleterious){
                        indexGenotype[0]=2;
                    } else {
                        indexGenotype[0]=1;
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

    public static void merge(String inputDir, String outDir){
        List<File> dirs=IOUtils.getDirListInDir(inputDir);
//        RowTableTool<String> table=new RowTableTool<>(vmapIIGroupFile);
//        Map<String,String> taxonGroupMap=table.getHashMap(0, 14);
//        String[] subDirs={"001_diploidCountMerge","002_tetraploidCountMerge","003_hexaploidCountMerge"};
//        File[] outFiles=new File[subDirs.length];
//        for (int i = 0; i < subDirs.length; i++) {
//            outFiles[i]=new File(outDir, subDirs[i]);
//            outFiles[i].mkdir();
//        }
//        Predicate<File> p_diploid= f->taxonGroupMap.get(f.getName()).equals("Ae.tauschii");
//        Predicate<File> p_tetraploid= f->taxonGroupMap.get(f.getName()).equals("Tetraploid");
//        Predicate<File> p_hexaploid= f->taxonGroupMap.get(f.getName()).equals("Hexaploid");
//        List<File> ldCl=dirs.stream().filter(p_diploid).collect(Collectors.toList());
        List<File> temp;
        RowTableTool<String> tableTool;
        String outName;
        for (int i = 0; i < dirs.size(); i++) {
            temp= IOUtils.getVisibleFileListInDir(dirs.get(i).getAbsolutePath());
            tableTool=new RowTableTool<>(temp.get(0).getAbsolutePath());
            outName= PStringUtils.fastSplit(temp.get(0).getName(), ".").get(1);
            for (int j = 1; j < temp.size(); j++) {
                tableTool.add(new RowTableTool<>(temp.get(j).getAbsolutePath()));
            }
            tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
            System.out.println("count merge finished: "+outName);
        }
    }

    public static void calculateExpectedSynNonDelPerTaxonPerSub(String countMergedDir,
                                                                String derivedProbabilityFile, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(countMergedDir);
        String[] taxonName= files.stream().map(File::getName)
                .map(str -> PStringUtils.fastSplit(str,".triad.txt").get(0)).toArray(String[]::new);
        String[] outFile= Arrays.stream(taxonName).map(s->s+".expected.txt").toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->calculateExpectedSynNonDelPerTaxonPerSub(files.get(e),
                taxonName[e], derivedProbabilityFile,new File(outDir, outFile[e])));
    }

    public static double[] getTaxonSynNonDelDerivedProbability(String taxon, String subgenome,
                                                                  String derivedProbabilityFile){
        RowTableTool<String> table=new RowTableTool<>(derivedProbabilityFile);
        Table<String,String,String> synTable=table.getTable(0, 1,7);
        Table<String,String,String> nonTable=table.getTable(0, 1,8);
        Table<String,String,String> delTable=table.getTable(0, 1,9);
        double[] synNonDel=new double[3];
        Arrays.fill(synNonDel, -1);
        synNonDel[0]=Float.parseFloat(synTable.get(taxon,subgenome));
        synNonDel[1]=Float.parseFloat(nonTable.get(taxon,subgenome));
        synNonDel[2]=Float.parseFloat(delTable.get(taxon,subgenome));
        return synNonDel;
    }

    public static void calculateExpectedSynNonDelPerTaxonPerSub(File mergedTaxonFile, String taxonName,
                                                                String derivedProbabilityFile,
                                                                File outFile){
        try (BufferedReader br = IOTool.getReader(mergedTaxonFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String header=br.readLine();
            bw.write(header+"\tAncestralNum\tExpectedSyn\tExpectedNonsyn\tExpectedDel\tp_lowerTail_syn" +
                    "\tp_upperTaill_syn\tp_lowerTail_non\tp_upperTaill_non\tp_lowerTail_del\tp_upperTaill_del");
            bw.newLine();
            StringBuilder sb=new StringBuilder();
            String line;
            List<String> temp;
            int synNum, nonNum, ancestralNum;
            int[] synNonDelDerivedNum=new int[3];
            Arrays.fill(synNonDelDerivedNum, -1);
            double expectedDerivedNum;
            String sub;
            BinomialDistribution binomialDistribution;
            double lowerTail, upperTail;
            List<String> strList=Collections.nCopies(9, "NA");
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                sub=temp.get(2).substring(8,9);
                synNum=Integer.parseInt(temp.get(3));
                synNonDelDerivedNum[0]=Integer.parseInt(temp.get(4));
                nonNum=Integer.parseInt(temp.get(5));
                synNonDelDerivedNum[1]=Integer.parseInt(temp.get(6));
                synNonDelDerivedNum[2]=Integer.parseInt(temp.get(8));
                ancestralNum=synNum+nonNum;
                sb.setLength(0);
                if (ancestralNum==0){
                    sb.append(String.join("\t", temp)).append("\t");
                    sb.append(ancestralNum).append("\t");
                    sb.append(String.join("\t",strList));
                    bw.write(sb.toString());
                    bw.newLine();
                    continue;
                }
                sb.append(String.join("\t", temp)).append("\t");
                sb.append(ancestralNum).append("\t");
                double[] synNonDelP=getTaxonSynNonDelDerivedProbability(taxonName, sub, derivedProbabilityFile);
                for (int i = 0; i < synNonDelP.length; i++) {
                    expectedDerivedNum=NumberTool.format(ancestralNum*synNonDelP[i], 5);
                    sb.append(expectedDerivedNum).append("\t");
                }
                for (int i = 0; i < synNonDelP.length; i++) {
                    binomialDistribution=new BinomialDistribution(ancestralNum,synNonDelP[i]);
                    lowerTail =binomialDistribution.cumulativeProbability(synNonDelDerivedNum[i]);
                    upperTail = (1-binomialDistribution.cumulativeProbability(synNonDelDerivedNum[i]))+binomialDistribution.probability(synNonDelDerivedNum[i]);
                    sb.append(NumberTool.format(lowerTail,5)).append("\t").append(NumberTool.format(upperTail, 5)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void normalizedTriadByAncestralNumAndAddExpectedRegion(String inputDir, String outDir){
        List<File> files=IOTool.getVisibleDir(inputDir);
        IntStream.range(0, files.size()).forEach(f->normalizedTriadByAncestralNumAndAddExpectedRegion(files.get(f), outDir));
    }

    public static void normalizedTriadByAncestralNumAndAddExpectedRegion(File inputFile, String outDir){
        String fileName=inputFile.getName().replaceAll("triad", "triad.normalized");
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(new File(outDir, fileName))) {
            bw.write("TriadID\tnormalizedNumDerivedInSynA" +
                    "\tnormalizedNumDerivedInSynB\tnormalizedNumDerivedInSynD\tsynExpectedRegion" +
                    "\tnormalizedNumDerivedInNonsynA" +
                    "\tnormalizedNumDerivedInNonsynB\tnormalizedNumDerivedInNonsynD\tnonsynExpectedRegion" +
                    "\tnormalizedNumDerivedInHGDeleteriousA\tnormalizedNumDerivedInHGDeleteriousB" +
                    "\tnormalizedNumDerivedInHGDeleteriousD\tdelExpectedRegion");
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
            double[] p_valueSyn, p_valueNon, p_valueDel;
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
                p_valueSyn=new double[6];
                p_valueNon=new double[6];
                p_valueDel=new double[6];
                Arrays.fill(p_valueSyn, -1);
                Arrays.fill(p_valueNon, -1);
                Arrays.fill(p_valueDel, -1);
                triadID=temp[0].get(0);
                for (int i = 0; i < cdsLen.length; i++) {
                    cdsLen[i]=Integer.parseInt(temp[i].get(1));
                    snpNum[i]=Integer.parseInt(temp[i].get(3))+Integer.parseInt(temp[i].get(5));
                    derivedSyn[i]=Double.parseDouble(temp[i].get(4));
                    derivedNonsyn[i]=Double.parseDouble(temp[i].get(6));
                    deleterious[i]=Double.parseDouble(temp[i].get(8));
                    if (temp[i].get(14).equals("NA"))continue;
                    p_valueSyn[2*i]=Double.parseDouble(temp[i].get(14));
                    p_valueSyn[2*i+1]=Double.parseDouble(temp[i].get(15));
                    p_valueNon[2*i]=Double.parseDouble(temp[i].get(16));
                    p_valueNon[2*i+1]=Double.parseDouble(temp[i].get(17));
                    p_valueDel[2*i]=Double.parseDouble(temp[i].get(18));
                    p_valueDel[2*i+1]=Double.parseDouble(temp[i].get(19));
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
                    region= ExpectedRegion.getExpectedRegion(p_valueSyn);
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
                    region=ExpectedRegion.getExpectedRegion(p_valueNon);
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
                    region=ExpectedRegion.getExpectedRegion(p_valueDel);
                }
                sb.append(region);
                bw.write(sb.toString());
                bw.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void retainTriad(String triadFile, String taxa_InfoDBFile, String inputDir, String outDir){
        List<File> files=IOTool.getVisibleDir(inputDir);
        IntStream.range(0, files.size()).parallel().forEach(e->retainTriad(triadFile, taxa_InfoDBFile, files.get(e), outDir));
    }

    private static void retainTriad(String triadFile, String taxa_InfoDBFile, File inputFile, String outDir){
        Triads triads =new Triads(triadFile);
        Map<String,String> taxaSubCharMap=RowTableTool.getMap(taxa_InfoDBFile, 0, 3);
        String taxonName=PStringUtils.fastSplit(inputFile.getName(),".").get(0);
        Ploidy ploidy=Ploidy.newInstanceFromSubChar(taxaSubCharMap.get(taxonName));
        RowTableTool<String> table=new RowTableTool<>(inputFile.getAbsolutePath());
        List<String> geneNames=table.getColumn(0);
        Collections.sort(geneNames);
        String[] threeName;
        int[] index=new int[ploidy.getSubgenomewNum()];
        String triadID, taxon, temp;
        taxon=PStringUtils.fastSplit(inputFile.getName(), ".").get(0);
        StringBuilder sb;
        BufferedWriter bw;
        int[] cdsLen;
        try {
            bw=IOTool.getWriter(new File(outDir, taxon + ".triad.txt.gz"));
            bw.write("TriadID\tcdsLen\tGeneName\tnumSyn\tnumDerivedInSyn\tnumHeterInSyn\tnumNonsyn" +
                    "\tnumDerivedInNonsyn\tnumHeterInNonsyn\tnumHGDeleterious\tnumDerivedInHGDeleterious\tnumHeterInHGDeleterious\tsubgenome");
            bw.newLine();
            String subgenome, geneName;
            boolean lessThanZero=false;
            for (int i = 0; i < triads.getTriadRecordNum(); i++) {
                threeName= triads.getTriadRecord(i).getTriadGeneNameArray();
                if (ploidy.getSubgenomewNum()>1){
                    for (int j = 0; j < index.length; j++) {
                        index[j]= Collections.binarySearch(geneNames, threeName[j]);
                    }
                }else {
                    index[0]=Collections.binarySearch(geneNames, threeName[2]);
                }
                for (int j = 0; j < index.length; j++) {
                    if (index[j] < 0){
                        lessThanZero=true;
                        break;
                    }else {
                        lessThanZero=false;
                    }
                }
                if (lessThanZero) continue;
                triadID= triads.getTraidID(i);
                cdsLen= triads.getCDSLen(i);
                for (int j = 0; j < index.length; j++) {
                    temp=String.join("\t", table.getRow(index[j]));
                    sb=new StringBuilder();
                    sb.append(triadID).append("\t").append(cdsLen[j]).append("\t");
                    sb.append(temp).append("\t");
                    geneName=table.getRow(index[j]).get(0);
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
        System.out.println("retainTriad finished: "+taxonName);
    }

    public static void normalizedTriadByAncestralNum(String inputDir,String outDir){
        List<File> files=IOTool.getVisibleDir(inputDir);
        IntStream.range(0, files.size()).forEach(f->normalizedTriadByAncestralNum(files.get(f), outDir));
    }

    public static void normalizedTriadByAncestralNum(File inputFile, String outDir){
        String fileName=inputFile.getName().replaceAll("triad", "triad.normalized");
        try (BufferedReader br = IOTool.getReader(inputFile);
             BufferedWriter bw =IOTool.getWriter(new File(outDir, fileName))) {
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
                    region= Standardization.getNearestPointIndex(normalizedDerivedSyn).getRegion();
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

    public static void retainTriadPseudohexaploid(String triadFile, String triadIDFile, List<File> inputFiles,
                                                  String outDir){
        RowTableTool<String> table=new RowTableTool<>(triadIDFile);
        List<String> triadIDList=table.getColumn(0);
        inputFiles.stream().parallel().forEach(f->retainTriadPseudohexaploid(triadFile, triadIDList, f, outDir));
    }

    private static void retainTriadPseudohexaploid(String triadFile, List<String> triadIDList, File inputFile, String outDir){
        Triads triads =new Triads(triadFile);
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
            bw=IOTool.getWriter(new File(outDir, taxon + ".triad.txt.gz"));
            bw.write("TriadID\tcdsLen\tGeneName\tnumSyn\tnumDerivedInSyn\tnumNonsyn\tnumDerivedInNonsyn" +
                    "\tnumHGDeleterious\tnumDerivedInHGDeleterious\tsubgenome");
            bw.newLine();
            String subgenome, geneName;
            for (int i = 0; i < triads.getTriadRecordNum(); i++) {
                if (!triadIDList.contains(triads.getTraidID(i))) continue;
                threeName= triads.getTriadRecord(i).getTriadGeneNameArray();
                for (int j = 0; j < indexABD.length-1; j++) {
                    indexABD[j]= Collections.binarySearch(geneNames, threeName[j]);
                }
                if (indexABD[0] < 0) continue;
                if (indexABD[1] < 0) continue;
//                if (indexABD[2] < 0) continue;
                triadID= triads.getTraidID(i);
                cdsLen= triads.getCDSLen(i);
                for (int j = 0; j < indexABD.length-1; j++) {
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

    public static void syntheticPseudohexaploid(String tetraploidDir, String diploidDir, String outDir,
                                                int numPseudohexaploid){
        List<File> diploidFiles=IOUtils.getVisibleFileListInDir(diploidDir);
        List<File> tetraploidFiles=IOUtils.getVisibleFileListInDir(tetraploidDir);
        List<int[]> combinations=new ArrayList<>();
        int[] combination;
        for (int i = 0; i < diploidFiles.size(); i++) {
            for (int j = 0; j < tetraploidFiles.size(); j++) {
                combination=new int[2];
                combination[0]=j;
                combination[1]=i;
                combinations.add(combination);
            }
        }
        int[] randomIndex=ArrayTool.getRandomNonrepetitionArray(numPseudohexaploid, 0, combinations.size());
        IntStream.range(0, randomIndex.length).parallel().forEach(e->syntheticPseudohexaploid(diploidFiles.get(combinations.get(randomIndex[e])[1]),
                tetraploidFiles.get(combinations.get(randomIndex[e])[0]), outDir));
    }

    private static void syntheticPseudohexaploid(File diploidFile, File tetraploidFile, String outDir){
        String diploidName=PStringUtils.fastSplit(diploidFile.getName(),".").get(0);
        String tetraploidName=PStringUtils.fastSplit(tetraploidFile.getName(), ".").get(0);
        File outFile=new File(outDir, tetraploidName+"."+diploidName+".triad.txt.gz");
        try (BufferedReader brDiploid = IOTool.getReader(diploidFile);
             BufferedReader brTetraploid = IOTool.getReader(tetraploidFile);
             BufferedWriter bw =IOTool.getWriter( outFile)) {
            String line, header;
            List<String> temp;
            header=brDiploid.readLine();
            brTetraploid.readLine();
            bw.write(header);
            bw.newLine();
            while ((line=brTetraploid.readLine())!=null){
                bw.write(line);
                bw.newLine();
                line=brTetraploid.readLine();
                bw.write(line);
                bw.newLine();
                line=brDiploid.readLine();
                bw.write(line);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}