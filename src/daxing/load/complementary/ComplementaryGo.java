package daxing.load.complementary;

import com.google.common.io.Files;
import daxing.common.*;
import daxing.individualIntrogression.P2;
import daxing.individualIntrogression.P3;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.GeneLoad;
import daxing.load.ancestralSite.IndividualChrLoad;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ComplementaryGo {

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/002_exonVCFAnno/002_exonAnnotationByDerivedSift/001_exonSNPAnnotation";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_exon/001_exonVCF";
        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String triadFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/triadGenes1.1.txt";
        String pgfFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/wheat_v1.1_Lulab.pgf";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftLoadComplementary";
        int blockGeneNum=20;
        go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, triadFile, blockGeneNum, pgfFile, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile, String triadFile,
                          int blockGeneNum, String pgfFile, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_hexaploidPseudohexaploid","004_triadsBlock",
                "005_mergeTriadsBlock","006_staticsValue"};
        for (String s : subdir) {
            new File(outDir, s).mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
        IntStream.range(0, exonVCFFiles.size()).forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxonOutDirMap, e+1, exonSNPAnnoDir));
        merge(new File(outDir, subdir[0]).getAbsolutePath(), new File(outDir, subdir[1]).getAbsolutePath());
//        syntheticPseudohexaploidHexaploid(new File(outDir, subdir[1]).getAbsolutePath(), taxa_InfoDBFile,
//                new File(outDir, subdir[2]).getAbsolutePath());
//        calculateLoadInfo(triadFile, pgfFile, blockGeneNum, new File(outDir, subdir[2]).getAbsoluteFile(), new File(outDir,
//                subdir[3]));
//        mergeTriadsBlockBySubspecies(new File(outDir, subdir[3]), new File(outDir, subdir[4]));
////        mergeTriadsBlockWithEightModel(new File(outDir, subdir[3]), new File(outDir, subdir[4]));
//        mergeAllIndividualTriadsBlock(new File(outDir, subdir[3]), pgfFile, new File(outDir, subdir[4]));
    }

    private static void go(File exonSNPAnnoFile, File exonVCFFile,
                           Map<String, File> taxonOutDirMap, int chr, String exonSNPAnnoDir){
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
            boolean isRefAlleleAncestral;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            byte[] indexGenotype; // syn nonsyn del
//            ExonSNPAnno exonSNPAnno=new ExonSNPAnno(exonSNPAnnoDir);
//            Table<Integer,Integer,Double> chrPosCorrRatioMap= exonSNPAnno.getCorrRatio();
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
                    if (genotype.equals("0/1")) continue;
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
        for (String taxonName : taxonNames) {
            taxonDir = new File(outDir, taxonName);
            taxonDir.mkdir();
            taxonOutDirMap.put(taxonName, taxonDir);
        }
        return taxonOutDirMap;
    }

    public static void merge(String inputDir, String outDir){
        List<File> dirs=IOUtils.getDirListInDir(inputDir);
        List<File> temp;
        RowTableTool<String> tableTool;
        String outName;
        for (File dir : dirs) {
            temp = IOUtils.getVisibleFileListInDir(dir.getAbsolutePath());
            tableTool = new RowTableTool<>(temp.get(0).getAbsolutePath());
            outName = PStringUtils.fastSplit(temp.get(0).getName(), ".").get(1);
            for (int j = 1; j < temp.size(); j++) {
                tableTool.add(new RowTableTool<>(temp.get(j).getAbsolutePath()));
            }
            tableTool.sortAsText("geneName");
            tableTool.write(new File(outDir, outName + ".txt.gz"), IOFileFormat.TextGzip);
            System.out.println("count merge finished: " + outName);
        }
    }

    public static void syntheticPseudohexaploidHexaploid(String inputDir, String vmapIIGroupFile, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        Map<String, String> taxonGroupMap=RowTableTool.getMap(vmapIIGroupFile,0,15);
        Map<String, List<File>> groupFileMap= files.stream().collect(Collectors.groupingBy(f->taxonGroupMap.get(PStringUtils.fastSplit(f.getName(), ".").get(0))));
        groupFileMap.remove("OtherHexaploid");
        groupFileMap.remove("OtherTetraploid");
        String[] subdirs={"001_pseudohexaploid","002_hexaploid"};
        File[] subdirFiles=new File[subdirs.length];
        for (int i = 0; i < subdirs.length; i++) {
            subdirFiles[i]=new File(outDir, subdirs[i]);
            subdirFiles[i].mkdir();
        }
        Set<String> keyset=groupFileMap.keySet();
        List<String> keyList=new ArrayList<>(keyset);
        List<File> valuesFile;
        List<File>[] p3Files=new List[P3.values().length];
        for (int i = 0; i < p3Files.length; i++) {
            p3Files[i]=new ArrayList<>();
        }
        String key;
        P3 p3;
        P2 p2;
        try {
            for (String s : keyList) {
                key = s;
                valuesFile = groupFileMap.get(key);
                if (key.equals("Landrace") || key.equals("Cultivar")) {
                    p2 = P2.valueOf(key);
                    for (File file : valuesFile) {
                        Files.copy(file, new File(subdirFiles[1], "hexaploid." + p2.getAbbreviation() + "." + file.getName()));
                    }
                } else {
                    p3 = P3.valueOf(PStringUtils.fastSplit(key, ".").get(0));
                    p3Files[p3.getIndex()].addAll(valuesFile);
                }
            }
            for (int i = 0; i < P3.values().length-1; i++) {
                p3=P3.newInstanceFrom(i);
                System.out.println("synthetic pseudohexaploid using "+p3.getAbbreviation()+" AT");
                syntheticPseudohexaploid(p3Files[i], p3Files[3], subdirFiles[0].getAbsolutePath(),
                        "pseudohexaploid."+p3.getAbbreviation()+"_AT");
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void syntheticPseudohexaploid(List<File> tetraploidFiles, List<File> diploidFiles,
                                                String pseudohexaploidOutDir, String outNamePrefix){
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
//        int[] randomIndex= ArrayTool.getRandomNonrepetitionArray(300, 0, combinations.size());
//        IntStream.range(0, randomIndex.length).parallel().forEach(e->syntheticPseudohexaploid(diploidFiles.get(combinations.get(randomIndex[e])[1]),
//                tetraploidFiles.get(combinations.get(randomIndex[e])[0]), pseudohexaploidOutDir, outNamePrefix));
        combinations.stream().parallel().forEach(e->syntheticPseudohexaploid(diploidFiles.get(e[1]),
                tetraploidFiles.get(e[0]), pseudohexaploidOutDir, outNamePrefix));
    }

//    private static void syntheticPseudohexaploid(List<File> tetraploidFiles, List<File> diploidFiles,
//                                                String pseudohexaploidOutDir, String outNamePrefix){
//        int numPseudohexaploid=tetraploidFiles.size()*diploidFiles.size();
//        syntheticPseudohexaploid(tetraploidFiles, diploidFiles, pseudohexaploidOutDir, outNamePrefix, numPseudohexaploid);
//        System.out.println("numPseudohexaploid: "+numPseudohexaploid+" ("+tetraploidFiles.size()+" "+diploidFiles.size()+")");
//    }

    private static void syntheticPseudohexaploid(File diploidFile, File tetraploidFile, String outDir, String outNamePrefix){
        String diploidName=PStringUtils.fastSplit(diploidFile.getName(),".").get(0);
        String tetraploidName=PStringUtils.fastSplit(tetraploidFile.getName(), ".").get(0);
        File outFile=new File(outDir, outNamePrefix+"."+tetraploidName+"_"+diploidName+".txt.gz");
        try (BufferedReader brDiploid = IOTool.getReader(diploidFile);
             BufferedReader brTetraploid = IOTool.getReader(tetraploidFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line, header;
            header=brDiploid.readLine();
            brTetraploid.readLine();
            bw.write(header);
            bw.newLine();
            List<String> lines=new ArrayList<>(80000);
            while ((line=brTetraploid.readLine())!=null){
                lines.add(line);
            }
            while ((line=brDiploid.readLine())!=null){
                lines.add(line);
            }
            Collections.sort(lines);
            for (String s : lines) {
                bw.write(s);
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void calculateLoadInfo(String triadsGeneFile, String pgfFile, int blockGeneNum,
                                         File individualLoadInfoFilesDir,
                                         File outDir){
        TriadsBlockUtils.writeTriadsBlock(triadsGeneFile, pgfFile, blockGeneNum,
                new File(outDir, "triadsBlock.txt.gz").getAbsolutePath());
        TriadsBlock[] triadsBlockArray=
                TriadsBlockUtils.readTriadBlock(new File(outDir, "triadsBlock.txt.gz").getAbsolutePath());
        List<File> dirs=IOUtils.getDirListInDir(individualLoadInfoFilesDir.getAbsolutePath());
        for (File dir : dirs) {
            List<File> individualLoadInfoFiles = IOUtils.getVisibleFileListInDir(dir.getAbsolutePath());
            String[] outFileNames = individualLoadInfoFiles.stream().map(File::getName)
                    .map(str -> str.replaceAll(".txt.gz", ".triadsBlock.txt.gz")).toArray(String[]::new);
            IntStream.range(0, individualLoadInfoFiles.size()).parallel().forEach(e -> calculateLoadInfo(triadsBlockArray,
                    individualLoadInfoFiles.get(e), new File(outDir,
                            outFileNames[e])));
        }
    }

    public static void calculateLoadInfo(TriadsBlock[] triadsBlockArray, File individualLoadInfoFile, File outFile){
        String taxonName=PStringUtils.fastSplit(individualLoadInfoFile.getName(),".").get(2);
        List<String> geneNameList=new ArrayList<>();
        List<int[]> loadInfoList=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(individualLoadInfoFile);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            br.readLine();
            String line;
            List<String> temp;
            int[] loadInfo;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                geneNameList.add(temp.get(0));
                loadInfo=new int[9];
                for (int i = 1; i < temp.size(); i++) {
                    loadInfo[i-1]=Integer.parseInt(temp.get(i));
                }
                loadInfoList.add(loadInfo);
            }
            List<String>[] blockGeneName;
            TIntArrayList[] abdIndexArray;
            int index;
            bw.write("TriadsBlock\tBlockGeneNum\tGenotypedGeneNum\t"+taxonName+"\tSlightlyModel\tStronglyModel");
            bw.newLine();
            IndividualTriadsBlockLoad individualTriadsBlockLoad;
            for (TriadsBlock triadsBlock : triadsBlockArray) {
                blockGeneName = triadsBlock.getBlockGeneName();
                abdIndexArray = new TIntArrayList[WheatLineage.values().length];
                for (int j = 0; j < abdIndexArray.length; j++) {
                    abdIndexArray[j] = new TIntArrayList();
                }
                for (int j = 0; j < blockGeneName.length; j++) {
                    for (int k = 0; k < blockGeneName[j].size(); k++) {
                        index = Collections.binarySearch(geneNameList, blockGeneName[j].get(k));
                        if (index < 0) continue;
                        abdIndexArray[j].add(index);
                    }
                }
                individualTriadsBlockLoad = new IndividualTriadsBlockLoad(triadsBlock, geneNameList,
                        loadInfoList, abdIndexArray);
                bw.write(individualTriadsBlockLoad.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * group by subspecies (according to file name)
     * @param triadsBlockInputDir triadsBlockInputDir
     * @param outDir AllIndividualTriadsBlockBySubspecies.txt.gz
     *               column: Triads, row: taxon
     */
    public static void mergeTriadsBlockBySubspecies(File triadsBlockInputDir, File outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadsBlockInputDir.getAbsolutePath());
        List<File> files1=files.stream().filter(file -> file.getName().startsWith("h") | file.getName().startsWith("p")).collect(Collectors.toList());
        Map<String, List<File>> typeFileListMap= files1.stream().collect(Collectors.groupingBy(file -> PStringUtils.fastSplit(file.getName(),".").get(1)));
        File triadsBlockFile= files.stream().filter(file -> file.getName().startsWith("triads")).collect(Collectors.toList()).get(0);
        String key, taxon;
        List<File> value;
        File file;
        File outFile=new File(outDir, "allIndividualTriadsBlockBySubspecies.txt.gz");
        StringBuilder headerSB=new StringBuilder();
        StringBuilder sb=new StringBuilder();
        List<String>[][] additiveDominanceSlightlyStrongly;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            List<String> triadIDList=RowTableTool.getColumnList(triadsBlockFile.getAbsolutePath(), 0);
            triadIDList.add(0,"Group");
            triadIDList.add(1,"AdditiveOrDominance");
            triadIDList.add(2,"SlightlyOrStrongly");
            headerSB.append("TaxonID").append("\t");
            headerSB.append(String.join("\t", triadIDList));
            bw.write(headerSB.toString());
            bw.newLine();
            String additiveOrDominance, slightlyOrStrongly;
            int count=0;
            System.out.println("start writing ...");
            for (Map.Entry<String,List<File>> entry : typeFileListMap.entrySet()){
                key=entry.getKey();
                value=entry.getValue();
                for (File item : value) {
                    file = item;
                    additiveDominanceSlightlyStrongly =
                            IndividualTriadsBlockLoad.getAdditiveDominanceSlightlyStronglyLoad(file);
                    taxon = PStringUtils.fastSplit(file.getName(), ".").get(2);
                    for (int j = 0; j < additiveDominanceSlightlyStrongly.length; j++) {
                        for (int k = 0; k < additiveDominanceSlightlyStrongly[j].length; k++) {
                            sb.setLength(0);
                            sb.append(taxon).append("\t").append(key).append("\t");
                            additiveOrDominance = j == 0 ? "additive" : "dominance";
                            slightlyOrStrongly = k == 0 ? "slightly" : "strongly";
                            sb.append(additiveOrDominance).append("\t").append(slightlyOrStrongly).append("\t");
                            sb.append(String.join("\t", additiveDominanceSlightlyStrongly[j][k]));
                            bw.write(sb.toString());
                            bw.newLine();
                        }
                    }
                    count++;
                    if (count % 200 == 0) {
                        System.out.println("writing " + count + " taxon");
                    }
                }
                System.out.println(key+" finished");
            }
            bw.flush();
            System.out.println("total "+count+" taxon has been write to "+outFile.getAbsolutePath());
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    public static void mergeTriadsBlockWithEightModel(File triadsBlockInputDir, File outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadsBlockInputDir.getAbsolutePath());
        List<File> files1=files.stream().filter(file -> file.getName().startsWith("h") | file.getName().startsWith("p")).collect(Collectors.toList());
        Map<String, List<File>> typeFileListMap= files1.stream().collect(Collectors.groupingBy(file -> PStringUtils.fastSplit(file.getName(),".").get(1)));
        File triadsBlockFile= files.stream().filter(file -> file.getName().startsWith("triads")).collect(Collectors.toList()).get(0);
        String key, taxon;
        List<File> value;
        File file;
        File outFile=new File(outDir, "mergedTriadsBlockWithEightModel.txt.gz");
        StringBuilder headerSB=new StringBuilder();
        StringBuilder sb=new StringBuilder();
        List<String>[] slightlyStronglyEightModel;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            List<String> triadIDList=RowTableTool.getColumnList(triadsBlockFile.getAbsolutePath(), 0);
            triadIDList.add(0,"Group");
            triadIDList.add(1,"SlightlyOrStrongly");
            headerSB.append("TaxonID").append("\t");
            headerSB.append(String.join("\t", triadIDList));
            bw.write(headerSB.toString());
            bw.newLine();
            String slightlyOrStrongly;
            int count=0;
            System.out.println("start writing ...");
            for (Map.Entry<String,List<File>> entry : typeFileListMap.entrySet()){
                key=entry.getKey();
                value=entry.getValue();
                for (File item : value) {
                    file = item;
                    taxon = PStringUtils.fastSplit(file.getName(), ".").get(2);
                    slightlyStronglyEightModel = IndividualTriadsBlockLoad.getSlightlyStronglyLoad(file);
                    for (int j = 0; j < slightlyStronglyEightModel.length; j++) {
                        slightlyOrStrongly = j == 0 ? "slightly" : "strongly";
                        sb.setLength(0);
                        sb.append(taxon).append("\t").append(key).append("\t").append(slightlyOrStrongly).append("\t");
                        sb.append(String.join("\t", slightlyStronglyEightModel[j]));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    count++;
                    if (count % 200 == 0) {
                        System.out.println("writing " + count + " taxon");
                    }
                }
                System.out.println(key+" finished");
            }
            bw.flush();
            System.out.println("total "+count+" taxon has been write to "+outFile.getAbsolutePath());
        } catch (IOException e) {
            e.printStackTrace();
        }

    }

    /**
     *
     * @param triadsBlockInputDir triadsBlockInputDir
     * @param pgfFile pgfFile
     * @param outDir allIndividualTriadsBlock.txt.gz
     *               row: triads column: taxon
     */
    public static void mergeAllIndividualTriadsBlock(File triadsBlockInputDir, String pgfFile, File outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(triadsBlockInputDir.getAbsolutePath());
        Predicate<File> pseudoHexaploid=file -> file.getName().startsWith("h") | file.getName().startsWith("p");
        List<File> files1=files.stream().filter(pseudoHexaploid).collect(Collectors.toList());
        File triadsBlockFile=files.stream().filter(pseudoHexaploid.negate()).collect(Collectors.toList()).get(0);
        Map<String, List<String>[]> triasBlockGeneNameListMap=getTriadsBlockGeneNameListMap(triadsBlockFile);
        File pseudohexaploidFile=new File(outDir, "pseudohexaploidInfo.txt.gz");
        File allIndividualTriadsBlockFile=new File(outDir, "allIndividualTriadsBlockBySubcontinent.txt.gz");
        List<String> triadsBlockIDList=RowTableTool.getColumnList(files1.get(0).getAbsolutePath(), 0, "\t");
        List<String> blockGeneNumList=RowTableTool.getColumnList(files1.get(0).getAbsolutePath(), 1, "\t");
        List<String> genotypedGeneNumList=RowTableTool.getColumnList(files1.get(0).getAbsolutePath(), 2, "\t");
        List<List<String>> columnList=new ArrayList<>(1500);
        TIntArrayList ifPseudohexaploidList=new TIntArrayList();
        List<String> taxonList=new ArrayList<>();
        List<String> taxonGroupList=new ArrayList<>();
        List<String> columnLoad;
        String taxon, taxonGroup;
        int ifPseudohexaploid, count=0;
        System.out.println(" start reading ...");
        for (File file : files1) {
            columnLoad = RowTableTool.getColumnList(file.getAbsolutePath(), 3, "\t");
            ifPseudohexaploid = PStringUtils.fastSplit(file.getName(), ".").get(0).equals("pseudohexaploid") ?
                    1 : 0;
            taxonGroup = PStringUtils.fastSplit(file.getName(), ".").get(1);
            taxon = PStringUtils.fastSplit(file.getName(), ".").get(2);
            ifPseudohexaploidList.add(ifPseudohexaploid);
            columnList.add(columnLoad);
            taxonGroupList.add(taxonGroup);
            taxonList.add(taxon);
            count++;
            if (count % 200 == 0) {
                System.out.println("reading " + count + " files");
            }
        }
        System.out.println("total "+ count+" files has been reading into memory");
        int cnt=0;
        PGF pgf=new PGF(pgfFile);
        pgf.sortGeneByName();
        try (BufferedWriter bw = IOTool.getWriter(allIndividualTriadsBlockFile);
             BufferedWriter bwInfo=IOTool.getWriter(pseudohexaploidFile)) {
            StringBuilder header=new StringBuilder();
            StringBuilder sb=new StringBuilder();
            header.append("TriadsBlockID\tBlockGeneNum\tGenotypedGeneNum\tRange\tCDSLen\t").append(String.join("\t",
                    taxonList));
            bw.write(header.toString());
            bw.newLine();
            String triadsBlockID;
            List<String>[] geneNameListArray;
            int geneIndex0, geneIndex, cdsLen0;
            ChrRange[] chrRangeABD;
            int[] cdsLenABD;
            Range range0, range;
            ChrRange chrRange0, chrRange;
            for (int i = 0; i < columnList.get(0).size(); i++) {
                sb.setLength(0);
                triadsBlockID=triadsBlockIDList.get(i);
                geneNameListArray=triasBlockGeneNameListMap.get(triadsBlockID);
                cdsLenABD=new int[3];
                chrRangeABD=new ChrRange[3];
                for (int j = 0; j < geneNameListArray.length; j++) {
                    geneIndex0=pgf.getGeneIndex(geneNameListArray[j].get(0));
                    range0=pgf.getGene(geneIndex0).getGeneRange();
                    chrRange0=ChrRange.changeToChrRange(range0);
                    cdsLen0=pgf.getGene(geneIndex0).getLongestTranscriptCDSLen();
                    for (int k = 1; k < geneNameListArray[j].size(); k++) {
                        geneIndex=pgf.getGeneIndex(geneNameListArray[j].get(k));
                        range=pgf.getGene(geneIndex).getGeneRange();
                        chrRange=ChrRange.changeToChrRange(range);
//                        if (range.getRangeStart() < range0.getRangeStart()){
//                            range0.setRangeStart(range.getRangeStart());
//                        }
//                        if (range.getRangeEnd() > range0.getRangeEnd()){
//                            range0.setRangeEnd(range.getRangeEnd());
//                        }
                        if (chrRange.getStart() < chrRange0.getStart()){
                            chrRange0.setStart(chrRange.getStart());
                        }
                        if (chrRange.getEnd() > chrRange0.getEnd()){
                            chrRange0.setEnd(chrRange.getEnd());
                        }
                        cdsLen0+=pgf.getGene(geneIndex).getLongestTranscriptCDSLen();
                    }
                    chrRangeABD[j]=chrRange0;
                    cdsLenABD[j]=cdsLen0;
                }
                sb.append(triadsBlockID).append("\t");
                sb.append(blockGeneNumList.get(i)).append("\t");
                sb.append(genotypedGeneNumList.get(i)).append("\t");
                for (ChrRange value : chrRangeABD) {
                    sb.append(value.getChr()).append(":").append(value.getStart()).append(",");
                    sb.append(value.getEnd() - 1).append(";");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (int k : cdsLenABD) {
                    sb.append(k).append(",");
                }
                sb.deleteCharAt(sb.length()-1).append("\t");
                for (List<String> strings : columnList) {
                    sb.append(strings.get(i)).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                cnt++;
            }
           bw.flush();
           header.setLength(0);
           header.append("TaxonID\tIfPseudohexaploid\tGroup");
           bwInfo.write(header.toString());
           bwInfo.newLine();
           String ifPseudo;
            for (int i = 0; i < ifPseudohexaploidList.size(); i++) {
                if (taxonGroupList.get(i).length() ==2) continue;
                sb.setLength(0);
                sb.append(taxonList.get(i)).append("\t");
                ifPseudo=ifPseudohexaploidList.get(i) == 1 ? "PseudoHexaploid" : "Hexaploid";
                sb.append(ifPseudo).append("\t");
                sb.append(taxonGroupList.get(i));
                bwInfo.write(sb.toString());
                bwInfo.newLine();
            }
            bwInfo.flush();
        } catch (IOException e) {
            System.out.println(cnt);
            e.printStackTrace();
        }
    }

    private static Map<String, List<String>[]> getTriadsBlockGeneNameListMap(File triadsBlockFile){
        Map<String, List<String>[]> triadsBlockGeneListMap = new HashMap<>();
        try (BufferedReader br = IOTool.getReader(triadsBlockFile)) {
            br.readLine();
            String line;
            List<String>[] value;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                value=new List[WheatLineage.values().length];
                for (int i = 0; i < value.length; i++) {
                    value[i]=PStringUtils.fastSplit(temp.get(4+i), ",");
                }
                triadsBlockGeneListMap.put(temp.get(0), value);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return triadsBlockGeneListMap;
    }
}
