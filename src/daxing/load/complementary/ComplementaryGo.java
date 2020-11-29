package daxing.load.complementary;

import com.google.common.io.Files;
import daxing.common.ArrayTool;
import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.individualIntrogression.P2;
import daxing.individualIntrogression.P3;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.GeneLoad;
import daxing.load.ancestralSite.IndividualChrLoad;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class ComplementaryGo {

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
        String[] subdir={"001_count","002_countMerge","003_hexaploidPseudohexaploid"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
//        IntStream.range(0, exonVCFFiles.size()).forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
//                taxonOutDirMap, e+1));
//        merge(new File(outDir, subdir[0]).getAbsolutePath(), new File(outDir, subdir[1]).getAbsolutePath());
        syntheticPseudohexaploidHexaploid(new File(outDir, subdir[1]).getAbsolutePath(), taxa_InfoDBFile,
                new File(outDir, subdir[2]).getAbsolutePath());
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
            tableTool.sortAsText("geneName");
            tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
            System.out.println("count merge finished: "+outName);
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
            for (int i = 0; i < keyList.size(); i++) {
                key=keyList.get(i);
                valuesFile=groupFileMap.get(key);
                if (key.equals("Landrace") || key.equals("Cultivar")){
                    p2=P2.valueOf(key);
                    for (int j = 0; j < valuesFile.size(); j++) {
                        Files.copy(valuesFile.get(j), new File(subdirFiles[1],"hexaploid."+p2.getAbbreviation()+"."+valuesFile.get(j).getName()));
                    }
                }else {
                    p3=P3.valueOf(PStringUtils.fastSplit(key,".").get(0));
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
                                                String pseudohexaploidOutDir, String outNamePrefix,
                                                int numPseudohexaploid){
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
        int[] randomIndex= ArrayTool.getRandomNonrepetitionArray(numPseudohexaploid, 0, combinations.size());
        IntStream.range(0, randomIndex.length).parallel().forEach(e->syntheticPseudohexaploid(diploidFiles.get(combinations.get(randomIndex[e])[1]),
                tetraploidFiles.get(combinations.get(randomIndex[e])[0]), pseudohexaploidOutDir, outNamePrefix));
    }

    private static void syntheticPseudohexaploid(List<File> tetraploidFiles, List<File> diploidFiles,
                                                String pseudohexaploidOutDir, String outNamePrefix){
        int numPseudohexaploid=tetraploidFiles.size()*diploidFiles.size();
        syntheticPseudohexaploid(tetraploidFiles, diploidFiles, pseudohexaploidOutDir, outNamePrefix, numPseudohexaploid);
        System.out.println("numPseudohexaploid: "+numPseudohexaploid+" ("+tetraploidFiles.size()+" "+diploidFiles.size()+")");
    }

    private static void syntheticPseudohexaploid(File diploidFile, File tetraploidFile, String outDir, String outNamePrefix){
        String diploidName=PStringUtils.fastSplit(diploidFile.getName(),".").get(0);
        String tetraploidName=PStringUtils.fastSplit(tetraploidFile.getName(), ".").get(0);
        File outFile=new File(outDir, outNamePrefix+"."+tetraploidName+"_"+diploidName+".triad.txt.gz");
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
            for (int i = 0; i < lines.size(); i++) {
                bw.write(lines.get(i));
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
