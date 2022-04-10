package daxing.v2.relate;

import com.google.common.base.Enums;
import com.google.common.base.Predicate;
import com.google.common.collect.Table;
import daxing.common.chrrange.ChrRange;
import daxing.common.factors.WheatLineage;
import daxing.common.genotype.VCF;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.common.vmap2Group.GroupBySubcontinent;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import xiaohan.utils.CollectionsCalculator;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.Field;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Misc {

    /**
     * haploid_phased
     * @param ancestralVcfDir ancestral vcf Dir   chr1A 1B 1D ...
     * @param ancestralAlleleDir
     * @param outDir
     */
    public static void transformToHaploid(String ancestralVcfDir, String ancestralAlleleDir, String outDir){
        List<File> files= IOTool.getFileListInDirEndsWith(ancestralVcfDir, ".vcf.gz");
        List<File> ancestralFile = IOTool.getFileListInDirEndsWith(ancestralAlleleDir, ".txt.gz");
        String[] outNames = files.stream().map(File::getName).map(s->s.replaceAll(".vcf.gz","_haploid.vcf.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                Table<String, String, String> chrPosAncTable = RowTableTool.getTable(ancestralFile.get(e).getAbsolutePath()
                        , 2);
                String line;
                List<String> temp, tem;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                String chrID, pos, refAllele, altAllel, ancestralAllele;
                while ((line = br.readLine()) != null) {
                    temp = PStringUtils.fastSplit(line);
                    chrID = temp.get(0);
                    pos = temp.get(1);
                    refAllele = temp.get(3);
                    altAllel = temp.get(4);
                    if (!chrPosAncTable.contains(chrID, pos)) continue;
                    ancestralAllele = chrPosAncTable.get(chrID, pos);
                    if (!(ancestralAllele.equals(refAllele) || ancestralAllele.equals(altAllel))) continue;
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0,3))).append("\t");
                    if (ancestralAllele.equals(refAllele)){
                        sb.append(refAllele).append("\t").append(altAllel).append("\t");
                    }else if (ancestralAllele.equals(altAllel)){
                        sb.append(altAllel).append("\t").append(refAllele).append("\t");
                    }
                    sb.append(String.join("\t",temp.subList(5,9))).append("\t");
                    for (int i = 9; i < temp.size(); i++) {
                        tem = PStringUtils.fastSplit(temp.get(i),":");
                        if (refAllele.equals(ancestralAllele)){
                            if (tem.get(0).equals("0/0")){
                                sb.append("0").append("\t");
                            }else if (tem.get(0).equals("0/1")){
                                sb.append(".").append("\t");
                            }else if (tem.get(0).equals("1/1")){
                                sb.append("1").append("\t");
                            }else if (tem.get(0).equals("./.")){
                                sb.append(".").append("\t");
                            }
                        }else if (altAllel.equals(ancestralAllele)){
                            if (tem.get(0).equals("0/0")){
                                sb.append("1").append("\t");
                            }else if (tem.get(0).equals("0/1")){
                                sb.append(".").append("\t");
                            }else if (tem.get(0).equals("1/1")){
                                sb.append("0").append("\t");
                            }else if (tem.get(0).equals("./.")){
                                sb.append(".").append("\t");
                            }
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                    bw.flush();
                }
            }catch (Exception exception){
                exception.printStackTrace();
            }
        });
    }

    /**
     * .haps
     * Chromosome number [integer]
     * SNP ID [string]
     * SNP position [integer]
     * Ancestral allele [char]
     * Alternative allele [char]
     *
     * haps 文件要求 Chromosome number [integer]
     * 1A:0, 1B:1, 1D:2
     * @param v2_haploid_Dir
     * @param hapsOutFileDir
     */
    public static void prepareHaps(String v2_haploid_Dir, String hapsOutFileDir){
        List<File> files = IOTool.getFileListInDirEndsWith(v2_haploid_Dir, ".vcf.gz");
        String[] outNames_hapsOutFileNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".haps" +
                ".gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bwHaps = IOTool.getWriter(new File(hapsOutFileDir, outNames_hapsOutFileNames[e]))) {
                String line;
                List<String> temp;
                Integer chrNumber;
                StringBuilder sb = new StringBuilder();
                List<String> taxaList;
                Map<String, Integer> refChrToNumMap = getRefChrToNumMap();
                while ((line=br.readLine()).startsWith("##"))continue;
                temp = PStringUtils.fastSplit(line);
                taxaList = temp.subList(9, temp.size());
                while ((line=br.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    // haps 文件要求 Chromosome number [integer]
                    // 1A:0, 1B:1, 1D:2
                     chrNumber = refChrToNumMap.get(temp.get(0));
                     sb.append(chrNumber).append(" ").append(temp.get(2)).append(" ");
                     sb.append(temp.get(1)).append(" ");
                     sb.append(temp.get(3)).append(" ");
                     sb.append(temp.get(4)).append(" ");
                    for (int i = 9; i < temp.size(); i++) {
                        sb.append(temp.get(i)).append(" ");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bwHaps.write(sb.toString());
                    bwHaps.newLine();
                }
                bwHaps.flush();
//                bwExample.write("ID_1\tID_2\tmissing\n0\t0\t0");
//                bwExample.newLine();
//                for (int i = 0; i < taxaList.size(); i++) {
//                    sb.setLength(0);
//                    sb.append(taxaList.get(i)).append("\t").append("NA").append("\t").append("0");
//                    bwExample.write(sb.toString());
//                    bwExample.newLine();
//                }
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    private static Map<String, Integer> getRefChrToNumMap(){
        String[] refChr = RefV1Utils.getChromosomes();
        Map<String, Integer> refChrToNumMap = new HashMap<>();
        for (int i = 0; i < refChr.length; i++) {
            refChrToNumMap.put(refChr[i], i);
        }
        return refChrToNumMap;
    }

    /**
     *
     * @param haploidDir ancestral_haploid_imputation
     * @param outDir ancestral_pseudodiploid
     */
    public static void syntheticPseudodiploidFromHaploid(String haploidDir,  String taxaInfo, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(haploidDir,".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",
                "_pseudodiploid.vcf.gz")).toArray(String[]::new);
        IntStream.range(0,files.size()).parallel().forEach(e->{
            long start = System.nanoTime();
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                WheatLineage subgenome = WheatLineage.valueOf(files.get(e).getName().substring(4,5));
                bw.write(getPseudoDiploidHeader(taxaInfo, line, subgenome));
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                EnumMap<GroupBySubcontinent, TIntArrayList> groupByContinentIndexArrayMap = getGroupByContinentMap(taxaInfo, line);
                EnumSet<GroupBySubcontinent> groupBySubcontinentEnumSet;
                TIntArrayList indexList;
                while ((line=br.readLine())!=null){
                    int sumgeno =0;
                    temp =PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    sb.append(String.join("\t",temp.subList(0,9))).append("\t");
                    groupBySubcontinentEnumSet = GroupBySubcontinent.getSubgenomeGroupByContinent(subgenome);
                    for (Map.Entry<GroupBySubcontinent,TIntArrayList> entry : groupByContinentIndexArrayMap.entrySet()){
                        if (!groupBySubcontinentEnumSet.contains(entry.getKey())) continue;
                        indexList = entry.getValue();
                        for (int i = 0; i < indexList.size(); i=i+2) {
                            if (i==indexList.size()-1){
                                sb.append(temp.get(indexList.get(i))).append("|").append(temp.get(indexList.get(0))).append("\t");
                                sumgeno= sumgeno+Integer.parseInt(temp.get(indexList.get(i)))+Integer.parseInt(temp.get(indexList.get(0)));
                            }else {
                                sb.append(temp.get(indexList.get(i))).append("|").append(temp.get(indexList.get(i+1))).append("\t");
                                sumgeno= sumgeno+Integer.parseInt(temp.get(indexList.get(i)))+Integer.parseInt(temp.get(indexList.get(i+1)));
                            }
                        }
                    }
                    if (sumgeno == 0) continue;
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                System.out.println(files.get(e).getName()+ " completed in "+ Benchmark.getTimeSpanSeconds(start)+ " " +
                        "seconds");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    /**
     *
     * @param taxaInfo
     * @param vcfLine
     * @return
     */
    public static EnumMap<GroupBySubcontinent, TIntArrayList> getGroupByContinentMap(String taxaInfo, String vcfLine){
        Map<String, String> taxaGroupbyContinentMap=RowTableTool.getMap(taxaInfo, 0, 36);
        EnumMap<GroupBySubcontinent, TIntArrayList> groupByContinentIndexArrayMap = new EnumMap<>(GroupBySubcontinent.class);
        for (GroupBySubcontinent group : GroupBySubcontinent.values()){
            groupByContinentIndexArrayMap.put(group, new TIntArrayList());
        }
        String group;
        GroupBySubcontinent groupByContinent;
        List<String> headList = PStringUtils.fastSplit(vcfLine);
        for (int i = 9; i < headList.size(); i++) {
            group = taxaGroupbyContinentMap.get(headList.get(i));
            if (group.equals("FTT")){
                group = "FT";
            }
            if(Enums.getIfPresent(GroupBySubcontinent.class, group).isPresent()){
                groupByContinent = GroupBySubcontinent.valueOf(group);
                groupByContinentIndexArrayMap.get(groupByContinent).add(i);
            }
        }
        return groupByContinentIndexArrayMap;
    }

    private static String getPseudoDiploidHeader(String taxaInfo, String vcfLine, WheatLineage subgenome){
        List<String> temp = PStringUtils.fastSplit(vcfLine);
        EnumMap<GroupBySubcontinent, TIntArrayList> groupByContinentIndexArrayMap = getGroupByContinentMap(taxaInfo, vcfLine);
        StringBuilder sb = new StringBuilder();
        sb.append(String.join("\t", temp.subList(0,9))).append("\t");
        TIntArrayList indexList;
        EnumSet<GroupBySubcontinent> groupBySubcontinentEnumSet = GroupBySubcontinent.getSubgenomeGroupByContinent(subgenome);
        for (Map.Entry<GroupBySubcontinent,TIntArrayList> entry : groupByContinentIndexArrayMap.entrySet()){
            if (!groupBySubcontinentEnumSet.contains(entry.getKey())) continue;
            indexList = entry.getValue();
            for (int i = 0; i < indexList.size(); i=i+2) {
                if (i==indexList.size()-1){
                    sb.append(temp.get(indexList.get(i))).append("|").append(temp.get(indexList.get(0))).append("\t");
                }else {
                    sb.append(temp.get(indexList.get(i))).append("|").append(temp.get(indexList.get(i+1))).append("\t");
                }
            }
        }
        sb.deleteCharAt(sb.length()-1);
        return sb.toString();
    }

    public static void filterIntrogressionRegion(String pseudoDiploidDir, String introgressionRegionDir, String outDir){
        Predicate<File> ab = file -> file.getName().substring(4,5).equals("A") || file.getName().substring(4,5).equals("B");
        List<File> vcfFiles = IOTool.getFileListInDirEndsWith(pseudoDiploidDir, ".vcf.gz").stream().filter(ab).collect(Collectors.toList());
        List<File> introressionFiles = IOTool.getFileListInDirEndsWith(introgressionRegionDir, ".txt.gz").stream().filter(ab).collect(Collectors.toList());
        String[] outNames = vcfFiles.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", "_filterIntrogressionRegion.vcf.gz")).toArray(String[]::new);
        IntStream.range(0, vcfFiles.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(vcfFiles.get(e));
                 BufferedReader brRegion = IOTool.getReader(introressionFiles.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                List<ChrRange> chrRangeList_introgressionRegion = new ArrayList<>();
                String line, refChr;
                int start, end, pos;
                ChrRange chrRange;
                List<String> temp;
                brRegion.readLine();
                while ((line=brRegion.readLine())!=null){
                    temp =PStringUtils.fastSplit(line);
                    refChr = temp.get(0);
                    start = Integer.parseInt(temp.get(1));
                    end = Integer.parseInt(temp.get(2));
                    chrRange = new ChrRange(refChr, start, end);
                    chrRangeList_introgressionRegion.add(chrRange);
                }
                Collections.sort(chrRangeList_introgressionRegion);
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    refChr = temp.get(0);
                    pos = Integer.parseInt(temp.get(1));
                    chrRange = new ChrRange(refChr, pos, pos+1);
                    int rangeIndex = Collections.binarySearch(chrRangeList_introgressionRegion, chrRange);
                    if (rangeIndex >= 0) continue;
                    bw.write(line);
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void main(String[] args) {
//        String ancestralVcfDir = "/data4/home/aoyue/vmap2/daxing/analysis/008_vmap2_1062_spelt/001_subsetVmap2.1/004_vmap2_ancestral/001_ancestralVCF_byChrID";
//        String ancestralAlleleDir = "/data4/home/aoyue/vmap2/daxing/analysis/001_vmap2.1Before20200525/019_ancestral/004_ancestral/005_secer_hv/001_byChrID";
//        String prepareRelateInputFile_outDir = "/data4/home/aoyue/vmap2/daxing/analysis/008_vmap2_1062_spelt/020_relate/001_prepareInputFile";
//        String[] subDirs = {"001_vmap2_RefToAnc_Haploid_byChrID","002_vmap2_RefToAnc_Haploid_byRefChr",
//                "003_vmap2_RefToAnc_Haploid_byRefChr_imputation","004_haps"};
//        for (int i = 0; i < subDirs.length; i++) {
//            new File(prepareRelateInputFile_outDir, subDirs[i]).mkdir();
//        }
//        transformToHaploid(ancestralVcfDir, ancestralAlleleDir, new File(prepareRelateInputFile_outDir, subDirs[0]).getAbsolutePath());
//        VCF.fastMergeVCFtoChr(new File(prepareRelateInputFile_outDir, subDirs[0]).getAbsolutePath(),
//                new File(prepareRelateInputFile_outDir, subDirs[1]).getAbsolutePath());
//        // beagle imputation haploid
//        prepareHaps(new File(prepareRelateInputFile_outDir, subDirs[2]).getAbsolutePath(), new File(prepareRelateInputFile_outDir,
//                subDirs[3]).getAbsolutePath());
//
//        // ------------------
//        // pseudoDiploid
//
//        String haploidDir = args[0];
//        String taxaInfo = args[1];
//        String outDir = args[2];
//        Misc.syntheticPseudodiploidFromHaploid(haploidDir, taxaInfo, outDir);
    }
}
