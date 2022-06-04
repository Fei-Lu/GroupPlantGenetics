package daxing.v2.TREE;

import com.google.common.base.Enums;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import daxing.common.chrrange.ChrRange;
import daxing.common.factors.WheatLineage;
import daxing.common.table.RowTableTool;
import daxing.common.utiles.IOTool;
import daxing.common.vmap2Group.GroupBySubcontinent;
import daxing.individualIntrogression.individual.Donor;
import daxing.individualIntrogression.pop2Indi.IndividualsFdDxyDonor;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang3.EnumUtils;
import org.apache.commons.lang3.StringUtils;
import pgl.infra.range.Range;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.*;
import java.text.NumberFormat;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Misc {

    public static void start(String rhoDir, String v2PosAlleleDir, int windowSize,
                             int stepSize,
                             String outDir){
        String[] subDirs={"001_v2_WindowStep","002_rhoWithPos"};
        for (int i = 0; i < subDirs.length; i++) {
            new File(outDir, subDirs[i]).mkdir();
        }
        Misc.extractSlidingWindow(v2PosAlleleDir, windowSize, stepSize, new File(outDir, subDirs[0]).getAbsolutePath());
        Misc.addPosToRho(rhoDir, new File(outDir, subDirs[0]).getAbsolutePath(),
                new File(outDir, subDirs[1]).getAbsolutePath());
//        Misc.outPutMissingRate(v2Dir, taxaInfoFile, new File(outDir, subDirs[2]).getAbsolutePath());
//        Misc.slidingWindowMissingRate(new File(outDir, subDirs[2]).getAbsolutePath(), windowSize, stepSize,
//                new File(outDir, subDirs[3]).getAbsolutePath());
//        Misc.transformToRefChr(new File(outDir, subDirs[1]).getAbsolutePath(), new File(outDir, subDirs[4]).getAbsolutePath());
    }

    /**
     * change individual introgression region to ./.
     * @param v2Dir
     * @param individualFdDxyDonorDir
     * @param outDir
     */
    public static void changeIntrogressionToMissing(String v2Dir, String individualFdDxyDonorDir, String outDir){
        List<File> files =IOTool.getFileListInDirEndsWith(v2Dir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",
                "_changeIntrogressionToMissing.vcf.gz")).toArray(String[]::new);
        IndividualsFdDxyDonor individualsFdDxyDonor = new IndividualsFdDxyDonor(individualFdDxyDonorDir);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            long start = System.nanoTime();
            System.out.println("start "+files.get(e).getName());
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw  =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, refChr;
                int refPos;
                List<String> temp, headerList;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                headerList = PStringUtils.fastSplit(line);
                bw.write(line);
                bw.newLine();
                StringBuilder sb = new StringBuilder();
                EnumSet<Donor> donorEnumSet;
                int cnt = 0;
                int count=0;
                while ((line=br.readLine())!=null){
                    cnt++;
                    temp =PStringUtils.fastSplit(line);
                    refChr = temp.get(0);
                    refPos = Integer.parseInt(temp.get(1));
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0,9))).append("\t");
                    int sumGenotype=0;
                    int genotype=0;
                    for (int i = 9; i < temp.size(); i++) {
                        donorEnumSet = individualsFdDxyDonor.getDonorFrom(headerList.get(i), refChr, refPos);
                        if (donorEnumSet.size()>0){
                            sb.append(".").append("\t");
                        }else if (temp.get(i).startsWith(".")){
                            sb.append(temp.get(i)).append("\t");
                        }else {
                            genotype = Integer.parseInt(temp.get(i));
                            sumGenotype+=genotype;
                            sb.append(temp.get(i)).append("\t");
                        }
                    }
                    if (sumGenotype == 0) continue;
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                    count++;
                    if (count% 1000000==0){
                        System.out.println("Written "+count+" variants");
                    }
                }
                System.out.println("Total "+count+"（"+cnt+")"+" variants had been written to "+new File(outDir,
                        outNames[e]).getName());
                System.out.println((cnt-count)+" variants had been removed due to introgression");
                System.out.println(files.get(e).getName()+ " completed in "+ Benchmark.getTimeSpanSeconds(start)+" " +
                        "seconds");
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void extractPopFromPhylip(String phylipDir, String taxaInfoFile, String outDir){
        List<File> fileList = IOTool.getFileListInDirEndsWith(phylipDir, ".phy");
        Map<String, String> taxaGroupbyContinentMap= RowTableTool.getMap(taxaInfoFile, 0, 36);
        IntStream.range(0,fileList.size()).forEach(e->{
            try (BufferedReader br = IOTool.getReader(fileList.get(e))) {
                String refChr = fileList.get(e).getName().substring(3,5);
                WheatLineage subgenome = WheatLineage.valueOf(refChr.substring(1));
                EnumSet<GroupBySubcontinent> enumSet = GroupBySubcontinent.getSubgenomeGroupByContinent(subgenome);
                EnumMap<GroupBySubcontinent, BufferedWriter> groupBWMap= new EnumMap<>(GroupBySubcontinent.class);
                BufferedWriter bw;
                for (GroupBySubcontinent group : enumSet){
                    bw = IOTool.getWriter(new File(outDir, "chr"+refChr+"_vmap2.1_"+group.name()+".phy"));
                    groupBWMap.put(group, bw);
                }
                String line;
                GroupBySubcontinent group;
                String[] temp;
                line=br.readLine();
                temp = StringUtils.split(line, " ");
                StringBuilder sb = new StringBuilder();
                Map<GroupBySubcontinent,Integer> groupCountMap=getTotalNumTaxa();
                for (Map.Entry<GroupBySubcontinent,BufferedWriter> entry : groupBWMap.entrySet()){
                    sb.setLength(0);
                    group = entry.getKey();
                    bw = entry.getValue();
                    sb.append(groupCountMap.get(group)).append(" ").append(temp[1]);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                while ((line=br.readLine())!=null){
                    temp = StringUtils.split(line.substring(0,20)," ");
                    if (! EnumUtils.isValidEnum(GroupBySubcontinent.class, taxaGroupbyContinentMap.get(temp[0]))) continue;
                    group=GroupBySubcontinent.valueOf(taxaGroupbyContinentMap.get(temp[0]));
                    bw = groupBWMap.get(group);
                    bw.write(line);
                    bw.newLine();
                }
                for (Map.Entry<GroupBySubcontinent,BufferedWriter> entry : groupBWMap.entrySet()){
                    bw = entry.getValue();
                    bw.flush();
                    bw.close();
                }
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });

    }

    private static Map<GroupBySubcontinent, Integer> getTotalNumTaxa(){
        Map<GroupBySubcontinent, Integer> map = new HashMap<>();
        map.put(GroupBySubcontinent.WE, 91);
        map.put(GroupBySubcontinent.DE,52);
        map.put(GroupBySubcontinent.FTT,50);
        map.put(GroupBySubcontinent.AT,36);
        map.put(GroupBySubcontinent.LR_AF,18);
        map.put(GroupBySubcontinent.LR_AM,26);
        map.put(GroupBySubcontinent.LR_CSA,19);
        map.put(GroupBySubcontinent.LR_WA,56);
        map.put(GroupBySubcontinent.LR_EU,58);
        map.put(GroupBySubcontinent.LR_EA,127);
        map.put(GroupBySubcontinent.CL,223);
        return map;
    }


    /**
     *
     * @param posAlleleDir v2 chr pos ref alt
     * @param windowSize snp number
     * @param stepSize snp num
     * @param outDir
     */
    public static void extractSlidingWindow(String posAlleleDir, int windowSize, int stepSize, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(posAlleleDir,".vcf.gz");
        String[] outFiles =
                files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".windowStep.txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outFiles[e]))) {
                String line;
                List<String> temp;
                List<ChrRange> chrRanges = new ArrayList<>();
                ChrRange chrRange;
                br.readLine();
                String chr = null;
                int pos;
                StringBuilder sb = new StringBuilder();
                Integer[] indexes_Inter = IntStream.range(0,stepSize).boxed().toArray(Integer[]::new);
                Arrays.sort(indexes_Inter,Collections.reverseOrder());
                int[] indexes = Arrays.stream(indexes_Inter).mapToInt(inter->inter).toArray();
                bw.write("Chr\tStart\tEnd");
                bw.newLine();
                while ((line=br.readLine()).startsWith("##")) continue;
                while ((line=br.readLine())!=null){
                    if (chrRanges.size()==windowSize){
                        sb.setLength(0);
                        sb.append(chrRanges.get(0).getChr()).append("\t");
                        sb.append(chrRanges.get(0).getStart()).append("\t");
                        sb.append(chrRanges.get(windowSize-1).getStart());
                        bw.write(sb.toString());
                        bw.newLine();
                        for (int i = 0; i < indexes.length; i++) {
                            chrRanges.remove(indexes[i]);
                        }
                    }
                    temp =PStringUtils.fastSplit(line.substring(0,30));
                    chr = temp.get(0);
                    pos = Integer.parseInt(temp.get(1));
                    chrRange = new ChrRange(chr, pos, pos+1);
                    chrRanges.add(chrRange);
                }
                // 6B多一行
                if (chr.equals("6B")){
                    int size = chrRanges.size();
                    sb.setLength(0);
                    sb.append(chrRanges.get(0).getChr()).append("\t");
                    sb.append(chrRanges.get(0).getStart()).append("\t");
                    sb.append(chrRanges.get(size-1).getStart());
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void addPosToRho(String rhoDir, String slidingWindowDir, String outDir){
        String[] refChrArray = RefV1Utils.getChromosomes();
        List<File> rhoFiles = IOTool.getFileListInDirEndsWith(rhoDir, "all.txt.gz");
        String[] outNames = rhoFiles.stream().map(File::getName).map(s -> s.replaceAll("_rhoPred_all.txt.gz","_RefChrPos.txt.gz")).toArray(String[]::new);
        List<File> files = IOTool.getFileListInDirEndsWith(slidingWindowDir, ".txt.gz");
        RowTableTool<String> rowTableToolPos, rowTableToolRho;
        String refChr;
        for (int i = 0; i < rhoFiles.size(); i++) {
            rowTableToolRho = new RowTableTool<>(rhoFiles.get(i).getAbsolutePath());
            refChr = rhoFiles.get(i).getName().substring(3,5);
            int refChrIndex = Arrays.binarySearch(refChrArray, refChr);
            rowTableToolPos = new RowTableTool<>(files.get(refChrIndex).getAbsolutePath());
            rowTableToolPos.addColumn("Rho", rowTableToolRho.getColumn(0));
            rowTableToolPos.write(new File(outDir, outNames[i]), IOFileFormat.TextGzip);
        }
    }

    public static void transformToRefChr(String rhoWithPosDir, String rhoWithRefChrPosDir){
        List<File> files =IOTool.getFileListInDirEndsWith(rhoWithPosDir, "pos.txt.gz");
        Map<String, List<File>> groupBySubcontinentListEnumMap =
                files.stream().collect(Collectors.groupingBy(file -> file.getName().substring(7).replaceAll("_rhoPred_all_pos.txt.gz","")));
        for (Map.Entry<String, List<File>> entry: groupBySubcontinentListEnumMap.entrySet()){
            List<File> fileList = entry.getValue();
            BufferedReader br1, br2;
            BufferedWriter bw;
            GroupBySubcontinent groupBySubcontinent = GroupBySubcontinent.valueOf(entry.getKey());
            String[] refChrArray = groupBySubcontinent.getGroupRefChr();
            try {
                for (int i = 0; i < fileList.size(); i=i+2) {
                    br1 = IOTool.getReader(fileList.get(i));
                    bw =IOTool.getWriter(new File(rhoWithRefChrPosDir, "chr"+refChrArray[i/2]+"_"+entry.getKey()+
                            "_RefChrPos.txt.gz"));
                    String line;
                    List<String> temp;
                    br1.readLine();
                    bw.write("Chr\tStart\tEnd\tRho");
                    bw.newLine();
                    short chrID;
                    int start, end, refStart, refEnd;
                    String refChr;
                    StringBuilder sb = new StringBuilder();
                    while ((line=br1.readLine())!=null){
                        temp =PStringUtils.fastSplit(line);
                        chrID = Short.parseShort(temp.get(0));
                        start = Integer.parseInt(temp.get(1));
                        end = Integer.parseInt(temp.get(2));
                        refChr = RefV1Utils.getChromosome(chrID, start);
                        refStart = RefV1Utils.getPosOnChromosome(chrID, start);
                        refEnd = RefV1Utils.getPosOnChromosome(chrID, end);
                        sb.setLength(0);
                        sb.append(refChr).append("\t").append(refStart).append("\t");
                        sb.append(refEnd).append("\t").append(temp.get(3));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    br1.close();
                    br2 = IOTool.getReader(fileList.get(i+1));
                    br2.readLine();
                    while ((line=br2.readLine())!=null){
                        temp =PStringUtils.fastSplit(line);
                        chrID = Short.parseShort(temp.get(0));
                        start = Integer.parseInt(temp.get(1));
                        end = Integer.parseInt(temp.get(2));
                        refChr = RefV1Utils.getChromosome(chrID, start);
                        refStart = RefV1Utils.getPosOnChromosome(chrID, start);
                        refEnd = RefV1Utils.getPosOnChromosome(chrID, end);
                        sb.setLength(0);
                        sb.append(refChr).append("\t").append(refStart).append("\t");
                        sb.append(refEnd).append("\t").append(temp.get(3));
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                    br2.close();
                    bw.flush();
                    bw.close();
                }
            } catch (Exception e) {
                e.printStackTrace();
            }
        }



    }

    public static void outPutMissingRate(String v2Dir, String taxaInfoFile, String outDir){
        List<File> files =IOTool.getFileListInDirEndsWith(v2Dir, ".vcf.gz");
        Multimap<String, String> groupTaxaMap= HashMultimap.create();
        Map<String, String> taxaGroupMap=RowTableTool.getMap(taxaInfoFile,0,36);
        for (Map.Entry<String, String> entry: taxaGroupMap.entrySet()){
            groupTaxaMap.put(entry.getValue(), entry.getKey());
        }
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz","_missingRate.txt.gz")).toArray(String[]::new);
        IntStream.range(0,files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, group;
                List<String> temp;
                String subgenome=RefV1Utils.getSubgenomeFromChrID(Integer.parseInt(files.get(e).getName().substring(3,6)));
                EnumSet<GroupBySubcontinent> enumSet = GroupBySubcontinent.getSubgenomeGroupByContinent(WheatLineage.valueOf(subgenome));
                EnumMap<GroupBySubcontinent,TIntArrayList> groupIndexesMap = new EnumMap<>(GroupBySubcontinent.class);
                for (GroupBySubcontinent groupBySubcontinent:enumSet){
                    groupIndexesMap.put(groupBySubcontinent, new TIntArrayList());
                }
                while ((line=br.readLine()).startsWith("##")) continue;
                temp =PStringUtils.fastSplit(line);
                for (int i = 9; i < temp.size(); i++) {
                    group = taxaGroupMap.get(temp.get(i));
                    if(!Enums.getIfPresent(GroupBySubcontinent.class, group).isPresent()) continue;
                    if (!enumSet.contains(GroupBySubcontinent.valueOf(group))) continue;
                    groupIndexesMap.get(GroupBySubcontinent.valueOf(group)).add(i);
                }
                StringBuilder sb = new StringBuilder();
                TIntArrayList indexList;
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setGroupingUsed(false);
                numberFormat.setMaximumFractionDigits(5);
                sb.append("ChrID\tPos\t");
                for (Map.Entry<GroupBySubcontinent,TIntArrayList> entry: groupIndexesMap.entrySet()){
                    sb.append(entry.getKey().name()).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp = PStringUtils.fastSplit(line);
                    sb.setLength(0);
                    sb.append(temp.get(0)).append("\t");
                    sb.append(temp.get(1)).append("\t");
                    for (Map.Entry<GroupBySubcontinent,TIntArrayList> entry : groupIndexesMap.entrySet()){
                        indexList = entry.getValue();
                        double missingCount=0;
                        for (int i = 0; i < indexList.size(); i++) {
                            if (temp.get(indexList.get(i)).equals("./.")){
                                missingCount++;
                            }
                        }
                        double missingRate = missingCount/indexList.size();
                        sb.append(numberFormat.format(missingRate)).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void slidingWindowMissingRate(String v2MissingRateDir, int windowSize, int stepSize, String outDir){
        List<File> files =IOTool.getFileListInDirEndsWith(v2MissingRateDir, ".txt.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp;
                line=br.readLine();
                temp = PStringUtils.fastSplit(line);
                TDoubleArrayList[] groupMissRateListArray = new TDoubleArrayList[temp.size()-2];
                List<Range> ranges = new ArrayList<>();
                StringBuilder sb = new StringBuilder();
                sb.append("ChrID\tStart\tEnd\tN\t");
                sb.append(String.join("\t", temp.subList(2,temp.size())));
                bw.write(sb.toString());
                bw.newLine();
                for (int i = 0; i < groupMissRateListArray.length; i++) {
                    groupMissRateListArray[i] = new TDoubleArrayList();
                }
                short chrID;
                int pos;
                NumberFormat numberFormat = NumberFormat.getInstance();
                numberFormat.setGroupingUsed(false);
                numberFormat.setMaximumFractionDigits(5);
                Integer[] indexes_Inter = IntStream.range(0,stepSize).boxed().toArray(Integer[]::new);
                Arrays.sort(indexes_Inter,Collections.reverseOrder());
                int[] indexes = Arrays.stream(indexes_Inter).mapToInt(inter->inter).toArray();
                while ((line=br.readLine())!=null){
                    if (ranges.size()==windowSize){
                        sb.setLength(0);
                        sb.append(ranges.get(0).chr).append("\t");
                        sb.append(ranges.get(0).start).append("\t");
                        sb.append(ranges.get(windowSize-1).start).append("\t");
                        sb.append(windowSize).append("\t");
                        for (int i = 0; i < groupMissRateListArray.length; i++) {
                            double missingRate_mean=(groupMissRateListArray[i].sum())/windowSize;
                            sb.append(numberFormat.format(missingRate_mean)).append("\t");
                        }
                        sb.deleteCharAt(sb.length()-1);
                        bw.write(sb.toString());
                        bw.newLine();
                        for (int index: indexes){
                            ranges.remove(index);
                            for (int i = 0; i < groupMissRateListArray.length; i++) {
                                groupMissRateListArray[i].removeAt(index);
                            }
                        }
                    }
                    temp = PStringUtils.fastSplit(line);
                    for (int i = 2; i < temp.size(); i++) {
                        groupMissRateListArray[i-2].add(Double.parseDouble(temp.get(i)));
                    }
                    chrID = Short.parseShort(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    ranges.add(new Range(chrID, pos, pos+1));
                }
                if (ranges.size() > 0){
                    sb.setLength(0);
                    sb.append(ranges.get(0).chr).append("\t");
                    sb.append(ranges.get(0).start).append("\t");
                    sb.append(ranges.get(ranges.size()-1).start).append("\t");
                    sb.append(ranges.size()).append("\t");
                    for (int i = 0; i < groupMissRateListArray.length; i++) {
                        double missingRate_mean=(groupMissRateListArray[i].sum())/ranges.size();
                        sb.append(numberFormat.format(missingRate_mean)).append("\t");
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

}
