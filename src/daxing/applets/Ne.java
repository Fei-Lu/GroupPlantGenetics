package daxing.applets;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Multimap;
import com.google.common.collect.Table;
import com.google.common.collect.TreeMultimap;
import daxing.common.*;
import daxing.common.vmap2Group.Group;
import daxing.common.vmap2Group.GroupType;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.util.CombinatoricsUtils;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Ne {

//    public static void main(String[] args) {
////        syntheticPseudoDiploid();
////        transformRefToAncestral();
////        bulidSMC();
////        mergeSMC();
//    }

    /**
     *
     * @param vmap2InputDir
     * @param taxaInfo_depthFile
     * @param outDir
     * @param group
     * @param subgenomeCombination
     * @param pseudoDiploidNum 是taxaInfo_depthFile文件中，个体数最少的亚群的组合数
     */
    public static void syntheticPseudoDiploid(String vmap2InputDir, String taxaInfo_depthFile, String outDir,
                                              Group group, SubgenomeCombination subgenomeCombination,
                                              int pseudoDiploidNum){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> files = IOUtils.getFileListInDirEndsWith(vmap2InputDir, "gz");
        List<File> subgenomeFiles=getFiles(files, subgenomeCombination);
        String[] subgenomeOutFileName= subgenomeFiles.stream().map(File::getName).
                map(s -> s.replaceAll(".vcf.gz", ".PseudoDiploid.vcf.gz")).toArray(String[]::new);
        Map<String,String> taxaGroupSubgenomeMap=getTaxaGroupSubgenomeMap(taxaInfo_depthFile, group, subgenomeCombination);
        IntStream.range(0, subgenomeFiles.size()).parallel().forEach(e->syntheticPseudoDiploid(group, pseudoDiploidNum,
                subgenomeFiles.get(e),taxaGroupSubgenomeMap, new File(outDir, subgenomeOutFileName[e])));
        System.out.println(subgenomeCombination.name()+" subgenome completed");
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void syntheticPseudoDiploid(Group group, int pseudoDiploidNumInPop, File vmap2File,
                                               Map<String, String> taxaGroupMap, File outFile){
        Set<String> groupSet=getValueSet(taxaGroupMap);
        List<String> groupList=new ArrayList<>(groupSet);
        Collections.sort(groupList);
        try (BufferedReader br = IOTool.getReader(vmap2File);
             BufferedWriter bw =IOTool.getWriter(outFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            TIntArrayList[] distinguishedLineageIndexArray= getDistinguishedLineageIndex(taxaGroupMap, line);
            Iterator<int[]>[] iteratorCombinations=new Iterator[distinguishedLineageIndexArray.length];
            for (int i = 0; i < distinguishedLineageIndexArray.length; i++) {
                iteratorCombinations[i]= CombinatoricsUtils.combinationsIterator(distinguishedLineageIndexArray[i].size(), 2);
            }
            temp=PStringUtils.fastSplit(line);
            StringBuilder sb=new StringBuilder();
            sb.append(String.join("\t", temp.subList(0,9))).append("\t");
            List<int[]>[] indexListArray=new List[distinguishedLineageIndexArray.length];
            for (int i = 0; i < indexListArray.length; i++) {
                indexListArray[i]=new ArrayList<>();
            }
            int[] indexArray;
            GroupType groupType;
            List<int[]> combinationTotalList;
            for (int i = 0; i < iteratorCombinations.length; i++) {
                groupType = GroupType.newInstanceFrom(groupList.get(i), group);
                Iterator<int[]> iterator=iteratorCombinations[i];
                combinationTotalList=new ArrayList<>();
                while (iterator.hasNext()){
                    int[] combination= iterator.next();
                    indexArray=new int[2];
                    indexArray[0]=distinguishedLineageIndexArray[i].get(combination[0]);
                    indexArray[1]=distinguishedLineageIndexArray[i].get(combination[1]);
                    combinationTotalList.add(indexArray);
                }
                Collections.shuffle(combinationTotalList);
                for (int j = 0; j < pseudoDiploidNumInPop; j++) {
                    indexArray=combinationTotalList.get(j);
                    sb.append(groupType.getGroupAbbrev()).append("_").append(temp.get(indexArray[0])).append("_").append(temp.get(indexArray[1])).append("\t");
                    indexListArray[i].add(indexArray);
                }
            }
            sb.deleteCharAt(sb.length()-1);
            bw.write(sb.toString());
            bw.newLine();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                sb.setLength(0);
                sb.append(String.join("\t", temp.subList(0,9))).append("\t");
                for (List<int[]> indexList : indexListArray){
                    for (int[] combination: indexList){
                        if (temp.get(combination[0]).startsWith("./.") || temp.get(combination[1]).startsWith("0/1")){
                            sb.append(".|.").append("\t");
                        }else {
                            sb.append(temp.get(combination[0]).substring(0,1)).append("|");
                            sb.append(temp.get(combination[1]).substring(0,1)).append(":99,99:99,99,99").append("\t");
                        }
                    }
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
            System.out.println(outFile.getName()+" had completed");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Map<String, String> getTaxaGroupSubgenomeMap(String taxaInfo_depthFile, Group group,
                                                         SubgenomeCombination subgenomeCombination){
        List<String> subgenomeGroupList=group.getGroup(subgenomeCombination);
        Map<String, String> taxaGroupSubgenomeMap=new HashMap<>();
        Map<String, String> taxaGroupMap=null;
        switch (group){
            case Subcontinent:
                taxaGroupMap=RowTableTool.getMap(taxaInfo_depthFile,0,2);
                break;
            case Subspecies:
                taxaGroupMap=RowTableTool.getMap(taxaInfo_depthFile,0,3);
                taxaGroupMap.remove("NA");
                break;
        }
        for (Map.Entry<String, String> entry : taxaGroupMap.entrySet()){
            if (!subgenomeGroupList.contains(entry.getValue())) continue;
            taxaGroupSubgenomeMap.put(entry.getKey(), entry.getValue());
        }
        return taxaGroupSubgenomeMap;
    }

    /**
     *
     * @param taxaGroupMap taxaGroupMap
     * @param vcfHeader vcfHeader
     * @return distinguishedLineageHeaderList and distinguishedLineageIndexList
     */
    private static TIntArrayList[] getDistinguishedLineageIndex(Map<String, String> taxaGroupMap, String vcfHeader){
        Set<String> groupSet=getValueSet(taxaGroupMap);
        List<String> groupList=new ArrayList<>(groupSet);
        Collections.sort(groupList);
        List<String> headerList=PStringUtils.fastSplit(vcfHeader);
//        List<String> distinguishedLineageHeaderList=new ArrayList<>();
        TIntArrayList[] distinguishedLineageIndexArray=new TIntArrayList[groupList.size()];
        for (int i = 0; i < distinguishedLineageIndexArray.length; i++) {
            distinguishedLineageIndexArray[i]=new TIntArrayList();
        }
        for (int i = 9; i < headerList.size(); i++) {
            if (!taxaGroupMap.containsKey(headerList.get(i))) continue;
            int groupIndex=Collections.binarySearch(groupList, taxaGroupMap.get(headerList.get(i)));
            distinguishedLineageIndexArray[groupIndex].add(i);
//            distinguishedLineageHeaderList.add(headerList.get(i));
        }
        return distinguishedLineageIndexArray;
    }

    private static Set<String> getValueSet(Map<String, String> taxaGroupMap){
        return taxaGroupMap.entrySet().stream().map(entry->entry.getValue()).collect(Collectors.toSet());
    }

    private static List<File> getFiles(List<File> files, SubgenomeCombination subgenomeCombination){
        int[] ab=WheatLineage.ablineage();
        int[] d=WheatLineage.D.getChrID();
        TIntArrayList abList=new TIntArrayList(ab);
        TIntArrayList dList = new TIntArrayList(d);
        switch (subgenomeCombination){
            case AB:
                return files.stream().filter(file -> abList.contains(Integer.parseInt(file.getName().substring(3,6)))).collect(Collectors.toList());
            case D:
                return files.stream().filter(file -> dList.contains(Integer.parseInt(file.getName().substring(3,6)))).collect(Collectors.toList());

        }
        return null;
    }

    private static Table<Integer, Integer, String> getAncestral(String ancestralFile){
        Table<Integer, Integer, String> table= HashBasedTable.create();
        try (BufferedReader br = IOTool.getReader(ancestralFile)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                table.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)), temp.get(2));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return table;
    }

    public static void transformRefToAncestral(String vcfDir, String ancestralDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> vcfFiles=IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> ancestralFiles=IOUtils.getVisibleFileListInDir(ancestralDir);
        String[] outNames=vcfFiles.stream().map(File::getName).map(s->s.replaceAll("\\.vcf.gz",".RefToAncestral" +
                ".vcf")).toArray(String[]::new);
        IntStream.range(0, ancestralFiles.size()).parallel().forEach(e->transformRefToAncestral(vcfFiles.get(e),
                ancestralFiles.get(e), new File(outDir, outNames[e])));
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void transformRefToAncestral(File vcfFile, File ancestralFile, File outFile){
        Table<Integer, Integer, String> ancestralTable=getAncestral(ancestralFile.getAbsolutePath());
        try (BufferedReader br = IOTool.getReader(vcfFile);
             BufferedWriter bw=IOTool.getWriter(outFile)) {
            String line, refAllele, altAllele, ancestralAllele;
            List<String> temp;
            int chr, pos;
            int total=0;
            int count=0;
            int ancestralCount=0;
            while ((line=br.readLine()).startsWith("##")){
                bw.write(line);
                bw.newLine();
            }
            bw.write(line);
            bw.newLine();
            String te;
            StringBuilder sb=new StringBuilder();
            while ((line=br.readLine())!=null){
                total++;
                temp=PStringUtils.fastSplit(line);
                chr=Integer.parseInt(temp.get(0));
                pos=Integer.parseInt(temp.get(1));
                refAllele=temp.get(3);
                altAllele=temp.get(4);
                ancestralAllele=ancestralTable.get(chr, pos);
                if (!ancestralTable.contains(chr, pos)) continue;
                ancestralCount++;
                if (refAllele.equals(ancestralAllele)){
                    count++;
                    bw.write(line);
                    bw.newLine();
                }else if (altAllele.equals(ancestralAllele)){
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0,3))).append("\t");
                    sb.append(ancestralAllele).append("\t").append(refAllele).append("\t");
                    sb.append(String.join("\t", temp.subList(5,9))).append("\t");
                    for (int i = 9; i < temp.size(); i++) {
                        if (temp.get(i).startsWith("0|0")){
                            te=temp.get(i);
                            sb.append("1|1").append(te.substring(3)).append("\t");
                        }else if (temp.get(i).startsWith("1|1")){
                            te=temp.get(i);
                            sb.append("0|0").append(te.substring(3)).append("\t");
                        }else if (temp.get(i).startsWith("0|1")){
                            te=temp.get(i);
                            sb.append("1|0").append(te.substring(3)).append("\t");
                        }else if (temp.get(i).startsWith("1|0")){
                            te=temp.get(i);
                            sb.append("0|1").append(te.substring(3)).append("\t");
                        }else {
                            sb.append(temp.get(i)).append("\t");
                        }
                    }
                    sb.deleteCharAt(sb.length()-1);
                    bw.write(sb.toString());
                    bw.newLine();
                }
            }
            bw.flush();
            System.out.println(vcfFile.getName()+" total sites: "+total+", ancestral site: "+ancestralCount+", " +
                    "refAllele equal " +
                    "ancestral allele sites: "+ count);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void bulidSMC(String vcfDir, String bedDir, String outDir, String taxaInfo_depthFile, Group group,
                                String outFileAB, String outFileD, String chrSizeFile, String logDir, int sampleSize){
//        vcfDir="/data1/home/daxing/vmap2.1Data/002_vmap2.1RefToAncestral";
//        bedDir="/data1/home/daxing/vmap2.1Data/001_vmap2.1_complementChr";
//        outDir="/Users/xudaxing/Desktop/out";
//        taxaInfo_depthFile="/Users/xudaxing/Desktop/groupSMC_R1.txt";
//        distTaxon="/Users/xudaxing/Desktop/popDistTaxonR1.txt";
//        outFileAB="/Users/xudaxing/Desktop/resAB.sh";
//        outFileD="/Users/xudaxing/Desktop/resD.sh";
//        chrSizeFile="/Users/xudaxing/Desktop/genome.txt";
//        String[] group_AB={"WildEmmer","DomesticatedEmmer","FreeThreshTetraploid","Landrace", "Cultivar"};
//        String[] group_D={"Ae.tauschii","Landrace", "Cultivar"};
//        logDir="/Users/xudaxing/Desktop/log";
        smc(vcfDir, bedDir, outDir,SubgenomeCombination.AB, taxaInfo_depthFile, group,
                outFileAB, chrSizeFile, logDir, sampleSize);
        smc(vcfDir, bedDir, outDir,SubgenomeCombination.D, taxaInfo_depthFile, group,
                outFileD, chrSizeFile, logDir, sampleSize);
    }

    /**
     *
     * @param vcfDir vcf dir
     * @param vcfComplementBedDir vcfComplementBedDir
     * @param subgenomeCombination
     * @param taxaInfo_DepthFile groupFile
     * DomesticatedEmmer	CItr14822
     * DomesticatedEmmer	CItr14824
     * DomesticatedEmmer	CItr14916
     * FreeThreshTetraploid	CItr14892
     * FreeThreshTetraploid	CItr7798
     * @param group group
     * @param shOutFile shOutFile
     */
    private static void smc(String vcfDir, String vcfComplementBedDir, String outDir,
                      SubgenomeCombination subgenomeCombination,
                           String taxaInfo_DepthFile, Group group, String shOutFile, String chrSizeFile, String logDir
            , int sampleSize){
        List<File> vcfFiles=IOUtils.getFileListInDirEndsWith(vcfDir, "gz");
        List<File> vcfComplementBedFiles=IOUtils.getFileListInDirEndsWith(vcfComplementBedDir, "gz");
        Predicate<File> subgenomeP = null;
        TIntArrayList subgenomeChrs = null;
        if (subgenomeCombination.equals(SubgenomeCombination.AB)){
            TIntArrayList ab=new TIntArrayList(WheatLineage.ablineage());
            subgenomeP=f->ab.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.ablineage());
        }else if (subgenomeCombination.equals(SubgenomeCombination.D)){
            TIntArrayList d=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
            subgenomeP=f->d.contains(StringTool.getNumFromString(f.getName()));
            subgenomeChrs=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        }else {
            System.out.println("error your parameter subgenome");
            System.exit(1);
        }
        List<File> subgenomeVcfFiles=vcfFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        List<File> subgenomeVcfComplementBedFiles=vcfComplementBedFiles.stream().filter(subgenomeP).collect(Collectors.toList());
        Multimap<String, String> popmap=getPopMultiMap(taxaInfo_DepthFile, group, sampleSize);
        Map<String, String> popDistTaxonMap=getPopDistTaxonMap(popmap);
        Map<Integer, Integer> chrSizeMap=new HashMap<>();
        try (BufferedReader bufferedReader1=IOTool.getReader(chrSizeFile);
             BufferedWriter bufferedWriter =IOTool.getWriter(shOutFile)) {
            String line;
            List<String> temp;
            while ((line=bufferedReader1.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                chrSizeMap.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)));
            }
            StringBuilder sb;
            List<String> individuals;
            String distTaxon;
            File outFile, logFile;
            GroupType groupType;
            List<String> groupList= group.getGroup(subgenomeCombination);
            for (int i = 0; i < subgenomeVcfFiles.size(); i++) {
                for (String s : groupList) {
                    individuals = new ArrayList<>(popmap.get(s));
                    distTaxon = popDistTaxonMap.get(s);
                    groupType = GroupType.newInstanceFrom(s, group);
                    sb = new StringBuilder();
                    sb.append("nohup smc++ vcf2smc --cores 1 -m ").append(subgenomeVcfComplementBedFiles.get(i)).append(" ");
                    sb.append(subgenomeVcfFiles.get(i)).append(" -l ").append(chrSizeMap.get(subgenomeChrs.get(i))).append(" -d ");
                    sb.append(distTaxon).append(" ").append(distTaxon).append(" ");
                    outFile = new File(outDir, subgenomeVcfFiles.get(i).getName().substring(0, 6) + "." + groupType + ".smc" +
                            ".gz");
                    sb.append(outFile.getAbsolutePath()).append(" ").append(subgenomeChrs.get(i)).append(" ");
                    sb.append(groupType).append(":");
                    for (String individual : individuals) {
                        sb.append(individual).append(",");
                    }
                    sb.deleteCharAt(sb.length()-1).append(" >");
                    logFile=new File(logDir, subgenomeVcfFiles.get(i).getName().substring(0, 6)+
                            "_"+groupType.getGroupAbbrev()+".log");
                    sb.append(logFile.getAbsolutePath()).append(" 2>&1");
                    bufferedWriter.write(sb.toString());
                    bufferedWriter.newLine();
                }
            }
            bufferedWriter.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *
     * @param taxaInfo_DistinguishedLineageFile
     * @param group
     * @param sampleSize vcf2smc亚群的样本大小，不能超过亚群的组合数
     * @return
     */
    private static Multimap<String, String> getPopMultiMap(String taxaInfo_DistinguishedLineageFile,
                                                           Group group, int sampleSize){
        int index= group.equals(Group.Subspecies) ? 3 : 2;
        Set<String> groupSet=RowTableTool.getColumnSet(taxaInfo_DistinguishedLineageFile, index);
        groupSet.remove("NA");
        List<String> groupList=new ArrayList<>(groupSet);
        Collections.sort(groupList);
        Multimap<String, String> res= TreeMultimap.create();
        try (BufferedReader br = IOTool.getReader(taxaInfo_DistinguishedLineageFile)) {
            String line;
            List<String> temp;
            br.readLine();
            Multimap<String, String> popMap=TreeMultimap.create();
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                if (temp.get(index).equals("NA")) continue;
                groupSet.add(temp.get(index));
                popMap.put(temp.get(index), temp.get(0));
            }
            List<String> individualsList;
            StringBuilder sb=new StringBuilder();
            List<int[]> combinationList;
            int[] combiantion;
            GroupType groupType;
            for (String pop:groupList){
                individualsList=new ArrayList<>(popMap.get(pop));
                Collections.sort(individualsList);
                groupType = GroupType.newInstanceFrom(pop, group);
                Iterator<int[]> iterator=CombinatoricsUtils.combinationsIterator(individualsList.size(), 2);
                combinationList=new ArrayList<>();
                while (iterator.hasNext()){
                    combinationList.add(iterator.next());
                }
                Collections.shuffle(combinationList);
                for (int i = 0; i < sampleSize; i++) {
                    combiantion=combinationList.get(i);
                    sb.setLength(0);
                    sb.append(groupType).append("_");
                    sb.append(individualsList.get(combiantion[0])).append("_").append(individualsList.get(combiantion[1]));
                    res.put(pop, sb.toString());
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return res;
    }

    private static Map<String,String> getPopDistTaxonMap(Multimap<String, String> popmap){
        Set<String> keySet=popmap.keySet();
        Map<String,String> res=new HashMap<>();
        List<String> valueList;
        for (String key:keySet){
            valueList=new ArrayList<>(popmap.get(key));
            Collections.shuffle(valueList);
            res.put(key, valueList.get(0));
        }
        return res;
    }

    /**
     * 将smc合并为RefChrSmc
     * @param inputDir
     * @param outDir
     * @param popNumAB
     * @param popNumD
     */
    public static void mergeSMC(String inputDir, String outDir, int popNumAB, int popNumD){
        System.out.println(DateTime.getDateTimeOfNow());
        long start=System.nanoTime();
        mergeSMC_AB(inputDir, outDir, popNumAB);
        mergeSMC_D(inputDir, outDir, popNumD);
        System.out.println("mergeSMC had completed in "+Benchmark.getTimeSpanHours(start)+ " hours");
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void mergeSMC_D(String inputDir, String outDir, int popNumD){
        long start=System.nanoTime();
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        TIntArrayList dindex=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> d=f->dindex.contains(StringTool.getNumFromString(f.getName()));
        List<File> dFiles=files.stream().filter(d).collect(Collectors.toList());
        TIntArrayList[] num=new TIntArrayList[2];
        num[0]=new TIntArrayList();
        num[1]=new TIntArrayList();
        int count=0;
        int index=0;
        boolean flag=false;
        for (int i = 0; i < dFiles.size(); i++) {
            count++;
            if (count>popNumD){
                flag=!flag;
                index=flag ? 1 : 0;
                count=1;
            }
            num[index].add(i);
        }
        IntStream.range(0, num[0].size()).parallel().forEach(e->mergeSMC(dFiles.get(num[0].get(e)), dFiles.get(num[1].get(e)), outDir));
        System.out.println("D subgenome had completed in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
    }

    private static void mergeSMC_AB(String inputDir, String outDir, int popNumAB){
        long start=System.nanoTime();
        List<File> files=IOUtils.getVisibleFileListInDir(inputDir);
        TIntArrayList dindex=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> d=f->dindex.contains(StringTool.getNumFromString(f.getName()));
        List<File> dFiles=files.stream().filter(d.negate()).collect(Collectors.toList());
        TIntArrayList[] num=new TIntArrayList[2];
        num[0]=new TIntArrayList();
        num[1]=new TIntArrayList();
        int count=0;
        int index=0;
        boolean flag=false;
        for (int i = 0; i < dFiles.size(); i++) {
            count++;
            if (count>popNumAB){
                flag=!flag;
                index=flag ? 1 : 0;
                count=1;
            }
            num[index].add(i);
        }
        IntStream.range(0, num[0].size()).parallel().forEach(e->mergeSMC(dFiles.get(num[0].get(e)), dFiles.get(num[1].get(e)), outDir));
        System.out.println("AB subgenome had completed in "+ Benchmark.getTimeSpanMinutes(start)+ " minutes");
    }

    private static void mergeSMC(File file1, File file2, String outDir){
        int chrID=StringTool.getNumFromString(file1.getName().substring(0,6));
        String filename=file1.getName().substring(6);
        String chr= RefV1Utils.getChromosome(chrID, 1);
        File outFile=new File(outDir, "chr" + chr + filename);
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            long numLine1=IOTool.getReader(file1).lines().count();
            int num1=0;
            BufferedReader br1=IOTool.getReader(file1);
            BufferedReader br2=IOTool.getReader(file2);
            String header=br1.readLine();
            num1++;
            bw.write(header);
            bw.newLine();
            String line;
            List<String> temp1, temp2;
            while ((line=br1.readLine())!=null && num1 < numLine1-1){
                num1++;
                bw.write(line);
                bw.newLine();
            }
            temp1=PStringUtils.fastSplit(line, " ");
            br1.close();
            br2.readLine();
            line=br2.readLine();
            temp2=PStringUtils.fastSplit(line, " ");
            int count;
            if (temp1.subList(1,temp1.size()).equals(temp2.subList(1,temp2.size()))){
                count=Integer.parseInt(temp1.get(0))+Integer.parseInt(temp2.get(0));
                bw.write(count+" "+String.join(" ", temp1.subList(1,4)));
                bw.newLine();
            }else {
                bw.write(line);
                bw.newLine();
            }
            while ((line=br2.readLine())!=null){
                bw.write(line);
                bw.newLine();
            }
            br2.close();
            bw.flush();
            System.out.println(outFile.getName()+ " had completed");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
