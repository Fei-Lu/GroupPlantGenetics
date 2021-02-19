package daxing.applets;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.WheatLineage;
import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.math3.util.CombinatoricsUtils;
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

public class Ne {

    static String[] groupBySubspecies={"Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid","Ae.tauschii",
            "LR_EU","LR_EA","Cultivar"};
    static String[] groupBySubspeciesAB={"Wild_emmer","Domesticated_emmer","Free_threshing_tetraploid",
            "LR_EU","LR_EA","Cultivar"};
    static String[] groupBySubspeciesD={"Ae.tauschii","LR_EU","LR_EA","Cultivar"};

    public enum Group{
        WE("Wild_emmer"),DE("Domesticated_emmer"),FT("Free_threshing_tetraploid"),AT("Ae.tauschii"),
        LR_EU("LR_EU"),LR_EA("LR_EA"), CL("Cultivar");

        String group;

        Group(String group) {
            this.group=group;
        }

        public static Group newInstanceFrom(String group){
            switch (group){
                case "Wild_emmer":
                    return WE;
                case "Domesticated_emmer":
                    return DE;
                case "Free_threshing_tetraploid":
                    return FT;
                case "Ae.tauschii":
                    return AT;
                case "LR_EU":
                    return LR_EU;
                case "LR_EA":
                    return LR_EA;
                case "Cultivar":
                    return CL;
            }
            return null;
        }
    }

    public static void syntheticPseudoDiploid(String vmap2InputDir, String taxaInfo_depthFile, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> files = IOUtils.getFileListInDirEndsWith(vmap2InputDir, "gz");
        List<File> abFiles=getFiles(files, "AB");
        List<File> dFiles=getFiles(files, "D");
        String[] abOutFile= abFiles.stream().map(File::getName).
                map(s -> s.replaceAll(".vcf.gz", ".DistinguishedLineages.vcf.gz")).toArray(String[]::new);
        String[] dOutFile= dFiles.stream().map(File::getName).
                map(s -> s.replaceAll(".vcf.gz", ".DistinguishedLineages.vcf.gz")).toArray(String[]::new);
        Map<String,String> taxaGroupABMap=getTaxaGroupMap(taxaInfo_depthFile, groupBySubspeciesAB);
        Map<String,String> taxaGroupDMap=getTaxaGroupMap(taxaInfo_depthFile, groupBySubspeciesD);
        IntStream.range(0, abFiles.size()).parallel().forEach(e->syntheticPseudoDiploid(abFiles.get(e),taxaGroupABMap,
                new File(outDir, abOutFile[e])));
        System.out.println("AB subgenome completed");
        System.out.println(DateTime.getDateTimeOfNow());
        IntStream.range(0, dFiles.size()).parallel().forEach(e->syntheticPseudoDiploid(dFiles.get(e),taxaGroupDMap,
                new File(outDir, dOutFile[e])));
        System.out.println("D subgemome completed");
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void syntheticPseudoDiploid(File vmap2File, Map<String, String> taxaGroupMap, File outFile){
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
            Group group;
            for (int i = 0; i < iteratorCombinations.length; i++) {
                group= Group.newInstanceFrom(groupList.get(i));
                Iterator<int[]> iterator=iteratorCombinations[i];
                while (iterator.hasNext()){
                    int[] combination= iterator.next();
                    indexArray=new int[2];
                    indexArray[0]=distinguishedLineageIndexArray[i].get(combination[0]);
                    indexArray[1]=distinguishedLineageIndexArray[i].get(combination[1]);
                    sb.append(group.name()).append("_").append(temp.get(indexArray[0])).append("_").append(temp.get(indexArray[1])).append("\t");
                    indexListArray[i].add(indexArray);
                }
            }
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
        } catch (IOException e) {
            e.printStackTrace();
        }
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

    private static Map<String, String> getTaxaGroupMap(String taxaInfo_depthFile, String[] groupBySubspecies) {
        List<String> group= Arrays.asList(groupBySubspecies);
        RowTableTool<String> tableTool= new RowTableTool<>(taxaInfo_depthFile);
        Predicate<List<String>> predicate= l->group.contains(l.get(2));
        tableTool.removeIf(predicate.negate());
        return tableTool.getHashMap(0,2);
    }

    /**
     *
     * @param files files
     * @param abOrD AB or D
     * @return files
     */
    private static List<File> getFiles(List<File> files, String abOrD){
        int[] ab=WheatLineage.ablineage();
        int[] d=WheatLineage.D.getChrID();
        TIntArrayList abList=new TIntArrayList(ab);
        TIntArrayList dList = new TIntArrayList(d);
        switch (abOrD){
            case "AB":
                return files.stream().filter(file -> abList.contains(Integer.parseInt(file.getName().substring(3,6)))).collect(Collectors.toList());
            case "D":
                return files.stream().filter(file -> dList.contains(Integer.parseInt(file.getName().substring(3,6)))).collect(Collectors.toList());

        }
        return null;
    }

    public static Table<Integer, Integer, String> getAncestral(String ancestralFile){
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
        String[] outNames=vcfFiles.stream().map(File::getName).map(s->s.replaceAll("\\.vcf",".RefToAncestral" +
                ".vcf")).toArray(String[]::new);
        IntStream.range(0, ancestralFiles.size()).parallel().forEach(e->transformRefToAncestral(vcfFiles.get(e),
                ancestralFiles.get(e), new File(outDir, outNames[e])));
        System.out.println(DateTime.getDateTimeOfNow());
    }

    public static void transformRefToAncestral(File vcfFile, File ancestralFile, File outFile){
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



}
