package daxing.individualIntrogression;

import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.NumberTool;
import daxing.common.RowTableTool;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

/**
 * 将个体load分为introgression区和nonintrogression区，introgression区的load是否低于nonintrogression区的load
 */
public class FdLoad {

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_deleteriousPosLib/002_exonAnnotationByDerivedSift/001_exonSNPAnnotationByChrID";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/003_exon/001_exonVCF";
        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String individualFdDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/005_introgression/006_fdResByIndividual/002_fdByIndividual";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftPerSite";
        go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, individualFdDir, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile, String individualFdDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_retainLandraceCutivarNotIncludeIndianDwarf", "004_addFd",
                "005_individualLoadFd"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
//        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
//                taxonOutDirMap, e+1));
//        merge(new File(outDir, subdir[0]).getAbsolutePath(), new File(outDir, subdir[1]).getAbsolutePath());
//        retainTreeValidatedLandraceCultivar(new File(outDir, subdir[1]).getAbsolutePath(), taxa_InfoDBFile,
//                new File(outDir, subdir[2]).getAbsolutePath());
        addFd(individualFdDir, new File(outDir, subdir[2]).getAbsolutePath(), new File(outDir, subdir[3]).getAbsolutePath());
//        summaryIndividualLoadFd(new File(outDir, subdir[3]).getAbsolutePath(),
//                new File(outDir, subdir[4]).getAbsolutePath(), IndividualChrPosLoad.LoadType.Del);
//        summaryIndividualLoadFd(new File(outDir, subdir[3]).getAbsolutePath(),
//                new File(outDir, subdir[4]).getAbsolutePath(), IndividualChrPosLoad.LoadType.Non);
//        summaryIndividualLoadFd(new File(outDir, subdir[3]).getAbsolutePath(),
//                new File(outDir, subdir[4]).getAbsolutePath(), IndividualChrPosLoad.LoadType.Syn);
        System.out.println(DateTime.getDateTimeOfNow());
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
            IndividualChrPosLoad[] taxonLoads=new IndividualChrPosLoad[taxonNames.size()];
            for (int i = 0; i < taxonLoads.length; i++) {
                taxonLoads[i]=new IndividualChrPosLoad(taxonNames.get(i), chr);
            }
            int pos, depth;
            String genotype;
            List<String> genotypeList, depthList;
            boolean isRefAlleleAncestral=true;
            boolean isSyn, isNonsyn, isDeleterious;
            byte genotypeByte;
            IndividualChrPosLoad.LoadType loadType;
            boolean ifHeter, ifHomozygousDerived;
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
                for (int i = 0; i < taxonNames.size(); i++) {
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
//                    if (genotype.equals("0/1")) continue;
                    depthList=PStringUtils.fastSplit(genotypeList.get(1),",");
                    depth=Integer.parseInt(depthList.get(0))+Integer.parseInt(depthList.get(1));
                    if (depth < 2) continue;
                    loadType= IndividualChrPosLoad.LoadType.newInstanceFrom(transcriptDB.getVariantType(chr, pos));
                    loadType= isDeleterious ? IndividualChrPosLoad.LoadType.Del : loadType;
                    genotypeByte= IndividualChrPosLoad.caculateGenotype(genotype, isRefAlleleAncestral);
                    ifHeter = genotypeByte==2 ? true : false;
                    ifHomozygousDerived = genotypeByte==1 ? true : false;
                    taxonLoads[i].addChrPos(chr, pos, loadType, ifHeter, ifHomozygousDerived);
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
        String outName, header, line;
        BufferedReader[] brs;
        BufferedWriter bw;
        try {
            for (int i = 0; i < dirs.size(); i++) {
                temp= IOUtils.getVisibleFileListInDir(dirs.get(i).getAbsolutePath());
                brs=new BufferedReader[temp.size()];
                outName= PStringUtils.fastSplit(temp.get(0).getName(), ".").get(1);
                for (int j = 0; j < brs.length; j++) {
                    brs[j]=IOTool.getReader(temp.get(j));
                }
                bw=IOTool.getWriter(new File(outDir, outName+".txt.gz"));
                header=brs[0].readLine();
                bw.write(header);
                bw.newLine();
                while ((line=brs[0].readLine())!=null){
                    bw.write(line);
                    bw.newLine();
                }
                brs[0].close();
                for (int j = 1; j < brs.length; j++) {
                    brs[j].readLine();
                    while ((line=brs[j].readLine())!=null){
                        bw.write(line);
                        bw.newLine();
                    }
                    brs[j].close();
                }
                bw.flush();
                bw.close();
                System.out.println("count merge finished: "+outName);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * retain TreeValidated Landrace Cultivar, do not include indian dwarf
     * @param countMergeDir
     * @param taxa_InfoDB
     * @param outDir
     */
    public static void retainTreeValidatedLandraceCultivar(String countMergeDir, String taxa_InfoDB, String outDir){
        List<File> files=IOUtils.getVisibleFileListInDir(countMergeDir);
        Map<String, String> taxaTreeValidatedGroupbySubspeciesMap=RowTableTool.getMap(taxa_InfoDB, 0, 15);
        Map<String, String> taxaFdIDMap=RowTableTool.getMap(taxa_InfoDB, 0, 23);
        Map<String, String> taxaFdBySubspeciesMap=RowTableTool.getMap(taxa_InfoDB, 0, 25);
        Predicate<File> landraceP= file -> taxaTreeValidatedGroupbySubspeciesMap.get(PStringUtils.fastSplit(file.getName(), ".").get(0)).equals("Landrace");
        Predicate<File> cultivarP= file -> taxaTreeValidatedGroupbySubspeciesMap.get(PStringUtils.fastSplit(file.getName(), ".").get(0)).equals("Cultivar");
        Predicate<File> indianDwarfP= file -> taxaFdBySubspeciesMap.get(PStringUtils.fastSplit(file.getName(),".").get(0)).equals("IndianDwarfWheat");
        Predicate<File> predicateLRCL=landraceP.or(cultivarP).and(indianDwarfP.negate());
        List<File> landraceCultivarFiles=files.stream().filter(predicateLRCL).collect(Collectors.toList());
        String[] outFileFdID= landraceCultivarFiles.stream().map(File::getName).map(f->taxaFdIDMap.get(PStringUtils.fastSplit(f,
                        ".").get(0))).toArray(String[]::new);
        try {
            for (int i = 0; i < landraceCultivarFiles.size(); i++) {
                Files.copy(landraceCultivarFiles.get(i).toPath(), new File(outDir,outFileFdID[i]+".del.load.perSite" +
                        ".txt.gz").toPath());
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     * add IfIntrogression
     * @param fdDir
     * @param lrClDir
     * @param outDir
     */
    public static void addFd(String fdDir, String lrClDir, String outDir){
        List<File> fdFiles=IOUtils.getVisibleFileListInDir(fdDir);
        Map<Range, String>[] fdFilesRangeP3Map=new Map[fdFiles.size()];
        for (int i = 0; i < fdFiles.size(); i++) {
            fdFilesRangeP3Map[i]=getIndividualRangeList(fdFiles.get(i).getAbsolutePath());
        }
        List<File> taxaLoadFiles=IOUtils.getVisibleFileListInDir(lrClDir);
        if (taxaLoadFiles.size()!=fdFiles.size()){
            System.out.println("check your fdDir and lrClDir");
            System.exit(1);
        }
        BufferedReader br;
        BufferedWriter bw;
        try {
            String line, outFile;
            List<String> temp;
            int chrID, pos, index;
            Range range;
            List<Range> ranges;
            Map<Range, String> fdRangeP3map;
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < taxaLoadFiles.size(); i++) {
                br=IOTool.getReader(taxaLoadFiles.get(i));
                outFile=taxaLoadFiles.get(i).getName().replaceAll(".csv.gz",".fdLoad.txt.gz");
                bw=IOTool.getWriter(new File(outDir, outFile));
                br.readLine();
                bw.write("Chr\tPos\tLoadType\tIfHeter\tIfHomozygousDerived\tIfIntrogression");
                bw.newLine();
                fdRangeP3map=fdFilesRangeP3Map[i];
                ranges=new ArrayList<>(fdRangeP3map.keySet());
                Collections.sort(ranges);
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    range=new Range(chrID, pos, pos+1);
                    index=Collections.binarySearch(ranges, range);
                    if (index < -2){
                        index=-index-2;
                        if (ranges.get(index).isOverlap(range) && ranges.get(index-1).isOverlap(range)){
                            sb.setLength(0);
                            temp.add("1");
                            sb.append(fdRangeP3map.get(ranges.get(index-1))).append("_").append(fdRangeP3map.get(ranges.get(index)));
                            temp.add(sb.toString());
                        }else if (ranges.get(index).isOverlap(range) && !(ranges.get(index-1).isOverlap(range))){
                            temp.add("1");
                            temp.add(fdRangeP3map.get(ranges.get(index)));
                        }else {
                            temp.add("0");
                            temp.add("."); //. no introgression
                        }
                    }else if (index == -2){
                        index=-index-2;
                        if (ranges.get(index).isOverlap(range)){
                            temp.add("1");
                            temp.add(fdRangeP3map.get(ranges.get(index)));
                        }else {
                            temp.add("0");
                            temp.add("."); //. no introgression
                        }
                    } else if (index == -1){
                        temp.add("0");
                        temp.add("."); //. no introgression
                    }else if (index > 0){
                        if (ranges.get(index-1).isOverlap(range)){
                            temp.add("1");
                            sb.setLength(0);
                            sb.append(fdRangeP3map.get(ranges.get(index-1))).append("_").append(fdRangeP3map.get(ranges.get(index)));
                            temp.add(sb.toString());
                        }else {
                            temp.add("1");
                            temp.add(fdRangeP3map.get(ranges.get(index)));
                        }
                    }else {
                        temp.add("1");
                        temp.add(fdRangeP3map.get(ranges.get(index)));
                    }
                    bw.write(String.join("\t", temp));
                    bw.newLine();
                }
                br.close();
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    private static Map<Range, String> getIndividualRangeList(String fdFile){
        List<Range> ranges=new ArrayList<>();
        List<String> p3List=new ArrayList<>();
        Map<Range, String> resMap=new HashMap<>();
        Range range;
        try (BufferedReader br = IOTool.getReader(fdFile)) {
            br.readLine();
            String line;
            List<String> temp;
            String chr;
            int refPos, refStart, refEnd, chrID, start, end;
            double fd;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line,",");
                fd=Double.parseDouble(temp.get(9));
                if ( fd < 1) continue;
                chr=temp.get(0);
                refPos=Integer.parseInt(temp.get(1));
                chrID= RefV1Utils.getChrID(chr, refPos);
                refStart=Integer.parseInt(temp.get(1));
                refEnd=Integer.parseInt(temp.get(2));
                start= RefV1Utils.getPosOnChrID(chr,refStart);
                end=RefV1Utils.getPosOnChrID(chr, refEnd);
                range=new Range(chrID, start, end+1);
                ranges.add(range);
                p3List.add(temp.get(11).substring(0,2));
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        range=ranges.get(0);
        String p3=p3List.get(0);
        for (int i = 1; i < ranges.size(); i++) {
            if (range.isOverlap(ranges.get(i)) && p3.equals(p3List.get(i))){
                range=new Range(range.chr, range.start, ranges.get(i).end);
            }else {
                resMap.put(range, p3);
                range=ranges.get(i);
                p3=p3List.get(i);
            }
        }
        resMap.put(range, p3);
        return resMap;
    }

    /**
     * outFile
     * Taxa	Sub	IfFd	Load
     * CL001	A	0	0.0449
     * CL001	A	1	0.0687
     * CL001	B	0	0.0558
     * CL001	B	1	0.0708
     * CL001	D	0	0.0708
     * CL001	D	1	0.0265
     * @param addFdDir
     * @param summaryOutDir
     * @param loadType
     */
    public static void summaryIndividualLoadFd(String addFdDir, String summaryOutDir, IndividualChrPosLoad.LoadType loadType){
        if (loadType.equals("Non")){
            summaryIndividualLoadFdNon(addFdDir, summaryOutDir);
        }else {
            List<File> files=IOTool.getVisibleDir(addFdDir);
            BufferedReader br;
            BufferedWriter bw;
            try {
                bw = IOTool.getWriter(new File(summaryOutDir, loadType.name()+".load.txt.gz"));
                bw.write("Taxa\tSub\tIfFd\tLoad");
                bw.newLine();
                String line, taxaName, sub, group;
                int chrID, ifFd, groupIndex, total;
                StringBuilder sb=new StringBuilder();
                List<String> temp;
                String[] groupArray={"A0","A1","B0","B1","D0","D1"};
                double[] loadSum;
                for (int i = 0; i < files.size(); i++) {
                    taxaName=PStringUtils.fastSplit(files.get(i).getName(), ".").get(0);
                    br=IOTool.getReader(files.get(i));
                    br.readLine();
                    loadSum=new double[groupArray.length];
                    total=0;
                    while ((line=br.readLine())!=null){
                        temp=PStringUtils.fastSplit(line);
                        if (!temp.get(2).equals(loadType.name())) continue;
                        total++;
                        chrID=Integer.parseInt(temp.get(0));
                        sub=RefV1Utils.getSubgenomeFromChrID(chrID);
                        ifFd=Integer.parseInt(temp.get(5));
                        sb.setLength(0);
                        sb.append(sub).append(ifFd);
                        group=sb.toString();
                        groupIndex=Arrays.binarySearch(groupArray, group);
                        loadSum[groupIndex]+=Double.parseDouble(temp.get(3))*0.5+Double.parseDouble(temp.get(4));
                    }
                    br.close();
                    for (int j = 0; j < groupArray.length; j++) {
                        sb.setLength(0);
                        group=groupArray[j];
                        sb.append(taxaName).append("\t").append(group.substring(0,1)).append("\t");
                        sb.append(group.substring(1,2)).append("\t").append(NumberTool.format(loadSum[j]/total, 4));
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
    }

    public static void summaryIndividualLoadFdNon(String addFdDir, String summaryOutDir){
        List<File> files=IOTool.getVisibleDir(addFdDir);
        BufferedReader br;
        BufferedWriter bw;
        try {
            bw = IOTool.getWriter(new File(summaryOutDir, "non.load.txt.gz"));
            bw.write("Taxa\tSub\tIfFd\tLoad");
            bw.newLine();
            String line, taxaName, sub, group;
            int chrID, ifFd, groupIndex, total;
            StringBuilder sb=new StringBuilder();
            List<String> temp;
            String[] groupArray={"A0","A1","B0","B1","D0","D1"};
            double[] loadSum;
            for (int i = 0; i < files.size(); i++) {
                taxaName=PStringUtils.fastSplit(files.get(i).getName(), ".").get(0);
                br=IOTool.getReader(files.get(i));
                br.readLine();
                loadSum=new double[groupArray.length];
                total=0;
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    if (temp.get(2).equals("Syn")) continue;
                    total++;
                    chrID=Integer.parseInt(temp.get(0));
                    sub=RefV1Utils.getSubgenomeFromChrID(chrID);
                    ifFd=Integer.parseInt(temp.get(5));
                    sb.setLength(0);
                    sb.append(sub).append(ifFd);
                    group=sb.toString();
                    groupIndex=Arrays.binarySearch(groupArray, group);
                    loadSum[groupIndex]+=Double.parseDouble(temp.get(3))*0.5+Double.parseDouble(temp.get(4));
                }
                br.close();
                for (int j = 0; j < groupArray.length; j++) {
                    sb.setLength(0);
                    group=groupArray[j];
                    sb.append(taxaName).append("\t").append(group.substring(0,1)).append("\t");
                    sb.append(group.substring(1,2)).append("\t").append(NumberTool.format(loadSum[j]/total, 4));
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

}
