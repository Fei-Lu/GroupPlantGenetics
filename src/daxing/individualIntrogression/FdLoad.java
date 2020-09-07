package daxing.individualIntrogression;

import daxing.common.DateTime;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import pgl.infra.range.Range;
import pgl.infra.table.RowTable;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.*;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Predicate;
import java.util.stream.Collectors;

public class FdLoad {

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_deleteriousPosLib/002_exonAnnotationByDerivedSift/001_exonSNPAnnotationByChrID";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/003_exon/001_exonVCF";
        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String triadFile="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/001_vmap2.1Before20200525/002_analysis/014_deleterious/triadGenes1.1_cdsLen_geneHC.txt";
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftPerSite";
        go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, triadFile, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile, String triadFile, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_retainTreeValidatedAndIndianDwarf"};
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
//        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
//                taxonOutDirMap, e+1));
//        merge(new File(outDir, subdir[0]).getAbsolutePath(), taxa_InfoDBFile, new File(outDir, subdir[1]).getAbsolutePath());
        retainTreeValidatedLandraceCultivar(new File(outDir, subdir[1]).getAbsolutePath(), taxa_InfoDBFile,
                new File(outDir, subdir[2]).getAbsolutePath());
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

    public static void merge(String inputDir, String taxa_InfoDB, String outDir){
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

    public static void addFd(String fdDir, String taxa_InfoDB,  String countMergeDir, String outDir){
        List<File> fdFiles=IOUtils.getVisibleFileListInDir(fdDir);
        List<Range>[] fdFilesRange=new List[fdFiles.size()];
        List<String> fdIDList=new ArrayList<>();
        for (int i = 0; i < fdFiles.size(); i++) {
            fdFilesRange[i]=getIndividualRangeList(fdFiles.get(i).getAbsolutePath());
            fdIDList.add(fdFiles.get(i).getName().substring(0, 5));
        }
        Map<String,String> taxaFdIDMap= RowTableTool.getMap(taxa_InfoDB, 0, 23);
        List<File> taxaFiles=IOUtils.getDirListInDir(countMergeDir);
        Predicate<File> inFdIDList=
                f -> fdIDList.contains(taxaFdIDMap.get(PStringUtils.fastSplit(f.getName(), ".").get(0)));
        List<File> hexaploidFiles=taxaFiles.stream().filter(inFdIDList).collect(Collectors.toList());


    }

    private static List<Range> getIndividualRangeList(String fdFile){
        List<Range> ranges=new ArrayList<>();
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
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        List<Range> res=new ArrayList<>();
        range=ranges.get(0);
        for (int i = 1; i < ranges.size(); i++) {
            if (range.isOverlap(ranges.get(i))){
                range=new Range(range.chr, range.start, ranges.get(i).end);
            }else {
                res.add(range);
                range=ranges.get(i);
            }
        }
        res.add(range);
        return res;
    }

}
