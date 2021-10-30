package daxing.individualIntrogression;

import com.google.common.collect.Table;
import com.google.common.collect.TreeBasedTable;
import daxing.common.utiles.DateTime;
import daxing.common.utiles.IOTool;
import daxing.common.factors.LoadType;
import daxing.common.table.RowTableTool;
import daxing.load.ancestralSite.ChrSNPAnnoDB;
import daxing.load.ancestralSite.SNPAnnotation;
import pgl.infra.pos.ChrPos;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * 将个体按site记录load
 */
public class IntrogressionDonorBurdenGo {

    public static void start(){
        String exonSNPAnnoDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/002_deleteriousPosLib/002_exonAnnotationByDerivedSift/001_exonSNPAnnotationByChrID";
        String exonVCFDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/003_exon/001_exonVCF";
        String taxa_InfoDB="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String popFdFile="";
        String individualFdDir="";
        SNPAnnotation.MethodCallDeleterious methodCallDeleterious = SNPAnnotation.MethodCallDeleterious.GERP;
        String outDir="/Users/xudaxing/Documents/deleteriousMutation/001_analysis/003_vmap2.1_20200628/004_deleterious/001_triadsSelection/003_derivedSiftPerSite";
        go(exonSNPAnnoDir, exonVCFDir, taxa_InfoDB, popFdFile, individualFdDir, methodCallDeleterious, outDir);
    }

    public static void go(String exonSNPAnnoDir, String exonVCFDir, String taxa_InfoDBFile, String popFdFile,
                          String individualFdDir, SNPAnnotation.MethodCallDeleterious methodCallDeleterious, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        String[] subdir={"001_count","002_countMerge","003_retainLandraceCutivar",
                "004_loadPerSiteGrid","005_introgressionDonorBurden"};
        Path loadPerSiteGridPath= Paths.get(new File(outDir, subdir[3]).getAbsolutePath(), "loadPerSiteGrid.txt.gz");
        Path introgressionDonorBurdenPath=Paths.get(new File(outDir, subdir[4]).getAbsolutePath(), "introgressionDonorBurden.txt.gz");
        for (int i = 0; i < subdir.length; i++) {
            new File(outDir, subdir[i]).mkdir();
        }
        List<File> exonAnnoFiles= IOUtils.getVisibleFileListInDir(exonSNPAnnoDir);
        List<File> exonVCFFiles= IOUtils.getVisibleFileListInDir(exonVCFDir);
        Map<String, File> taxonOutDirMap=getTaxonOutDirMap(taxa_InfoDBFile, new File(outDir, subdir[0]).getAbsolutePath());
        IntStream.range(0, exonVCFFiles.size()).parallel().forEach(e->go(exonAnnoFiles.get(e), exonVCFFiles.get(e),
                taxonOutDirMap, e+1, methodCallDeleterious));
        merge(new File(outDir, subdir[0]).getAbsolutePath(), new File(outDir, subdir[1]).getAbsolutePath());
        retainTreeValidatedLandraceCultivar(new File(outDir, subdir[1]).getAbsolutePath(), taxa_InfoDBFile,
                new File(outDir, subdir[2]).getAbsolutePath());
        IntrogressionDonorBurdenGo.mergeLoadSiteRecordToGrid(new File(outDir, subdir[2]).getAbsolutePath(),
                loadPerSiteGridPath.toString());
        PopulationIndividualFd.writeWindowSize(popFdFile,individualFdDir, taxa_InfoDBFile, loadPerSiteGridPath.toString(),
                introgressionDonorBurdenPath.toString());
        System.out.println(DateTime.getDateTimeOfNow());
    }

    private static void go(File exonSNPAnnoFile, File exonVCFFile,
                           Map<String, File> taxonOutDirMap, int chr, SNPAnnotation.MethodCallDeleterious methodCallDeleterious){
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
            LoadType loadType;
            boolean ifHeter, ifHomozygousDerived;
            Optional<LoadType> optionalLoadType;
            while ((line=br.readLine())!=null){
                temp=PStringUtils.fastSplit(line);
                pos=Integer.parseInt(temp.get(1));
                if (!transcriptDB.contain(chr, pos)) continue;
                if (!transcriptDB.hasAncestral(chr, pos)) continue;
                isSyn=transcriptDB.isSyn(chr, pos);
                isNonsyn=transcriptDB.isNonsyn(chr, pos);
                if (!(isSyn || isNonsyn)) continue;
                isDeleterious=transcriptDB.isDeleterious(chr, pos, methodCallDeleterious);
                isRefAlleleAncestral=transcriptDB.isRefAlleleAncestral(chr, pos);
                for (int i = 0; i < taxonNames.size(); i++) {
                    genotypeList=PStringUtils.fastSplit(temp.get(i+9), ":");
                    genotype=genotypeList.get(0);
                    if (genotype.equals("./.")) continue;
//                    if (genotype.equals("0/1")) continue;
                    depthList=PStringUtils.fastSplit(genotypeList.get(1),",");
                    depth=Integer.parseInt(depthList.get(0))+Integer.parseInt(depthList.get(1));
                    if (depth < 2) continue;
                    optionalLoadType=LoadType.getInstanceFromString(transcriptDB.getVariantType(chr, pos));
                    optionalLoadType.orElseThrow(IllegalArgumentException::new);
                    loadType= optionalLoadType.get();
                    loadType= isDeleterious ? LoadType.Del : loadType;
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
        Map<String, String> taxaTreeValidatedGroupbySubspeciesMap=RowTableTool.getMap(taxa_InfoDB, 0, 10);
        Map<String, String> taxaFdIDMap=RowTableTool.getMap(taxa_InfoDB, 0, 35);
        Predicate<File> landraceP= file -> taxaTreeValidatedGroupbySubspeciesMap.get(PStringUtils.fastSplit(file.getName(), ".").get(0)).equals("Landrace");
        Predicate<File> cultivarP= file -> taxaTreeValidatedGroupbySubspeciesMap.get(PStringUtils.fastSplit(file.getName(), ".").get(0)).equals("Cultivar");
        Predicate<File> predicateLRCL=landraceP.or(cultivarP);
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
     * 将个体load和fd以vcf的形式输出
     * @param inputAddFdDir
     * @param outFile
     */
    public static void mergeLoadSiteRecordToGrid(String inputAddFdDir, String outFile){
        List<File> files=IOUtils.getVisibleFileListInDir(inputAddFdDir);
        long start=System.nanoTime();
        Table<Integer, Integer, String> chrPosLoadTypeTable= TreeBasedTable.create();
        BufferedReader br;
        try {
            String line, loadType;
            List<String> temp;
            int chrID, pos;
            for (int i = 0; i < files.size(); i++) {
                long start1=System.nanoTime();
                br=IOTool.getReader(files.get(i));
                br.readLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    loadType=temp.get(2);
                    chrPosLoadTypeTable.put(chrID, pos, loadType);
                }
                br.close();
                System.out.println(files.get(i).getName()+ " extract completed in "+Benchmark.getTimeSpanMilliseconds(start1)+ " ms");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("extract completed in "+ Benchmark.getTimeSpanSeconds(start)+ " s");
        IndividualLoadToGrid individualLoadToGrid =new IndividualLoadToGrid(chrPosLoadTypeTable);
        for (int i = 0; i < files.size(); i++) {
            individualLoadToGrid.addIndividual(files.get(i).getAbsolutePath());
        }
        individualLoadToGrid.write(outFile);
    }

    private static class ChrPosLoadType extends ChrPos{

        LoadType loadType;

        public ChrPosLoadType(int chrID, int pos, String loadType){
            super((short) chrID, pos);
            this.loadType=LoadType.valueOf(loadType);
        }

        public LoadType getLoadType() {
            return loadType;
        }
    }

    private static class IndividualLoadToGrid {

        List<ChrPosLoadType> chrPosLoadTypes;
        List<LoadSiteRecord[]> individualSummary;
        List<String> taxonNames;

        private IndividualLoadToGrid(Table<Integer, Integer, String> chrPosLoadTypeTable){
            List<ChrPosLoadType> chrPosLoadTypes=new ArrayList<>(chrPosLoadTypeTable.size());
            ChrPosLoadType chrPosLoadType;
            Map<Integer, Map<Integer, String>> rowMapMap=chrPosLoadTypeTable.rowMap();
            for (Map.Entry<Integer, Map<Integer, String>> entry1 : rowMapMap.entrySet()){
                for (Map.Entry<Integer, String> entry2 : entry1.getValue().entrySet()){
                    chrPosLoadType=new ChrPosLoadType(entry1.getKey(), entry2.getKey(), entry2.getValue());
                    chrPosLoadTypes.add(chrPosLoadType);
                }
            }
            Collections.sort(chrPosLoadTypes);
            this.chrPosLoadTypes=chrPosLoadTypes;
            individualSummary=new ArrayList<>();
            taxonNames=new ArrayList<>();
        }

        private int getSNPNum(){
            return this.chrPosLoadTypes.size();
        }

        private int getTaxonNum(){
            return this.individualSummary.size();
        }

        private int binarySearch(ChrPos chrPos){
            return Collections.binarySearch(chrPosLoadTypes, chrPos);
        }

        private ChrPosLoadType getChrPosLoadType(int snpIndex){
            return this.chrPosLoadTypes.get(snpIndex);
        }

        private LoadSiteRecord getIndividualLoadFdRecord(int snpIndex, int taxonIndex){
            return this.individualSummary.get(taxonIndex)[snpIndex];
        }

        private void addIndividual(String individualLoadFdFile){
            String taxaName=new File(individualLoadFdFile).getName().substring(0, 7);
            LoadSiteRecord[] loadSiteRecords =new LoadSiteRecord[this.getSNPNum()];
            Arrays.fill(loadSiteRecords, LoadSiteRecord.getDefault());
            try (BufferedReader br = IOTool.getReader(individualLoadFdFile)) {
                br.readLine();
                String line;
                List<String> temp;
                int chrID, pos, index;
                ChrPosLoadType chrPosLoadType;
                int[] ifHeterIfHomozygousDerivedIfIntrogression;
                LoadSiteRecord loadSiteRecord;
                while ((line= br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chrID=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    chrPosLoadType=new ChrPosLoadType(chrID, pos, temp.get(2));
                    index=binarySearch(chrPosLoadType);
                    ifHeterIfHomozygousDerivedIfIntrogression=new int[2];
                    ifHeterIfHomozygousDerivedIfIntrogression[0]=Integer.parseInt(temp.get(3));
                    ifHeterIfHomozygousDerivedIfIntrogression[1]=Integer.parseInt(temp.get(4));
                    loadSiteRecord =new LoadSiteRecord(ifHeterIfHomozygousDerivedIfIntrogression);
                    loadSiteRecords[index]= loadSiteRecord;
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
            individualSummary.add(loadSiteRecords);
            taxonNames.add(taxaName);
        }

        private void write(String outFile){
            try (BufferedWriter bw = IOTool.getWriter(outFile)) {
                StringBuilder sb=new StringBuilder();
                sb.append("Chr\tPos\tLoadType\t").append(String.join("\t", this.taxonNames));
                bw.write(sb.toString());
                bw.newLine();
                int snpNum=this.getSNPNum();
                int taxonNum=this.getTaxonNum();
                ChrPosLoadType chrPosLoadType;
                LoadSiteRecord loadSiteRecord;
                for (int i = 0; i < snpNum; i++) {
                    chrPosLoadType=this.getChrPosLoadType(i);
                    sb.setLength(0);
                    sb.append(chrPosLoadType.getChromosome()).append("\t").append(chrPosLoadType.getPosition()).append("\t");
                    sb.append(chrPosLoadType.getLoadType()).append("\t");
                    for (int j = 0; j < taxonNum; j++) {
                        loadSiteRecord =this.getIndividualLoadFdRecord(i, j);
                        sb.append(loadSiteRecord.toString()).append("\t");
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
    }

    private static class LoadSiteRecord {

        int[] ifHeterIfHomozygousDerivedIfIntrogression;

        public static final LoadSiteRecord DEFAULT = getDefault();

        private static final LoadSiteRecord getDefault(){
            int[] ifHeterIfHomozygousDerivedIfIntrogression={-1,-1};
            return new LoadSiteRecord(ifHeterIfHomozygousDerivedIfIntrogression);
        }

        private LoadSiteRecord(int[] ifHeterIfHomozygousDerivedIfIntrogression){
            this.ifHeterIfHomozygousDerivedIfIntrogression=ifHeterIfHomozygousDerivedIfIntrogression;
        }

        @Override
        public String toString() {
            StringBuilder sb=new StringBuilder();
            for (int i = 0; i < ifHeterIfHomozygousDerivedIfIntrogression.length; i++) {
                sb.append(ifHeterIfHomozygousDerivedIfIntrogression[i]).append(",");
            }
            sb.deleteCharAt(sb.length()-1);
            return sb.toString();
        }
    }
}
