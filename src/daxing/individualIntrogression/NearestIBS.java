package daxing.individualIntrogression;

import daxing.common.*;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOFileFormat;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class NearestIBS {

    /**
     * 并行处理每个染色体
     * @param vmap2VCFDir
     * @param taxa_InfoDB
     * @param fdResDir
     */
    public static void calculateNearestFd(String vmap2VCFDir, String taxa_InfoDB, String fdResDir, String fdOutDir,
                                          int index){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> vmap2FileList= IOUtils.getFileListInDirEndsWith(vmap2VCFDir, "gz");
        Predicate<File> hidden=File::isHidden;
        List<File> vmap2Files=vmap2FileList.stream().filter(hidden.negate()).collect(Collectors.toList());
        String chr;
        TIntArrayList abChrs=new TIntArrayList(WheatLineage.ablineage());
        TIntArrayList dChrs=new TIntArrayList(WheatLineage.valueOf("D").getChrID());
        Predicate<File> abP=file -> abChrs.contains(StringTool.getNumFromString(file.getName().substring(3,6)));
        Predicate<File> dP=file -> dChrs.contains(StringTool.getNumFromString(file.getName().substring(3,6)));
        List<File> abFiles=vmap2Files.stream().filter(abP).sorted().collect(Collectors.toList());
        List<File> dFiles=vmap2Files.stream().filter(dP).sorted().collect(Collectors.toList());
        List<File> fdFiles=IOUtils.getVisibleFileListInDir(fdResDir);
        List<File>[] chrFdABFiles=new List[14];
        List<File>[] chrFdDFiles=new List[7];
        for (int i = 0; i < chrFdABFiles.length; i++) {
            chrFdABFiles[i]=new ArrayList<>();
        }
        for (int i = 0; i < chrFdDFiles.length; i++) {
            chrFdDFiles[i]=new ArrayList<>();
        }
        List<String> chrAB= WheatLineage.abLineage();
        List<String> chrD=WheatLineage.valueOf("D").getChr();
        int chrABIndex=-1, chrDIndex=-1;
        for (int i = 0; i < fdFiles.size(); i++) {
            chr=fdFiles.get(i).getName().substring(3,5);
            chrABIndex= Collections.binarySearch(chrAB, chr);
            if (chrABIndex > -1){
                chrFdABFiles[chrABIndex].add(fdFiles.get(i));
            }else {
                chrDIndex=Collections.binarySearch(chrD, chr);
                chrFdDFiles[chrDIndex].add(fdFiles.get(i));
            }
        }
        RowTableTool<String> taxonTable=new RowTableTool<>(taxa_InfoDB);
        Map<String, String> taxonMap= taxonTable.getHashMap(35,0);
        GenotypeGrid genotypeGridA, genotypeGridB, genotypeGrid;
//        System.out.println("----------- Start calculate: "+chrAB.get(index)+" -----------");
//        genotypeGridA=new GenotypeGrid(abFiles.get(2*index).getAbsolutePath(), GenoIOFormat.VCF_GZ);
//        genotypeGridB=new GenotypeGrid(abFiles.get(2*index+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
//        genotypeGrid= GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
//        genotypeGrid.sortByTaxa();
//        calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdABFiles[index], fdOutDir);
//        System.out.println("----------- finished: "+chrAB.get(index)+" -----------");
        System.out.println("----------- Start calculate: "+chrD.get(index)+" -----------");
        genotypeGridA=new GenotypeGrid(dFiles.get(2*index).getAbsolutePath(), GenoIOFormat.VCF_GZ);
        genotypeGridB=new GenotypeGrid(dFiles.get(2*index+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
        genotypeGrid= GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
        genotypeGrid.sortByTaxa();
        calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdDFiles[index], fdOutDir);
        System.out.println("----------- finished: "+chrD.get(index)+" -----------");
//        for (int i = 0; i < chrFdABFiles.length; i++) {
//            System.out.println("----------- Start calculate: "+chrAB.get(i)+" -----------");
//            genotypeGridA=new GenotypeGrid(abFiles.get(2*i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
//            genotypeGridB=new GenotypeGrid(abFiles.get(2*i+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
//            genotypeGrid= GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
//            genotypeGrid.sortByTaxa();
//            calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdABFiles[i], fdOutDir);
//            System.out.println("----------- finished: "+chrAB.get(i)+" -----------");
//        }
//        for (int i = 0; i < chrFdDFiles.length; i++) {
//            System.out.println("----------- Start calculate: "+chrD.get(i)+" -----------");
//            genotypeGridA=new GenotypeGrid(dFiles.get(2*i).getAbsolutePath(),GenoIOFormat.VCF_GZ);
//            genotypeGridB=new GenotypeGrid(dFiles.get(2*i+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
//            genotypeGrid=GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
//            genotypeGrid.sortByTaxa();
//            calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdDFiles[i],fdOutDir);
//            System.out.println("----------- finished: "+chrD.get(i)+" -----------");
//        }
        System.out.println(DateTime.getDateTimeOfNow());
    }

    /**
     * 在每个染色体内部，并行处理每个p2
     * @param genotypeGrid
     * @param taxonMap
     * @param chrFdFiles
     * @param outDir
     */
    private static void calculateNearestFdCByTaxon(GenotypeGrid genotypeGrid, Map<String
            , String> taxonMap, List<File> chrFdFiles, String outDir){
        genotypeGrid.sortByTaxa();
        String p2, p3;
        Set<String> p2Set=new HashSet<>(), p3Set=new HashSet<>();
        for (int i = 0; i < chrFdFiles.size(); i++) {
            p2=chrFdFiles.get(i).getName().substring(14,21);
            p3=chrFdFiles.get(i).getName().substring(22,29);
            p2Set.add(p2);
            p3Set.add(p3);
        }
        List<String> p2List=new ArrayList<>(p2Set);
        List<String> p3List=new ArrayList<>(p3Set);
        Collections.sort(p2List);
        Collections.sort(p3List);
        List<File>[] taxonFdFiles=new List[p2List.size()];
        for (int i = 0; i < taxonFdFiles.length; i++) {
            taxonFdFiles[i]=new ArrayList<>();
        }
        int p2Index;
        for (int i = 0; i < chrFdFiles.size(); i++) {
            p2=chrFdFiles.get(i).getName().substring(14,21);
            p2Index=Collections.binarySearch(p2List, p2);
            taxonFdFiles[p2Index].add(chrFdFiles.get(i));
        }
        String[] outNames= IntStream.range(0, p2List.size()).boxed().map(e->
                                   chrFdFiles.get(0).getName().substring(0,13)+p2List.get(e)+".txt").toArray(String[]::new);
        for (int i = 0; i < taxonFdFiles.length; i++) {
            findNearestIBSWindow(genotypeGrid, p3List, taxonMap, taxonFdFiles[i], new File(outDir, outNames[i]));
        }
    }

    /**
     * 在每个染色体内部，处理每个p2种质对应的所有p3种质
     * @param chrGenotypeGrid
     * @param p3List
     * @param taxonMap
     * @param p3FdFilesPerP2
     */
    private static void findNearestIBSWindow(GenotypeGrid chrGenotypeGrid, List<String> p3List,
                                             Map<String, String> taxonMap, List<File> p3FdFilesPerP2,
                                             File outFile){
        long start= System.nanoTime();
        Collections.sort(p3FdFilesPerP2);
        List<ChrRange> chrRanges=extractWindow(p3FdFilesPerP2.get(0));
        String p2 =p3FdFilesPerP2.get(0).getName().substring(14,21);
        int p2Index, p3Index, startPosVCF, endPosVCF, startIndex, endIndex;
        short startChrID, endChrID;
        float ibs;
        String p3;
        p2Index= chrGenotypeGrid.getTaxonIndex(taxonMap.get(p2));
        List<RowTableTool<String>> p3Tables=getP3TablesPerP2(p3FdFilesPerP2);
        List<String> temp;
        TIntArrayList maxFdIndexList;
        double maxFd;
        double[] p3FdArray;
        TDoubleArrayList ibsListOfMaxFd;
        MaxFdMiniIBSP3[] maxFdMiniIBSP3Array=new MaxFdMiniIBSP3[chrRanges.size()];
        for (int i = 0; i < maxFdMiniIBSP3Array.length; i++) {
            maxFdMiniIBSP3Array[i]= new MaxFdMiniIBSP3(Double.MIN_VALUE, Double.MIN_VALUE, -1, new ArrayList<>());
        }
        List<String> p3s;
        for (int i = 0; i < chrRanges.size(); i++) {
            startPosVCF=chrRanges.get(i).getVCFStart();
            endPosVCF=chrRanges.get(i).getVCFEnd()-1;
            startChrID=chrRanges.get(i).getVCFStartChrID();
            endChrID=chrRanges.get(i).getVCFEndChrID();
            startIndex=chrGenotypeGrid.getSiteIndex(startChrID, startPosVCF);
            endIndex=chrGenotypeGrid.getSiteIndex(endChrID, endPosVCF);
            if (startIndex < 0 || endIndex < 0) continue;
            p3FdArray=new double[p3List.size()];
            Arrays.fill(p3FdArray,-1);
            for (int j = 0; j < p3List.size(); j++) {
                temp=p3Tables.get(j).getRow(i);
                if (ifFdChangeTo0(temp.get(8), temp.get(9))) continue;
                p3FdArray[j]=Double.parseDouble(temp.get(9));
            }
            maxFd= Arrays.stream(p3FdArray).max().getAsDouble();
            if (maxFd <= 0) continue;
            maxFdIndexList=new TIntArrayList();
            for (int j = 0; j < p3FdArray.length; j++) {
                if (p3FdArray[j]==maxFd){
                    maxFdIndexList.add(j);
                }
            }
            ibsListOfMaxFd=new TDoubleArrayList();
            for (int j = 0; j < maxFdIndexList.size(); j++) {
                p3=taxonMap.get(p3List.get(maxFdIndexList.get(j)));
                p3Index=chrGenotypeGrid.getTaxonIndex(p3);
                ibs=chrGenotypeGrid.getIBSDistance(p2Index, p3Index, startIndex, endIndex);
                ibsListOfMaxFd.add(ibs);
            }
            double miniIBS=ibsListOfMaxFd.min();
            p3s = new ArrayList<>();
            for (int j = 0; j < ibsListOfMaxFd.size(); j++) {
                if (ibsListOfMaxFd.get(j) != miniIBS) continue;
                p3s.add(p3List.get(maxFdIndexList.get(j)));
            }
            int snpCount= chrGenotypeGrid.getSubGenotypeTableBySite(IntStream.range(startIndex, endIndex+1).toArray()).getSiteNumber();
            maxFdMiniIBSP3Array[i].set(maxFd, miniIBS, snpCount, p3s);
        }
        StringBuilder sb=new StringBuilder();
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chr\tStart\tEnd\tMaxFd\tMiniIBS\tSNPCount\tP3");
            bw.newLine();
            for (int i = 0; i < chrRanges.size(); i++) {
                sb.setLength(0);
                sb.append(chrRanges.get(i).getChr()).append("\t").append(chrRanges.get(i).getStart()).append("\t");
                sb.append(chrRanges.get(i).getEnd()-1).append("\t");
                sb.append(maxFdMiniIBSP3Array[i].toString());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("****** "+p2+" finished in "+String.format("%.4f", Benchmark.getTimeSpanSeconds(start)) +
                " seconds ******");
    }

    public static boolean ifFdChangeTo0(String D, String fd){
        if (ifDChangeTo0(D)) return true;
        if (!StringTool.isNumeric(fd)) return true;
        if (Double.parseDouble(D) <= 0) return true;
        if (Double.parseDouble(fd) < 0) return true;
        if (Double.parseDouble(fd) > 1) return true;
        return false;
    }

    public static boolean ifDChangeTo0(String D){
        if (!StringTool.isNumeric(D)) return true;
        return false;
    }

    public static class MaxFdMiniIBSP3{
        double fd;
        double ibs;
        int snpCount;
        List<String> p3List;

        public MaxFdMiniIBSP3(double fd, double ibs, int snpCount, List<String> p3List){
            this.fd=fd;
            this.ibs=ibs;
            this.snpCount=snpCount;
            this.p3List=p3List;
        }

        public void set(double fd, double ibs, int snpCount, List<String> p3List){
            this.fd=fd;
            this.ibs=ibs;
            this.snpCount=snpCount;
            this.p3List=p3List;
        }

        @Override
        public String toString() {
            NumberFormat numberFormat =NumberFormat.getInstance();
            numberFormat.setGroupingUsed(false);
            numberFormat.setMaximumFractionDigits(5);
            StringBuilder sb = new StringBuilder();
           if (this.ibs==Double.MIN_VALUE || this.fd==Double.MIN_VALUE){
               sb.append("NA").append("\t").append("NA").append("\t").append("NA").append("\t").append("NA");
           }else {
               sb.append(fd).append("\t").append(numberFormat.format(ibs)).append("\t");
               sb.append(snpCount).append("\t");
               sb.append(String.join(",", p3List));
           }
           return sb.toString();
        }
    }

    /**
     * merge all chromosome of an individual
     * @param fdChrDir
     * @param outDir
     */
    public static void merge(String fdChrDir, String outDir){
        System.out.println(DateTime.getDateTimeOfNow());
        List<File> files=IOUtils.getVisibleFileListInDir(fdChrDir);
        Map<String, List<File>> taxaMap=
                files.stream().collect(Collectors.groupingBy(file -> file.getName().substring(13,20)));
        String taxaName;
        List<File> taxaChrFiles;
        RowTableTool<String> table0, table;
        for(Map.Entry<String,List<File>> entry : taxaMap.entrySet()){
            taxaName=entry.getKey();
            taxaChrFiles=entry.getValue();
            table0=new RowTableTool<>(taxaChrFiles.get(0).getAbsolutePath());
            for (int i = 1; i < taxaChrFiles.size(); i++) {
                table=new RowTableTool<>(taxaChrFiles.get(i).getAbsolutePath());
                table0.add(table);
            }
            table0.write(new File(outDir, taxaName+"_vmap2.1.csv.gz"), IOFileFormat.TextGzip);
        }
        System.out.println(DateTime.getDateTimeOfNow());
    }

    /**
     * 提取window
     * @param fdres
     * @return
     */
    public static List<ChrRange> extractWindow(File fdres){
        List<ChrRange> chrRanges=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(fdres)) {
            String line;
            List<String> temp;
            ChrRange chrRange;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line, ",");
                chrRange=new ChrRange(temp.get(0), Integer.parseInt(temp.get(1)), Integer.parseInt(temp.get(2))+1);
                chrRanges.add(chrRange);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return chrRanges;
    }

    private static List<RowTableTool<String>> getP3TablesPerP2(List<File> p3FdFilesPerP2){
        List<RowTableTool<String>> p3tables=new ArrayList<>();
        for (int i = 0; i < p3FdFilesPerP2.size(); i++) {
            p3tables.add(new RowTableTool<>(p3FdFilesPerP2.get(i).getAbsolutePath(),","));
        }
        return p3tables;
    }

    public static void start() {
        String vmapIIVCFDir="/Users/xudaxing/Desktop/fdTest/001_vmap2.1_0.001vcf";
        String taxaInfoFile="/Users/xudaxing/Documents/deleteriousMutation/002_vmapII_taxaGroup/taxa_InfoDB.txt";
        String fdResDir="/Users/xudaxing/Desktop/fdTest/003_fdRes";
        String outDir="/Users/xudaxing/Desktop/fdTest/004_fdOut";
        int index=0;
        NearestIBS.calculateNearestFd(vmapIIVCFDir, taxaInfoFile, fdResDir, outDir, index);
    }
}
