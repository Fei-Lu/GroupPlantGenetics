package daxing.individualIntrogression;

import daxing.common.*;
import gnu.trove.list.array.TIntArrayList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.utils.Benchmark;
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

public class IntrogressionByIndividual {

    /**
     * 并行处理每个染色体
     * @param vmap2VCFDir
     * @param taxa_InfoDB
     * @param fdResDir
     */
    public static void calculateNearestFd(String vmap2VCFDir, String taxa_InfoDB, String fdResDir, String fdOutDir){
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
        Map<String, String> taxonMap= taxonTable.getHashMap(23,0);
        GenotypeGrid genotypeGridA, genotypeGridB, genotypeGrid;
        for (int i = 0; i < chrFdABFiles.length; i++) {
            System.out.println("----------- Start calculate: "+chrAB.get(i)+" -----------");
            genotypeGridA=new GenotypeGrid(abFiles.get(2*i).getAbsolutePath(),GenoIOFormat.VCF_GZ);
            genotypeGridB=new GenotypeGrid(abFiles.get(2*i+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
            genotypeGrid=GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
            genotypeGrid.sortByTaxa();
            calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdABFiles[i], fdOutDir);
            System.out.println("----------- finished: "+chrAB.get(i)+" -----------");
        }
        for (int i = 0; i < chrFdDFiles.length; i++) {
            System.out.println("----------- Start calculate: "+chrD.get(i)+" -----------");
            genotypeGridA=new GenotypeGrid(dFiles.get(2*i).getAbsolutePath(),GenoIOFormat.VCF_GZ);
            genotypeGridB=new GenotypeGrid(dFiles.get(2*i+1).getAbsolutePath(),GenoIOFormat.VCF_GZ);
            genotypeGrid=GenotypeOperation.mergeGenotypesBySite(genotypeGridA,genotypeGridB);
            genotypeGrid.sortByTaxa();
            calculateNearestFdCByTaxon(genotypeGrid, taxonMap, chrFdDFiles[i],fdOutDir);
            System.out.println("----------- finished: "+chrD.get(i)+" -----------");
        }
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
            p2=chrFdFiles.get(i).getName().substring(14,19);
            p3=chrFdFiles.get(i).getName().substring(20,25);
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
            p2=chrFdFiles.get(i).getName().substring(14,19);
            p2Index=Collections.binarySearch(p2List, p2);
            taxonFdFiles[p2Index].add(chrFdFiles.get(i));
        }
        String[] outNames= IntStream.range(0, p2List.size()).boxed().map(e->
                                   chrFdFiles.get(0).getName().substring(0,14)+p2List.get(e)+".csv").toArray(String[]::new);
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
        String p2 =p3FdFilesPerP2.get(0).getName().substring(14,19);
        int p2Index, p3Index, startPosVCF, endPosVCF, startIndex, endIndex;
        short startChrID, endChrID;
        float ibsMini, ibs;
        int[] windowP3Index=new int[chrRanges.size()];
        Arrays.fill(windowP3Index, -1);
        String p3;
        p2Index= chrGenotypeGrid.getTaxonIndex(taxonMap.get(p2));
        float[] ibsArray;
        List<RowTableTool<String>> p3Tables=getP3TablesPerP2(p3FdFilesPerP2);
        List<String> temp;
        for (int i = 0; i < chrRanges.size(); i++) {
            startPosVCF=chrRanges.get(i).getVCFStart();
            endPosVCF=chrRanges.get(i).getVCFEnd();
            startChrID=chrRanges.get(i).getVCFStartChrID();
            endChrID=chrRanges.get(i).getVCFEndChrID();
            startIndex=chrGenotypeGrid.getSiteIndex(startChrID, startPosVCF);
            endIndex=chrGenotypeGrid.getSiteIndex(endChrID, endPosVCF);
            if (startIndex < 0 || endIndex < 0) continue;
            ibsMini=2;
            ibsArray=new float[p3List.size()];
            Arrays.fill(ibsArray, ibsMini);
            int miniIBSP3Index=-1;
            for (int j = 0; j < p3List.size(); j++) {
                temp=p3Tables.get(j).getRow(i);
                if (!StringTool.isNumeric(temp.get(8))) continue;
                if (!StringTool.isNumeric(temp.get(9))) continue;
                if (Double.parseDouble(temp.get(8))<0)continue; // D < 0 continue
                if (Double.parseDouble(temp.get(8))>1)continue; // D > 1 continue
                if (Double.parseDouble(temp.get(9))<0)continue; // fd < 0 continue
                if (Double.parseDouble(temp.get(9))>1)continue; // fd > 1 continue
                p3=taxonMap.get(p3List.get(j));
                p3Index=chrGenotypeGrid.getTaxonIndex(p3);
                ibs=chrGenotypeGrid.getIBSDistance(p2Index, p3Index, startIndex, endIndex);
                ibsArray[j]=ibs;
                if (ibs < ibsMini){
                    miniIBSP3Index=j;
                    ibsMini=ibs;
                }
            }
            windowP3Index[i]=miniIBSP3Index;
        }
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("chr,start,end,mid,sites,sitesUsed,ABBA,BABA,D,fd,fdM,MiniIBSP3");
            bw.newLine();
            List<String> line;
            for (int i = 0; i < chrRanges.size(); i++) {
                p3Index=windowP3Index[i];
                if (p3Index < 0) continue;
                line=p3Tables.get(p3Index).getRow(i);
                line.add(p3List.get(p3Index));
                bw.write(String.join(",",line));
                bw.write("\t");
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("****** "+p2+" finished in "+String.format("%.4f", Benchmark.getTimeSpanMinutes(start)) +
                " minutes ******");
    }

    /**
     * 提取window
     * @param fdres
     * @return
     */
    private static List<ChrRange> extractWindow(File fdres){
        List<ChrRange> chrRanges=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(fdres)) {
            String line;
            List<String> temp;
            ChrRange chrRange;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line, ",");
                chrRange=new ChrRange(temp.get(0), Integer.parseInt(temp.get(1)), Integer.parseInt(temp.get(2)));
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
}
