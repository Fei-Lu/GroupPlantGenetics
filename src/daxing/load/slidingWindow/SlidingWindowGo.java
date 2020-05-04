package daxing.load.slidingWindow;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import daxing.common.Triad;
import daxing.load.neutralSiteLoad.DynamicSNPGenotypeDB;
import daxing.load.neutralSiteLoad.IndividualChrLoad;
import daxing.load.neutralSiteLoad.SNPGenotype;
import pgl.infra.table.RowTable;
import pgl.infra.utils.*;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class SlidingWindowGo {

    public static void getDerivedCount(String vcfDir, String ancestralDir, String pgfFile, String triadGeneNameFile,
                                       String outdir, String vmapIIGroupFile){
        RowTableTool<String> table=new RowTableTool<>(triadGeneNameFile);
        List<String> geneNames=table.getColumn(2);
        GeneWindowDB geneWindowDB =new GeneWindowDB(geneNames, pgfFile);
        List<File> vcfFiles= IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> ancestralFiles=IOUtils.getVisibleFileListInDir(ancestralDir);
        Map<String,File> vmapIITaxonoutDirMap=getTaxonOutDirMap(vmapIIGroupFile, outdir);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(vcfFiles.size(), 10);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList= Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e-> getDerivedCount(vcfFiles.get(e), ancestralFiles.get(e),
                    geneWindowDB.getGeneWindow(e+1), vmapIITaxonoutDirMap));
        }
    }

    private static void getDerivedCount(File vcfFile, File ancestralFile, List<GeneWindowDB.GeneWindow> geneWindows,
                                        Map<String,File> taxonOutDirMap){
        long start=System.nanoTime();
        Table<Integer,Integer, Character> ancestralTable=getAncestral(ancestralFile);
        String[] geneName=new String[geneWindows.size()];
        for (int i = 0; i < geneName.length; i++) {
            geneName[i]=geneWindows.get(i).getGeneName();
        }
        try (BufferedReader br = IOTool.getReader(vcfFile)) {
            String line;
            List<String> temp, taxonNames;
            while ((line=br.readLine()).startsWith("##")){}
            temp=PStringUtils.fastSplit(line);
            taxonNames=temp.subList(9, temp.size());
            IndividualChrLoad[] individualChrLoads=new IndividualChrLoad[taxonNames.size()];
            for (int i = 0; i < individualChrLoads.length; i++) {
                individualChrLoads[i]=new IndividualChrLoad(taxonNames.get(i), geneName);
            }
            int chr = -1, pos;
            char refBase, altBase, ancestralBase;
            DynamicSNPGenotypeDB dynamicSNPGenotypeDB=new DynamicSNPGenotypeDB();
            SNPGenotype snpGenotype=null;
            int[][] derivedCount;
            String outFileName=null;
            for (int i = 0; i < geneWindows.size(); i++) {
                while ((line=br.readLine())!=null){
                    temp= PStringUtils.fastSplit(line);
                    chr=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    if (!ancestralTable.contains(chr, pos)) continue;
                    refBase=temp.get(3).charAt(0);
                    altBase=temp.get(4).charAt(0);
                    ancestralBase=ancestralTable.get(chr, pos);
                    if ((ancestralBase!=refBase) && (ancestralBase!=altBase)) continue;
                    if (pos < geneWindows.get(i).getStartOfWindowCDSRange()) continue;
                    if (geneWindows.get(i).containCDSPos(chr, pos)){
                        dynamicSNPGenotypeDB.addSNPGenotype(SNPGenotype.getSNPGenotype(line, ancestralBase, temp.subList(9, temp.size())));
                    }else {
                        snpGenotype=SNPGenotype.getSNPGenotype(line, ancestralBase, temp.subList(9, temp.size()));
                        break;
                    }
                }
                if (dynamicSNPGenotypeDB.size()==0) continue;
                derivedCount=dynamicSNPGenotypeDB.countAllTaxonDerived();
                for (int j = 0; j < derivedCount.length; j++) {
                    individualChrLoads[j].addGeneDerivedCount(i, derivedCount[j]);
                }
                if (i<=geneWindows.size()-2){
                    dynamicSNPGenotypeDB.retainAll(geneWindows.get(i+1));
                    if (geneWindows.get(i+1).containCDSPos(snpGenotype.getChromosome(), snpGenotype.getPosition())){
                        dynamicSNPGenotypeDB.addSNPGenotype(snpGenotype);
                    }
                }
            }
            for (int j = 0; j < individualChrLoads.length; j++) {
                outFileName="chr"+PStringUtils.getNDigitNumber(3, chr)+"."+taxonNames.get(j)+".txt.gz";
                individualChrLoads[j].write(taxonOutDirMap.get(taxonNames.get(j)), outFileName);
            }
            System.out.println("chr"+PStringUtils.getNDigitNumber(3, chr)+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Table<Integer, Integer, Character> getAncestral(File ancestralFile){
        Table<Integer, Integer, Character> ancestralTable= HashBasedTable.create();
        try (BufferedReader br = IOTool.getReader(ancestralFile)) {
            br.readLine();
            String line;
            List<String> temp;
            while ((line=br.readLine())!=null){
                temp= PStringUtils.fastSplit(line);
                ancestralTable.put(Integer.parseInt(temp.get(0)), Integer.parseInt(temp.get(1)), temp.get(2).charAt(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return ancestralTable;
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

    public static void merge(String inputDir, String triadFile, String outDir){
        List<File> dirs=IOTool.getVisibleDir(inputDir);
        int[][] indices= PArrayUtils.getSubsetsIndicesBySubsetSize(dirs.size(), 21);
        for (int i = 0; i < indices.length; i++) {
            Integer[] subLibIndices = new Integer[indices[i][1]-indices[i][0]];
            for (int j = 0; j < subLibIndices.length; j++) {
                subLibIndices[j] = indices[i][0]+j;
            }
            List<Integer> integerList=Arrays.asList(subLibIndices);
            integerList.parallelStream().forEach(e->merge(dirs.get(e), triadFile, outDir));
        }
    }

    private static void merge(File inputDir, String triadFile, String outDir){
        Triad triad=new Triad(triadFile);
        List<File> files= IOUtils.getVisibleFileListInDir(inputDir.getAbsolutePath());
        RowTableTool<String> tableTool=new RowTableTool<>(files.get(0).getAbsolutePath());
        String outName= PStringUtils.fastSplit(files.get(0).getName(), ".").get(1);
        for (int j = 1; j < files.size(); j++) {
            tableTool.add(new RowTableTool<>(files.get(j).getAbsolutePath()));
        }
        List<String> triadIDList=new ArrayList<>();
        String geneName, triadID;
        for (int j = 0; j < tableTool.getRowNumber(); j++) {
            geneName=tableTool.getCell(j, 0);
            triadID=triad.getTraidID(geneName);
            triadIDList.add(triadID);
        }
        tableTool.insertColumn("TriadID", 0, triadIDList);
        Comparator<List<String>> c=Comparator.comparing(l->l.get(0));
        Comparator<List<String>> cc=c.thenComparing(l->l.get(1).substring(8,9));
        tableTool.sortBy(cc);
        tableTool.write(new File(outDir, outName+".txt.gz"), IOFileFormat.TextGzip);
    }
}
