package daxing.load.neutralSiteLoad;

import com.google.common.collect.HashBasedTable;
import com.google.common.collect.Table;
import daxing.common.IOTool;
import daxing.common.RowTableTool;
import pgl.infra.table.RowTable;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.IOUtils;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;

public class Go {

    public static void getDerivedCount(String vcfDir, String ancestralDir, String pgfFile, String triadGeneNameFile,
                                       String outdir, String vmapIIGroupFile){
        RowTableTool<String> table=new RowTableTool<>(triadGeneNameFile);
        List<String> geneNames=table.getColumn(2);
        GenesDB genesDB =new GenesDB(geneNames, pgfFile);
        List<File> vcfFiles= IOUtils.getVisibleFileListInDir(vcfDir);
        List<File> ancestralFiles=IOUtils.getVisibleFileListInDir(ancestralDir);
        Map<String,File> vmapIITaxonoutDirMap=getTaxonOutDirMap(vmapIIGroupFile, outdir);
        IntStream.range(0, vcfFiles.size()).forEach(e->getDerivedCount(vcfFiles.get(e), ancestralFiles.get(e),
                genesDB.getChrGeneRange(e+1), vmapIITaxonoutDirMap));
    }

    public static void getDerivedCount(File vcfFile, File ancestralFile, GenesDB.GeneRange[] geneRanges
            , Map<String,File> taxonOutDirMap){
        long start=System.nanoTime();
        Table<Integer,Integer, Character> ancestralTable=getAncestral(ancestralFile);
        String[] geneName=new String[geneRanges.length];
        for (int i = 0; i < geneRanges.length; i++) {
            geneName[i]=geneRanges[i].geneName;
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
            int[] derivedCount;
            String outFileName=null;
            for (int i = 0; i < geneRanges.length; i++) {
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    chr=Integer.parseInt(temp.get(0));
                    pos=Integer.parseInt(temp.get(1));
                    if (!ancestralTable.contains(chr, pos)) continue;
                    refBase=temp.get(3).charAt(0);
                    altBase=temp.get(4).charAt(0);
                    ancestralBase=ancestralTable.get(chr, pos);
                    if ((ancestralBase!=refBase) && (ancestralBase!=altBase)) continue;
                    if (pos < geneRanges[i].start) continue;
                    if (pos < geneRanges[i].end){
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
                if (i<=geneRanges.length-2){
                    dynamicSNPGenotypeDB.retainAll(geneRanges[i+1]);
                    if (geneRanges[i+1].isContain(snpGenotype.getChromosome(), snpGenotype.getPosition())){
                        dynamicSNPGenotypeDB.addSNPGenotype(snpGenotype);
                    }
                }
            }
            for (int j = 0; j < individualChrLoads.length; j++) {
                outFileName="chr"+PStringUtils.getNDigitNumber(3, chr)+"."+taxonNames.get(j)+".txt";
                individualChrLoads[j].write(taxonOutDirMap.get(taxonNames.get(j)), outFileName);
            }
            System.out.println(outFileName+" completed in "+ Benchmark.getTimeSpanMinutes(start)+" minutes");
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
}
