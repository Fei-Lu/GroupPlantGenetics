package daxing;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import daxing.common.ChrRange;
import daxing.common.IOTool;
import gnu.trove.list.array.TIntArrayList;
import it.unimi.dsi.fastutil.chars.CharArrayList;
import it.unimi.dsi.fastutil.chars.CharList;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeGrid;
import pgl.infra.dna.genot.GenotypeOperation;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.PStringUtils;
import pgl.infra.utils.wheat.RefV1Utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.IntStream;

public class Start {

    public static void main(String[] args) {

    }

    public static void checkWithOldVmap2(String vmap2_OldDir, String vmap2_newDir, String taxaInfo_1062,
                                         String taxaInfo_643,
                                         String outDir){
        List<File> oldFiles = IOTool.getFileListInDirEndsWith(vmap2_OldDir, "vcf.gz");
        List<File> newFiles = IOTool.getFileListInDirEndsWith(vmap2_newDir, "vcf.gz");
        String[] outNames = newFiles.stream().map(File::getName).map(s -> s.replaceAll("vcf.gz","txt")).toArray(String[]::new);
        GenotypeGrid genotypeGrid_old, genotypeGrid_new, genotypeGridMerged;
        GenotypeTable genotypeTable_old_subset, genotypeTable_new_subset;
        Set<BiSNP> snpSet_old, snpSet_new;
        List<BiSNP> snpList_old, snpList_new, snpList_intersection;
        int[] posIndices_old, posIndices_new;
        short chrID;
        int taxonIndex1, taxonIndex2;
        List<String> vcfIDList;
        Multimap<String, String> taxaVcfIDMap= getTaxonVcfIDMap_1062(taxaInfo_1062);
        Map<String, String> taxonTaxonMap=getTaxaTaxaMap_643(taxaInfo_643);
        String[] taxa_643;
        double ibsDistance;
        String taxon;
        TIntArrayList posList;
        CharList charList;
        try {
            BufferedWriter bw;
            StringBuilder sb =new StringBuilder();
            TIntArrayList heteroSiteIndexes_old, heteroSiteIndexes_new;
            for (int i = 0; i < oldFiles.size(); i++) {
                bw = IOTool.getWriter(new File(outDir, outNames[i]));
                bw.write("VcfID\tChrID\tIbsDistance");
                bw.newLine();
                chrID = Short.parseShort(oldFiles.get(i).getName().substring(3,6));
                snpSet_old = getSNPSet(oldFiles.get(i).getAbsolutePath());
                snpSet_new = getSNPSet(newFiles.get(i).getAbsolutePath());
                snpList_old = new ArrayList<>(snpSet_old);
                snpList_new = new ArrayList<>(snpSet_new);
                Collections.sort(snpList_old);
                Collections.sort(snpList_new);
                posList = new TIntArrayList();
                charList = new CharArrayList();
                for (int j = 0; j < snpList_old.size(); j++) {
                    posList.add(snpList_old.get(j).getPosition());
                    charList.add(snpList_old.get(j).getAlternativeAlleleBase());
                }
                snpList_intersection = new ArrayList<>();
                for (int j = 0; j < snpList_new.size(); j++) {
                    int index = posList.binarySearch(snpList_new.get(j).getPosition());
                    if (index < 0) continue;
                    if (!(charList.getChar(index)==snpList_new.get(j).getAlternativeAlleleBase())){
                        System.out.println(snpList_new.get(j).getAlternativeAlleleBase());
                        System.out.println(charList.getChar(index));
                    }
                    snpList_intersection.add(snpList_new.get(j));
                }
                System.out.println("Intersection SNPs counts is "+snpList_intersection.size());
                Collections.sort(snpList_intersection);
                genotypeGrid_old = new GenotypeGrid(oldFiles.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                genotypeGrid_new = new GenotypeGrid(newFiles.get(i).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                posIndices_old = new int[snpList_intersection.size()];
                posIndices_new = new int[snpList_intersection.size()];
                for (int j = 0; j < snpList_intersection.size(); j++) {
                    posIndices_old[j] = genotypeGrid_old.getSiteIndex(chrID, snpList_intersection.get(j).getPosition());
                    posIndices_new[j] = genotypeGrid_new.getSiteIndex(chrID, snpList_intersection.get(j).getPosition());
                }
                genotypeTable_old_subset = genotypeGrid_old.getSubGenotypeTableBySite(posIndices_old);
                genotypeTable_new_subset = genotypeGrid_new.getSubGenotypeTableBySite(posIndices_new);
                genotypeTable_old_subset.sortBySite();
                genotypeTable_new_subset.sortBySite();
                genotypeGridMerged = GenotypeOperation.mergeGenotypesByTaxon((GenotypeGrid) genotypeTable_old_subset,
                        (GenotypeGrid) genotypeTable_new_subset);
                genotypeGridMerged.sortByTaxa();
                taxa_643 = genotypeGrid_old.getTaxaNames();
                for (String taxonName : taxa_643){
                    if (!taxonTaxonMap.containsKey(taxonName)) continue;
                    taxon = taxonTaxonMap.get(taxonName);
                    vcfIDList = new ArrayList<>(taxaVcfIDMap.get(taxon));
                    taxonIndex1 = genotypeGridMerged.getTaxonIndex(taxonName);
                    for (String vcfID: vcfIDList){
                        taxonIndex2 = genotypeGridMerged.getTaxonIndex(vcfID);
                        ibsDistance = genotypeGridMerged.getIBSDistance(taxonIndex1, taxonIndex2);
                        sb.setLength(0);
                        sb.append(vcfID).append("\t").append(chrID).append("\t").append(ibsDistance);
                        bw.write(sb.toString());
                        bw.newLine();
                    }
                }
                bw.flush();
                bw.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Set<String> getPosSet(String vcfFile){
        Set<String> posSet= new HashSet<>();
        try (BufferedReader br = IOTool.getReader(vcfFile)) {
            String line;
            List<String> temp;
            while ((line=br.readLine()).startsWith("##")) continue;
            while ((line = br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                posSet.add(temp.get(1));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return posSet;
    }

    private static Set<BiSNP> getSNPSet(String vcfFile){
        Set<BiSNP> posSet= new HashSet<>();
        int count=0;
        try (BufferedReader br = IOTool.getReader(vcfFile)) {
            String line;
            char ref, alt;;
            List<String> temp;
            BiSNP biSNP;
            short chr;
            int pos;
            while ((line=br.readLine()).startsWith("##")) continue;
            while ((line = br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                chr = Short.parseShort(temp.get(0));
                pos = Integer.parseInt(temp.get(1));
                ref = temp.get(3).charAt(0);
                alt = temp.get(4).charAt(0);
                biSNP = new BiSNP(chr, pos, ref, alt, "");
                posSet.add(biSNP);
                count++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("Total "+count+" SNPs in "+vcfFile);
        return posSet;
    }

    private static Map<String, String> getTaxaTaxaMap_643(String taxaInfo_643){
        Map<String, String> taxaTaxaMap = new HashMap<>();
        try (BufferedReader br = IOTool.getReader(taxaInfo_643)) {
            String line, taxa;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (temp.get(0).equals("PI466959")) continue;
                if (temp.get(0).equals("CS-2017")){
                    taxa = "CS_8X";
                    taxaTaxaMap.put(temp.get(0), taxa);
                    continue;
                }
                if (temp.get(0).equals("CS-2018")){
                    taxa = "CS_3X";
                    taxaTaxaMap.put(temp.get(0), taxa);
                    continue;
                }
                taxaTaxaMap.put(temp.get(0), temp.get(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return taxaTaxaMap;
    }

    private static Multimap<String, String> getTaxonVcfIDMap_1062(String taxaInfo_1062){
        Multimap<String, String> taxaVcfIDMap= HashMultimap.create();
        try (BufferedReader br = IOTool.getReader(taxaInfo_1062)) {
            String line;
            List<String> temp;
            br.readLine();
            while ((line=br.readLine())!=null){
                temp = PStringUtils.fastSplit(line);
                if (temp.get(28).equals("0")) continue;
                taxaVcfIDMap.put(temp.get(2), temp.get(0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return taxaVcfIDMap;
    }

    public static void getMissingRateHeterozygosityByTaxon(String vmap2Dir, String csID, String outFile){
        List<File> inputFiles = IOTool.getFileListInDirEndsWith(vmap2Dir, "vcf.gz");
        BufferedWriter bw = IOTool.getWriter(outFile);
        try {
            bw.write("Taxa\tChrID\tIBSDistanceTOCS\tMissingRate\tHeterozygosity");
            bw.newLine();
            GenotypeGrid genotypeGrid;
            String[] taxonNames;
            int chrID, taxonIndex, missingNumberByTaxon, nonMissingNumberByTaxon, csIndex;
            double missingRate, heterozygosity, ibsDistanceToCS;;
            StringBuilder sb = new StringBuilder();
            for (int e = 0; e < inputFiles.size(); e++) {
                genotypeGrid= new GenotypeGrid(inputFiles.get(e).getAbsolutePath(), GenoIOFormat.VCF_GZ);
                taxonNames=genotypeGrid.getTaxaNames();
                chrID = Integer.parseInt(inputFiles.get(e).getName().substring(3,6));
                csIndex = genotypeGrid.getTaxonIndex(csID);
                for (String taxon : taxonNames){
                    taxonIndex = genotypeGrid.getTaxonIndex(taxon);
                    missingNumberByTaxon = genotypeGrid.getMissingNumberByTaxon(taxonIndex);
                    nonMissingNumberByTaxon = genotypeGrid.getNonMissingNumberByTaxon(taxonIndex);
                    missingRate = ((double)missingNumberByTaxon)/(missingNumberByTaxon+nonMissingNumberByTaxon);
                    heterozygosity = genotypeGrid.getHeterozygousProportionByTaxon(taxonIndex);
                    ibsDistanceToCS = genotypeGrid.getIBSDistance(taxonIndex, csIndex);
                    sb.setLength(0);
                    sb.append(taxon).append("\t");
                    sb.append(chrID).append("\t");
                    sb.append(ibsDistanceToCS).append("\t");
                    sb.append(missingRate).append("\t").append(heterozygosity);
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

    public static void recombination(String inputFile1, String inputFile2, String outFile){
        try (BufferedReader br1 = IOTool.getReader(inputFile1);
             BufferedReader br2 = IOTool.getReader(inputFile2);
             BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chromosome\tPosition(bp)\tRate(cM/Mb)\tMap(cM)");
            bw.newLine();
            String line1, line2;
            List<String> temp1, temp2;
            br1.readLine();
            br2.readLine();
            String refChr;
            int refPos, start, end;
            int chrID, pos;
            double geneticsPos, recombinationRate;
            List<ChrRange> chrRangeList = new ArrayList<>();
            DoubleList recombinationRateList = new DoubleArrayList();
            ChrRange chrRange;
            StringBuilder sb = new StringBuilder();
            while ((line2= br2.readLine())!=null){
                temp2 = PStringUtils.fastSplit(line2);
                refChr = temp2.get(0).substring(3,5);
                start = Integer.parseInt(temp2.get(1));
                end = Integer.parseInt(temp2.get(2));
                recombinationRate = Double.parseDouble(temp2.get(4));
                chrRange = new ChrRange(refChr, start, end+1);
                chrRangeList.add(chrRange);
                recombinationRateList.add(recombinationRate);
            }
            while((line1= br1.readLine())!=null){
                temp1 = PStringUtils.fastSplit(line1);
                refChr = temp1.get(1).substring(3,5);
                refPos = Integer.parseInt(temp1.get(2));
                geneticsPos= Double.parseDouble(temp1.get(3));
                chrRange = new ChrRange(refChr, refPos, refPos+1);
                int hit = Collections.binarySearch(chrRangeList, chrRange);
                int index =  hit < 0 ? -hit-2 : hit;
                recombinationRate = recombinationRateList.getDouble(index);
                chrID = RefV1Utils.getChrID(refChr, refPos);
                pos = RefV1Utils.getPosOnChrID(refChr, refPos);
                sb.setLength(0);
                sb.append(chrID).append("\t").append(pos).append("\t").append(recombinationRate).append("\t");
                sb.append(geneticsPos);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static void biAllele(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, "vcf.gz");
        String[] outNames=
                files.stream().map(File::getName).map(s -> s.replaceAll("2.0.vcf.gz","2.0.vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line, subLine;
                List<String> temp;
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write("#CHROM\t"+line.substring(5));
                bw.newLine();
                int cnt=0;
                Set<String> altSet=new HashSet<>();
                while ((line=br.readLine())!=null){
                    subLine=line.substring(0,40);
                    temp= PStringUtils.fastSplit(subLine);
                    if (temp.get(4).length() > 1){
                        altSet.add(temp.get(4));
                    }
                    if (temp.get(4).equals("A") || temp.get(4).equals("T") ||
                            temp.get(4).equals("C") || temp.get(4).equals("G") ||
                            temp.get(4).equals("I") || temp.get(4).equals("D")){
                        bw.write(line);
                        bw.newLine();
                        cnt++;
                    }else {
                        altSet.add(temp.get(4));
                    }
                }
                bw.flush();
//                altSet.stream().forEach(System.out::println);
                System.out.println(new File(outDir, outNames[e]).getName()+  " have "+cnt+" biVariants");
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }

    public static void correctMaf(String vmap2Dir, String outDir){
        List<File> files =IOTool.getFileListInDirEndsWith(vmap2Dir, ".vcf.gz");
        String[] outNames=
                files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz", ".vcf")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw =IOTool.getWriter(new File(outDir, outNames[e]))) {
                String line;
                List<String> temp,tem,te;
                double maf;
                StringBuilder sb = new StringBuilder();
                while ((line=br.readLine()).startsWith("##")){
                    bw.write(line);
                    bw.newLine();
                }
                bw.write(line);
                bw.newLine();
                while ((line=br.readLine())!=null){
                    temp=PStringUtils.fastSplit(line);
                    tem=PStringUtils.fastSplit(temp.get(7), ";");
                    te=PStringUtils.fastSplit(tem.get(6), "=");
                    maf = Double.parseDouble(te.get(1));
                    maf = maf > 0.5 ? 1-maf : maf;
                    sb.setLength(0);
                    sb.append(String.join("\t", temp.subList(0, 7))).append("\t");
                    sb.append(String.join(";", tem.subList(0, 6))).append(";");
                    sb.append(te.get(0)).append("=").append(maf).append("\t");
                    sb.append(String.join("\t", temp.subList(8, temp.size())));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
            } catch (IOException ioException) {
                ioException.printStackTrace();
            }
        });
    }


    public static void extract(String inputDir, String outDir){
        List<File> files = IOTool.getFileListInDirEndsWith(inputDir, ".vcf.gz");
        String[] outNames = files.stream().map(File::getName).map(s -> s.replaceAll(".vcf.gz",".txt.gz")).toArray(String[]::new);
        IntStream.range(0, files.size()).parallel().forEach(e->{
            String line, subLine, ref, alt;
            List<String> temp, tem, te, t;
            String[] acgt={"A","C","G","T"};
            int[] countACGT;
            double maf;
            String chr;
            int refIndex, altIndex, minorCount, majorCount, chrID, pos, refPos, majorIndex, minorIndex;
            int refAlleleCount, altAlleleCount;
            try (BufferedReader br = IOTool.getReader(files.get(e));
                 BufferedWriter bw = IOTool.getWriter(new File(outDir, outNames[e]))) {
                while ((line = br.readLine()).startsWith("##")) continue;
                StringBuilder sb= new StringBuilder();
                while ((line=br.readLine())!=null){
                    subLine= line.substring(0,150);
                    temp=PStringUtils.fastSplit(subLine);
                    ref = temp.get(3);
                    alt = temp.get(4);
                    tem= PStringUtils.fastSplit(temp.get(7), ";");
                    te=PStringUtils.fastSplit(tem.get(6), "=");
                    maf = Double.parseDouble(te.get(1));
                    te = PStringUtils.fastSplit(tem.get(4), "=");
                    t = PStringUtils.fastSplit(te.get(1),  ",");
                    refAlleleCount = Integer.parseInt(t.get(0)) * 2 + Integer.parseInt(t.get(1));
                    altAlleleCount = Integer.parseInt(t.get(2)) * 2 + Integer.parseInt(t.get(1));
                    refIndex = Arrays.binarySearch(acgt, ref);
                    altIndex = Arrays.binarySearch(acgt, alt);
                    minorCount = (int)Math.round(199*maf);
                    majorCount = 199-minorCount;
                    minorIndex = refAlleleCount > altAlleleCount ? altIndex : refIndex;
                    majorIndex = refAlleleCount > altAlleleCount ? refIndex : altIndex;
                    countACGT = new int[acgt.length];
                    countACGT[minorIndex]=minorCount;
                    countACGT[majorIndex]=majorCount;
                    sb.setLength(0);
                    chrID = Integer.parseInt(temp.get(0));
                    pos = Integer.parseInt(temp.get(1));
                    chr= RefV1Utils.getChromosome(chrID, pos);
                    refPos = RefV1Utils.getPosOnChromosome(chrID, pos);
                    sb.append("chr").append(chr).append("\t").append(refPos-1).append("\t").append(refPos).append("\t");
                    for (int i = 0; i < countACGT.length; i++) {
                        sb.append(countACGT[i]).append(",");
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