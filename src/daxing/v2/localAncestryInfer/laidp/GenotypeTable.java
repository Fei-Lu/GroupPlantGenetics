package daxing.v2.localAncestryInfer.laidp;

import com.google.common.primitives.Doubles;
import daxing.common.bisnp.SNP;
import daxing.common.chrrange.ChrPos;
import daxing.common.utiles.ArrayTool;
import daxing.common.utiles.IOTool;
import gnu.trove.iterator.TIntIterator;
import gnu.trove.set.TIntSet;
import gnu.trove.set.hash.TIntHashSet;
import it.unimi.dsi.fastutil.doubles.DoubleArrayList;
import it.unimi.dsi.fastutil.doubles.DoubleList;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.text.NumberFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.stream.IntStream;

/**
 * 记得改BitSet loop
 * instance of this class is for haploid
 */
public class GenotypeTable {

    String[] taxa = null;
    SNP[] snps = null;

    /**
     * Bit genotype by site, the first dimension is site; the second dimension is genotype(haploid), and missing
     */
    BitSet[][] genoSite;

    /**
     * Bit genotype by taxon, the first dimension is taxon; the second dimension is genotype(haploid), and missing
     */
    BitSet[][] genoTaxon;

    private static final int BLOCK_SIZE = 4000; // for parallel calculation, such as alternative allele frequency

    public GenotypeTable(String haploidGenotypeFile){
        try {
            BufferedReader br = IOTool.getReader(haploidGenotypeFile);
            String temp;
            while ((temp = br.readLine()).startsWith("##")) {}
            List<String> l;
            l = PStringUtils.fastSplit(temp);
            this.taxa = new String[l.size()-9];
            for (int i = 9; i < l.size(); i++) {
                taxa[i-9] = l.get(i);
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<BlockLines>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            List<String> lines = new ArrayList<>();
            StringBuilder sb = new StringBuilder();
            NumberFormat numberFormat= NumberFormat.getNumberInstance();
            numberFormat.setGroupingUsed(true);
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()%BlockLines.blockSize == 0) {
                    BlockLines gsb = new BlockLines(lines, startIndex);
                    Future<BlockLines> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=lines.size();
                    lines = new ArrayList<>();
                }
                siteCount++;
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(haploidGenotypeFile);
                    System.out.println(sb);
                }
            }
            br.close();
            if (lines.size() != 0) {
                BlockLines gsb = new BlockLines(lines, startIndex);
                Future<BlockLines> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new SNP[siteCount];
            for (Future<BlockLines> blockVCFFuture : resultList) {
                BlockLines block = blockVCFFuture.get();
                for (int j = 0; j < block.getActBlockSize(); j++) {
                    this.snps[block.getStartIndex() + j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex() + j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(numberFormat.format(this.getSiteNumber())).append(" SNPs are in ").append(haploidGenotypeFile).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
    }

    private GenotypeTable(String[] taxa, SNP[] snps, BitSet[][] genoSite, BitSet[][] genoTaxon){
        this.taxa=taxa;
        this.snps=snps;
        this.genoSite=genoSite;
        this.genoTaxon=genoTaxon;
    }

    public static GenotypeTable newInstanceFromGenoTaxon(String[] taxa, SNP[] snps, BitSet[][] genoTaxon){
       GenotypeTable genotypeTable = new GenotypeTable(taxa, snps, null, genoTaxon);
       genotypeTable.transposeTaxonToSite();
       return genotypeTable;
    }

    public SNP[] getSnps() {
        return snps;
    }

    public BitSet[][] getGenoSite() {
        return genoSite;
    }

    public BitSet[][] getGenoTaxon() {
        return genoTaxon;
    }

    public String[] getTaxa() {
        return taxa;
    }

    public int getTaxaNumber() {
        return this.taxa.length;
    }

    public int getSiteNumber() {
        return this.snps.length;
    }

    public int getPosition(int siteIndex) {
        return this.snps[siteIndex].getPos();
    }

    public int getTaxonIndex(String taxon) {
        return ArrayUtils.indexOf(this.taxa, taxon);
//        return Arrays.binarySearch(taxa, taxon);
    }
    
    public int[] getTaxaIndices(List<String> taxaList){
        int[] res = new int[taxaList.size()];
        for (int i = 0; i < taxaList.size(); i++) {
            res[i] = this.getTaxonIndex(taxaList.get(i));
        }
        return res;
    }

    public int[][] getTaxaIndices(List<String>[] taxaList){
        int[][] res = new int[taxaList.length][];
        for (int i = 0; i < taxaList.length; i++) {
            res[i] = this.getTaxaIndices(taxaList[i]);
        }
        return res;
    }

    public int getSiteIndex(String chromosome, int position) {
        ChrPos query = new ChrPos(chromosome, position);
        return Arrays.binarySearch(this.snps, query);
    }

    private void transposeSiteToTaxon () {
        long start = System.nanoTime();
        genoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < genoTaxon.length; i++) {
            genoTaxon[i]=new BitSet[2];
            for (int j = 0; j < genoTaxon[i].length; j++) {
                genoTaxon[i][j]=new BitSet(this.getSiteNumber());
            }
        }
        for (int i = 0; i < genoSite.length; i++) {
            for (int j = 0; j < genoSite[i].length; j++) {
                for (int k = 0; k < genoSite[i][j].size(); k++) {
                    if (genoSite[i][j].get(k)){
                        genoTaxon[k][j].set(i);
                    }
                }
            }
        }
        System.out.println("Transpose genoSite to genoTaxon takes " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
    }

    public void transposeTaxonToSite () {
        long start = System.nanoTime();
        genoSite = new BitSet[this.getSiteNumber()][];
        for (int i = 0; i < genoSite.length; i++) {
            genoSite[i]=new BitSet[2];
            for (int j = 0; j <genoSite[i].length; j++) {
                genoSite[i][j]=new BitSet(this.getTaxaNumber());
            }
        }
        for (int i = 0; i < genoTaxon.length; i++) {
            for (int j = 0; j < genoTaxon[i].length; j++) {
                for (int k = 0; k < genoTaxon[i][j].size(); k++) {
                    if (genoTaxon[i][j].get(k)){
                        genoSite[k][j].set(i);
                    }
                }
            }
        }
        System.out.println("Transpose genoTaxon to genoSite takes " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
    }

    public boolean isMissing(int siteIndex, int taxonIndex) {
        return this.genoSite[siteIndex][1].get(taxonIndex);
    }

    public boolean isAlternativeAllele(int siteIndex, int taxonIndex){
        return this.genoSite[siteIndex][0].get(taxonIndex);
    }

    public byte getAlleleByte(int siteIndex, int taxonIndex){
        if(this.isMissing(siteIndex, taxonIndex)) return AlleleEncoder.genotypeMissingByte;
        byte referenceAlleleByte = this.getSnps()[siteIndex].getReferenceAlleleByte();
        byte alternativeAlleleByte = this.getSnps()[siteIndex].getAlternativeAlleleByte();
        return isAlternativeAllele(siteIndex, taxonIndex) ? alternativeAlleleByte : referenceAlleleByte;
    }

    public char getAlleleBase(int siteIndex, int taxonIndex){
        byte alleleByte = this.getAlleleByte(siteIndex, taxonIndex);
        return AlleleEncoder.getAlleleBaseFromByte(alleleByte);
    }

    /**
     *
     * @param threadsNum threadsNum
     * @param taxaIndices dim1 is different populations, dim2 is different taxon
     * @return alternative allele frequency, dim1 is different populations, dim2 is variants.
     * -1 means sample in one population is all missing
     */
    public double[][] calculateAltAlleleFrequency(int threadsNum, int[][] taxaIndices) {
        BitSet[] pop_bitSet = new BitSet[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            pop_bitSet[i] = new BitSet();
            for (int j = 0; j < taxaIndices[i].length; j++) {
                pop_bitSet[i].set(taxaIndices[i][j]);
            }
        }
        int numVariants = this.getSiteNumber();
        double[][] alts = new double[taxaIndices.length][numVariants];
        int numBlock = (numVariants + BLOCK_SIZE - 1) / BLOCK_SIZE;
        ExecutorService executorService = Executors.newFixedThreadPool(threadsNum);
        List<Future<Void>> futures = new ArrayList<>();
        for (int blockIndex = 0; blockIndex < numBlock; blockIndex++) {
            int startIndex = blockIndex * BLOCK_SIZE;
            int endIndex = Math.min((blockIndex + 1) * BLOCK_SIZE, numVariants);
            futures.add(executorService.submit(() -> {
                int count1, missing, total;
                BitSet[][] genoSite = this.getGenoSite();
                BitSet[] bitSets, clonedBitSets;
                int bitSetsLen;
                for (int variantsIndex = startIndex; variantsIndex < endIndex; variantsIndex++) {
                    bitSets = genoSite[variantsIndex];
                    bitSetsLen = genoSite[variantsIndex].length;
                    for (int i = 0; i < taxaIndices.length; i++) {
                        clonedBitSets = new BitSet[bitSetsLen];
                        for (int j = 0; j < bitSets.length; j++) {
                            clonedBitSets[j] = (BitSet) bitSets[j].clone();
                        }
                        clonedBitSets[0].and(pop_bitSet[i]);
                        clonedBitSets[1].and(pop_bitSet[i]);
                        count1 = clonedBitSets[0].cardinality();
                        missing =  clonedBitSets[1].cardinality();
                        total = pop_bitSet[i].cardinality();
                        if (missing == total){
                            // all sample in these pop is missing
                            alts[i][variantsIndex] = -1;
                        }else {
                            alts[i][variantsIndex] = (double) count1 / (total - missing);
                        }
                    }
                }
                return null;
            }));
        }
        for (Future<Void> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                executorService.shutdown();
                throw new RuntimeException(e);
            }
        }
        executorService.shutdown();
        return alts;
    }

    /**
     *
     * @param threadsNum threadsNum
     * @param taxaIndices dim1 is different populations, dim2 is different taxon
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @return derived allele frequency, dim1 is different populations, dim2 is variants.
     * -1 mean missing.
     * Missing include two cases:
     * case1: sample in one population is all missing
     * case2: ancestral allele is missing
     */
    public double[][] calculateDaf(int threadsNum, int[][] taxaIndices, BitSet[] ancestralAlleleBitSet){
        double[][] alts = this.calculateAltAlleleFrequency(threadsNum, taxaIndices);
        double[][] dafs = new double[alts.length][alts[0].length];
        for (int i = 0; i < alts.length; i++) {
            System.arraycopy(alts[i], 0, dafs[i], 0, alts[i].length);
        }
        for (int i = ancestralAlleleBitSet[0].nextSetBit(0); i >= 0; i = ancestralAlleleBitSet[0].nextSetBit(i+1)) {
            for (int j = 0; j < dafs.length; j++) {
                // missing
                if (alts[j][i] < 0){
                    dafs[j][i] = -1;
                }else {
                    dafs[j][i] = 1 - alts[j][i];
                }
            }
        }

        // site which ancestral allele is missing will be filling with -1
        for (int i = ancestralAlleleBitSet[1].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[1].nextSetBit(i+1)) {
            for (int j = 0; j < dafs.length; j++) {
                dafs[j][i] = -1;
            }
        }
        return dafs;
    }

    /**
     *
     * @param taxonIndex taxa belong to outGroup population
     * @return ancestralAlleleBitSet
     * dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     */
    public BitSet[] getAncestralAlleleFromTaxa(int[] taxonIndex){
        BitSet[] bitSets = new BitSet[2];
        BitSet ancestralAllele = (BitSet) this.genoTaxon[taxonIndex[0]][0].clone();
        for (int i = 1; i < taxonIndex.length; i++) {
            ancestralAllele.and(this.genoTaxon[taxonIndex[i]][0]);
        }
        BitSet missing = (BitSet) this.genoTaxon[taxonIndex[0]][1].clone();
        for (int i = 1; i < taxonIndex.length; i++) {
            missing.or(this.genoTaxon[taxonIndex[i]][1]);
        }
        bitSets[0] = ancestralAllele;
        bitSets[1] = missing;
        return bitSets;
    }

    /**
     * @param windowSize window size of sliding window, the unit is variants, not bp
     * @param stepSize step size of sliding window, the unit is variants, not bp
     * @return windowStartIndex array, site start index of all windows
     */
    public static int[] getWindowStartIndex(int windowSize, int stepSize, int numVariants){
        int numberWindow = (numVariants + stepSize - 1) / stepSize;
        IntList windowList = new IntArrayList();
//        int[] window = new int[numberWindow];
        int windowStartIndex;
        for (int windowIndex = 0; windowIndex < numberWindow; windowIndex++) {
            windowStartIndex = Math.min(windowIndex * stepSize, numVariants);
//            window[windowIndex] = windowStartIndex;
            if (windowStartIndex + windowSize > numVariants){
                windowList.add(windowStartIndex);
                break;
            }
            windowList.add(windowStartIndex);
        }
        return windowList.toIntArray();
    }

    /**
     *
     * @param pop_taxonIndex dim1 is different populations, dim2 is different taxa
     * @param windowStartIndexArray site start index of all windows
     * @param windowSize window size of sliding window, the unit is variants, not bp
     * @return pairwise dxy between all populations in all windows
     * dim1 is different window, same index as windowStartIndexArray,
     * dim2 is pairwise dxy between all populations
     */
    public double[][] calculateDxy(int[][] pop_taxonIndex, int[] windowStartIndexArray, int windowSize) {
        int pop_num = pop_taxonIndex.length;
        int variantsNum = this.getSiteNumber();
//        SNP snp_start, snp_end;
//        int[] windowLen = new int[windowStartIndexArray.length];
//        for (int i = 0; i < windowLen.length; i++) {
//            snp_start = this.snps[windowStartIndexArray[i]];
//            snp_end = this.snps[Math.min(windowStartIndexArray[i]+windowSize, this.getSiteNumber())-1];
//            windowLen[i] = snp_end.getPos() - snp_start.getPos() + 1;
//        }
        double[][] dxyArrays = new double[windowStartIndexArray.length][pop_num*(pop_num-1)/2];
        int hapCountA, hapCountB, hapCountAB;
        double[] cumDifference;
        BitSet bitSet, missing, subBitset;
        int k =0;
        int[] variantsNum_window;
        for (int popIndexA = 0; popIndexA < pop_num -1; popIndexA++) {
            for (int popIndexB = popIndexA+1; popIndexB < pop_num; popIndexB++) {
                hapCountA = pop_taxonIndex[popIndexA].length;
                hapCountB = pop_taxonIndex[popIndexB].length;
                cumDifference = new double[windowStartIndexArray.length];
                variantsNum_window = new int[windowStartIndexArray.length];
                for (int taxonA : pop_taxonIndex[popIndexA]){
                    for (int taxonB : pop_taxonIndex[popIndexB]){
                        bitSet = (BitSet) this.genoTaxon[taxonA][0].clone();
                        missing = (BitSet)this.genoTaxon[taxonA][1].clone();
                        missing.or(this.genoTaxon[taxonB][1]);
                        bitSet.xor(this.genoTaxon[taxonB][0]);
                        bitSet.andNot(missing);
                        for (int i = 0; i < windowStartIndexArray.length; i++) {
                            subBitset = bitSet.get(windowStartIndexArray[i], Math.min(windowStartIndexArray[i]+windowSize, variantsNum));
                            cumDifference[i] += subBitset.cardinality();
                            variantsNum_window[i] = Math.min(windowStartIndexArray[i]+windowSize, variantsNum) - windowStartIndexArray[i];
                        }
                    }
                }
                hapCountAB = hapCountA*hapCountB;
                for (int i = 0; i < cumDifference.length; i++) {
//                    dxyArrays[i][k] = cumDifference[i]/hapCountAB/windowLen[i];
                    dxyArrays[i][k] = cumDifference[i]/hapCountAB/variantsNum_window[i];
                }
                k++;
            }
        }
        return dxyArrays;
    }

    public double[][] calculatePairwiseDxy(int[] pop_taxonIndex_native,
                                           int[][] pop_taxonIndex_introgressed, int[] windowStartIndexArray, int windowSize) {
        int[][] pop_taxonIndex_native_introgressed = new int[1 + pop_taxonIndex_introgressed.length][];
        pop_taxonIndex_native_introgressed[0] = new int[pop_taxonIndex_native.length];
        System.arraycopy(pop_taxonIndex_native, 0, pop_taxonIndex_native_introgressed[0], 0, pop_taxonIndex_native.length);
        for (int i = 0; i < pop_taxonIndex_introgressed.length; i++) {
            pop_taxonIndex_native_introgressed[i+1] = new int[pop_taxonIndex_introgressed[i].length];
            System.arraycopy(pop_taxonIndex_introgressed[i], 0, pop_taxonIndex_native_introgressed[i + 1], 0,
                    pop_taxonIndex_introgressed[i].length);
        }
        return this.calculateDxy(pop_taxonIndex_native_introgressed, windowStartIndexArray, windowSize);
    }

    /**
     *
     * @param pop_taxonIndex_admixed taxon index of admixed population
     * @param pop_taxonIndex_native  taxon index of native population
     * @param pop_taxonIndex_introgressed taxon index of all introgressed population
     * @param windowStartIndexArray site start index of all windows
     * @param windowSize window size of sliding window, the unit is variants, not bp
     * @return dxy between all admixed taxon and all source (native + introgressed) populations
     * dim1 is different window, same index as windowStartIndexArray,
     * dim2 is taxon index of admixed population
     * dim3 is population index of native_introgressed populations
     */
    public double[][][] calculateAdmixedDxy(int[] pop_taxonIndex_admixed, int[] pop_taxonIndex_native,
                                             int[][] pop_taxonIndex_introgressed,
                                             int[] windowStartIndexArray, int windowSize){
        int[][] pop_taxonIndex_native_introgressed = new int[1 + pop_taxonIndex_introgressed.length][];
        pop_taxonIndex_native_introgressed[0] = new int[pop_taxonIndex_native.length];
        System.arraycopy(pop_taxonIndex_native, 0, pop_taxonIndex_native_introgressed[0], 0, pop_taxonIndex_native.length);
        for (int i = 0; i < pop_taxonIndex_introgressed.length; i++) {
            pop_taxonIndex_native_introgressed[i+1] = new int[pop_taxonIndex_introgressed[i].length];
            System.arraycopy(pop_taxonIndex_introgressed[i], 0, pop_taxonIndex_native_introgressed[i + 1], 0,
                    pop_taxonIndex_introgressed[i].length);
        }
        int variantsNum = this.getSiteNumber();
        double[][][] dxyArrays = new double[windowStartIndexArray.length][pop_taxonIndex_admixed.length][pop_taxonIndex_native_introgressed.length];
        SNP snp_start, snp_end;
        int[] windowLen = new int[windowStartIndexArray.length];
        for (int i = 0; i < windowLen.length; i++) {
            snp_start = this.snps[windowStartIndexArray[i]];
            snp_end = this.snps[Math.min(windowStartIndexArray[i]+windowSize, this.getSiteNumber())-1];
            windowLen[i] = snp_end.getPos() - snp_start.getPos() + 1;
        }
        int hapCount_admixed, hapCount_native_introgressed, hapCount;
        double[] cumDifference;
        BitSet bitSet, missing, subBitset;
        int popNum = pop_taxonIndex_native_introgressed.length;
        for (int popTaxonIndex = 0; popTaxonIndex < popNum; popTaxonIndex++) {
            for (int admixedTaxonIndex = 0; admixedTaxonIndex < pop_taxonIndex_admixed.length; admixedTaxonIndex++) {
                hapCount_native_introgressed = pop_taxonIndex_native_introgressed[popTaxonIndex].length;
                hapCount_admixed = 1;
                cumDifference = new double[windowStartIndexArray.length];
                for (int taxon: pop_taxonIndex_native_introgressed[popTaxonIndex]){
                    bitSet = (BitSet) this.genoTaxon[taxon][0].clone();
                    missing = (BitSet) this.genoTaxon[taxon][1].clone();
                    missing.or(this.genoTaxon[admixedTaxonIndex][1]);
                    bitSet.xor(this.genoTaxon[admixedTaxonIndex][0]);
                    bitSet.andNot(missing);
                    for (int i = 0; i < windowStartIndexArray.length; i++) {
                        subBitset = bitSet.get(windowStartIndexArray[i], Math.min(windowStartIndexArray[i]+windowSize, variantsNum));
                        cumDifference[i] += subBitset.cardinality();
                    }
                }
                hapCount = hapCount_native_introgressed*hapCount_admixed;
                for (int i = 0; i < cumDifference.length; i++) {
                    dxyArrays[i][admixedTaxonIndex][popTaxonIndex] = cumDifference[i]/hapCount/windowLen[i];
                }
            }
        }
        return dxyArrays;
    }

    /**
     *
     * @param threadsNum threadsNum
     * @param pop_taxonIndex dim1 is different populations, dim2 is different taxa
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @param ifZsocre if calculate_pattersonD_f zscore of pattersonD_f
     * @return pattersonD_f, with or without zscore
     */
    public double[] calculatePattersonD_f(int threadsNum, int[][] pop_taxonIndex, BitSet[] ancestralAlleleBitSet,
                                          boolean ifZsocre){
        double[][] dafs = this.calculateDaf(threadsNum, pop_taxonIndex, ancestralAlleleBitSet);
        int[] p3_taxaIndices = pop_taxonIndex[2];
        int p3_taxaIndices_len = p3_taxaIndices.length;
        int[][] p3_ab_taxaIndices = new int[2][];
        for (int i = 0; i < p3_ab_taxaIndices.length; i++) {
            p3_ab_taxaIndices[i] = new int[p3_taxaIndices_len/2];
        }
        System.arraycopy(p3_taxaIndices, 0, p3_ab_taxaIndices[0], 0, p3_taxaIndices_len/2);
        System.arraycopy(p3_taxaIndices, p3_taxaIndices_len/2, p3_ab_taxaIndices[1], 0, p3_taxaIndices_len/2);
        double[][] dafs_p3_ab = this.calculateDaf(threadsNum, p3_ab_taxaIndices, ancestralAlleleBitSet);
        double[] d_f = GenotypeTable.calculate_pattersonD_f(dafs, dafs_p3_ab, Integer.MAX_VALUE, Integer.MAX_VALUE);
        if (!ifZsocre) return d_f;
        int totalVariants = dafs[0].length;
        // default， split to 20
        double[] d_f_zscore = GenotypeTable.getJackknife_pattersonD_f_zscore(d_f, dafs, dafs_p3_ab, totalVariants/20, totalVariants);
        double[] res = new double[d_f.length+d_f_zscore.length];
        System.arraycopy(d_f, 0, res, 0, d_f.length);
        System.arraycopy(d_f_zscore, 0, res, d_f.length, d_f_zscore.length);
        return res;
    }

    /**
     * @param d_f pattersonD and f
     * @param dafs derived allele frequency, dim1 is different populations, dim2 is variants.
     *             length of dim1 is 3
     * @param dafs_p3_ab dim1 is p3a and p3b, dim2 is taxa indices
     * @param block_size linkage disequilibrium decays to background levels when using block_size variants
     * @param totalVariants total variants in genotype
     * @return pattersonD_f_zscore
     */
    public static double[] getJackknife_pattersonD_f_zscore(double[] d_f, double[][] dafs, double[][] dafs_p3_ab,
                                                            int block_size, int totalVariants){
        int[] windowStartIndex = GenotypeTable.getWindowStartIndex(block_size/2,block_size, totalVariants);
        int[] randoms = ArrayTool.getRandomNonrepetitionArray(windowStartIndex.length, 0, windowStartIndex.length);
        double[][] jackknife_D_f = new double[2][randoms.length];
        for (int i = 0; i < randoms.length; i++) {
            jackknife_D_f[0][i] = GenotypeTable.calculate_pattersonD_f(dafs, dafs_p3_ab, windowStartIndex[randoms[i]],
                    block_size)[0];
            jackknife_D_f[1][i] = GenotypeTable.calculate_pattersonD_f(dafs, dafs_p3_ab, windowStartIndex[randoms[i]],
                    block_size)[1];
        }
        DescriptiveStatistics stats;
        double standardDeviation;
//        double standError;
        double[] d_f_zscore = new double[4];
        for (int i = 0; i < 2; i++) {
            stats = new DescriptiveStatistics(jackknife_D_f[i]);
            standardDeviation = stats.getStandardDeviation();
//            standError = standardDeviation/Math.sqrt(randoms.length);
//            d_f_zscore[i] = d_f[i]/standError;
            d_f_zscore[i] = d_f[i]/standardDeviation;
            d_f_zscore[i+2] = standardDeviation;
        }
        return d_f_zscore;
    }

    /**
     * @param dafs derived allele frequency, dim1 is different populations, dim2 is variants.
     *             length of dim1 is 3
     * @param randomStartSiteIndex random site index when doing Jackknife
     * @param block_size linkage disequilibrium decays to background levels when using block_size variants
     * @return pattersonD_f
     */
    private static double[] calculate_pattersonD_f(double[][] dafs, double[][] dafs_p3_ab, int randomStartSiteIndex, int block_size){
        DoubleList abbaList = new DoubleArrayList();
        DoubleList babaList = new DoubleArrayList();
        DoubleList abba_p3abList = new DoubleArrayList();
        DoubleList baba_p3abList = new DoubleArrayList();
        for (int i = 0; i < dafs[0].length; i++) {
            if (dafs[0][i] < 0 || dafs[1][i] < 0 || dafs[2][i] < 0) continue;
            if (dafs_p3_ab[0][i] < 0 || dafs_p3_ab[1][i] < 0) continue;

            // for jackknife
            if (i >= randomStartSiteIndex && i < (randomStartSiteIndex+block_size)) continue;

            abbaList.add((1-dafs[0][i]) * dafs[1][i] * dafs[2][i]);
            babaList.add(dafs[0][i] * (1-dafs[1][i]) * dafs[2][i]);
            abba_p3abList.add((1-dafs[0][i]) * dafs_p3_ab[0][i] * dafs_p3_ab[1][i]);
            baba_p3abList.add(dafs[0][i] * (1-dafs_p3_ab[0][i]) * dafs_p3_ab[1][i]);
        }
        double abba=0, baba=0, abba_p3ab=0, baba_p3ab=0;
        for (int i = 0; i < abbaList.size(); i++) {
            abba += abbaList.getDouble(i);
            baba += babaList.getDouble(i);
            abba_p3ab += abba_p3abList.getDouble(i);
            baba_p3ab += baba_p3abList.getDouble(i);
        }
        double[] res = new double[2];
        res[0] = (abba-baba)/(abba+baba);
        res[1] = (abba-baba)/(abba_p3ab-baba_p3ab);
        return res;
    }


    /**
     *
     * @param threadsNum threadsNum
     * @param dafs_native derived allele frequency of native population
     * @param dafs_admixed derived allele frequency of each admixed taxon
     * @param dafs_introgressed derivedd allele frequency of all introgressed populations
     * @param windowStartIndexArray site start index of all windows
     * @param windowSize window size of sliding window, the unit is variants, not bp
     * @param variantsNum total variants in genotype
     * @return fd_windows_admixed_introgressed,
     * dim1 is window index
     * dim2 is admixed taxon
     * dim3 is introgressed population
     */
    public static double[][][] calculate_fd(int threadsNum, double[] dafs_native, double[][] dafs_admixed,
                                            double[][] dafs_introgressed, int[] windowStartIndexArray, int windowSize,
                                            int variantsNum){
        int windowStartIndex, windowEndIndex;
        List<Future<Void>> futures= new ArrayList<>();
        ExecutorService executorService = Executors.newFixedThreadPool(threadsNum);
        double[][][] fd_windows_admixed_introgressed = new double[windowStartIndexArray.length][dafs_admixed.length][dafs_introgressed.length];
        for (int i = 0; i < windowStartIndexArray.length; i++) {
            windowStartIndex = windowStartIndexArray[i];
            windowEndIndex = Math.min(windowStartIndex+windowSize, variantsNum);
            int windowIndexFinal = i;
            int windowStartFinalIndex = windowStartIndex;
            int windowEndFinalIndex = windowEndIndex;
            futures.add(executorService.submit(()->{
                DoubleList[][] abbaList = new DoubleList[dafs_admixed.length][dafs_introgressed.length];
                DoubleList[][] babaList = new DoubleList[dafs_admixed.length][dafs_introgressed.length];
                DoubleList[][] abba_pDList = new DoubleList[dafs_admixed.length][dafs_introgressed.length];
                DoubleList[][] baba_pDList = new DoubleList[dafs_admixed.length][dafs_introgressed.length];
                for (int j = 0; j < dafs_admixed.length; j++) {
                    for (int k = 0; k < dafs_introgressed.length; k++) {
                        abbaList[j][k] = new DoubleArrayList();
                        babaList[j][k] = new DoubleArrayList();
                        abba_pDList[j][k] = new DoubleArrayList();
                        baba_pDList[j][k] = new DoubleArrayList();
                    }
                }
                boolean ifVariantContinue= false;
                double daf_native, daf_admixed, daf_introgressed, daf_pD;
                double abba, baba, abba_pD, baba_pD;
                for (int variantIndex = windowStartFinalIndex; variantIndex < windowEndFinalIndex; variantIndex++) {
                    if (dafs_native[variantIndex] < 0){
                        ifVariantContinue = true;
                    }
                    for (double[] value : dafs_admixed) {
                        if (value[variantIndex] < 0) {
                            ifVariantContinue = true;
                            break;
                        }
                    }
                    for (double[] doubles : dafs_introgressed) {
                        if (doubles[variantIndex] < 0) {
                            ifVariantContinue = true;
                            break;
                        }
                    }
                    if (ifVariantContinue) continue;
                    for (int admixedTaxonIndex = 0; admixedTaxonIndex < dafs_admixed.length; admixedTaxonIndex++) {
                        for (int introgressedPopIndex = 0; introgressedPopIndex < dafs_introgressed.length; introgressedPopIndex++) {
                            daf_native = dafs_native[variantIndex];
                            daf_admixed = dafs_admixed[admixedTaxonIndex][variantIndex];
                            daf_introgressed = dafs_introgressed[introgressedPopIndex][variantIndex];
                            daf_pD = daf_admixed > 0 ? daf_admixed : daf_introgressed;
                            abba = (1-daf_native) * daf_admixed * daf_introgressed;
                            baba = daf_native * (1-daf_admixed) * daf_introgressed;
                            abba_pD = daf_pD * daf_pD* (1-daf_native);
                            baba_pD = daf_pD * (1-daf_pD) * daf_native;
                            abbaList[admixedTaxonIndex][introgressedPopIndex].add(abba);
                            babaList[admixedTaxonIndex][introgressedPopIndex].add(baba);
                            abba_pDList[admixedTaxonIndex][introgressedPopIndex].add(abba_pD);
                            baba_pDList[admixedTaxonIndex][introgressedPopIndex].add(baba_pD);
                        }
                    }
                }
                double[][] abbaSum = new double[dafs_admixed.length][dafs_introgressed.length];
                double[][] babaSum = new double[dafs_admixed.length][dafs_introgressed.length];
                double[][] abba_pDSum = new double[dafs_admixed.length][dafs_introgressed.length];
                double[][] baba_pDSum = new double[dafs_admixed.length][dafs_introgressed.length];
                for (int admixedTaxonIndex = 0; admixedTaxonIndex < abbaList.length; admixedTaxonIndex++) {
                    for (int introgressedPopIndex = 0; introgressedPopIndex < abbaList[admixedTaxonIndex].length; introgressedPopIndex++) {
                        for (int j = 0; j < abbaList[admixedTaxonIndex][introgressedPopIndex].size(); j++) {
                            abbaSum[admixedTaxonIndex][introgressedPopIndex]+=abbaList[admixedTaxonIndex][introgressedPopIndex].getDouble(j);
                            babaSum[admixedTaxonIndex][introgressedPopIndex]+=babaList[admixedTaxonIndex][introgressedPopIndex].getDouble(j);
                            abba_pDSum[admixedTaxonIndex][introgressedPopIndex]+=abba_pDList[admixedTaxonIndex][introgressedPopIndex].getDouble(j);
                            baba_pDSum[admixedTaxonIndex][introgressedPopIndex]+=baba_pDList[admixedTaxonIndex][introgressedPopIndex].getDouble(j);
                        }
                    }
                }
                double[][] pattersonD = new double[dafs_admixed.length][dafs_introgressed.length];
                double[][] fd = new double[dafs_admixed.length][dafs_introgressed.length];
                double difference_ABBA_BABA, sum_ABBA_BABA, difference_ABBA_BABA_pD;
                for (int admixedTaxonIndex = 0; admixedTaxonIndex < abbaSum.length; admixedTaxonIndex++) {
                    for (int introgressedPopIndex = 0; introgressedPopIndex < abbaSum[admixedTaxonIndex].length; introgressedPopIndex++) {
                        difference_ABBA_BABA = abbaSum[admixedTaxonIndex][introgressedPopIndex] - babaSum[admixedTaxonIndex][introgressedPopIndex];
                        sum_ABBA_BABA = abbaSum[admixedTaxonIndex][introgressedPopIndex] + babaSum[admixedTaxonIndex][introgressedPopIndex];
                        difference_ABBA_BABA_pD = abba_pDSum[admixedTaxonIndex][introgressedPopIndex] - baba_pDSum[admixedTaxonIndex][introgressedPopIndex];
                        if (sum_ABBA_BABA == 0 || difference_ABBA_BABA_pD == 0){
                            pattersonD[admixedTaxonIndex][introgressedPopIndex] = 0;
                            fd[admixedTaxonIndex][introgressedPopIndex] = 0;
                        }else {
                            pattersonD[admixedTaxonIndex][introgressedPopIndex] = difference_ABBA_BABA/sum_ABBA_BABA;
                            fd[admixedTaxonIndex][introgressedPopIndex] = difference_ABBA_BABA/difference_ABBA_BABA_pD;
                        }

                        if (pattersonD[admixedTaxonIndex][introgressedPopIndex] <=0 || fd[admixedTaxonIndex][introgressedPopIndex] < 0){
                            fd[admixedTaxonIndex][introgressedPopIndex] = 0;
                        }
                        if (fd[admixedTaxonIndex][introgressedPopIndex] > 1){
                            fd[admixedTaxonIndex][introgressedPopIndex] = 1;
                        }
                    }
                }
                fd_windows_admixed_introgressed[windowIndexFinal] = fd;
                return null;
            }));
        }

        for (Future<Void> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                executorService.shutdown();
                throw new RuntimeException(e);
            }
        }
        executorService.shutdown();
        try {
            if (!executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.MILLISECONDS)) {
                // task not been completed
                executorService.shutdownNow();
            }
        } catch (InterruptedException e) {
            executorService.shutdownNow();
        }
        return fd_windows_admixed_introgressed;

    }

    /**
     *
     * @param fd_windows_admixed fd_windows_admixed_introgressed
     *                                        dim1 is window index，
     *                                        dim2 is admixed taxon，
     *                                        dim3 is introgressed population
     * @param dxy_pairwise_nativeIntrogressed pairwise dxy between all source population (native + introgressed
     *                                        population)
     *                                        dim1 is different window, same index as windowStartIndexArray,
     *                                        dim2 is pairwise dxy between all populations
     * @param dxy_windows_admixed_nativeIntrogressed dxy between all admixed taxon and all source (native + introgressed) populations
     *                                               dim1 is different window, same index as windowStartIndexArray,
     *                                               dim2 is taxon index of admixed population,
     *                                               dim3 is population index of native_introgressed populations
     * @return grid source, dim1 is admixed taxon index, dim2 is window index
     */
    public static int[][] calculateSource(double[][][] fd_windows_admixed,
                                          double[][] dxy_pairwise_nativeIntrogressed,
                                          double[][][] dxy_windows_admixed_nativeIntrogressed, double[] f_upperLimit){
        int windowNum = fd_windows_admixed.length;
        int admixedTaxaNum = fd_windows_admixed[0].length;
        int introgressedPopNum = fd_windows_admixed[0][0].length;
        int[][] gridSource = new int[admixedTaxaNum][windowNum];
        for (int[] ints : gridSource) {
            Arrays.fill(ints, 1);
        }

        double[][] fd_threshold = new double[admixedTaxaNum][introgressedPopNum];
        double[] fd_taxon_introgressedPop;
        for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaNum; admixedTaxonIndex++) {
            for (int introgressedPopIndex = 0; introgressedPopIndex < introgressedPopNum; introgressedPopIndex++) {
                fd_taxon_introgressedPop = new double[windowNum];
                for (int i = 0; i < windowNum; i++) {
                    fd_taxon_introgressedPop[i] = fd_windows_admixed[i][admixedTaxonIndex][introgressedPopIndex];
                }
                fd_threshold[admixedTaxonIndex][introgressedPopIndex] = StatUtils.percentile(fd_taxon_introgressedPop
                        , (1-f_upperLimit[introgressedPopIndex])*100);
            }
        }
        double dxy_min, dxy_p1p2, dxy_p2p3, fd, fd_thresh;
        Source source;
        for (int windowIndex = 0; windowIndex < windowNum; windowIndex++) {
            dxy_min = Double.MAX_VALUE;
            for (int i = 0; i < dxy_pairwise_nativeIntrogressed[windowIndex].length; i++) {
                dxy_min = Math.min(dxy_pairwise_nativeIntrogressed[windowIndex][i], dxy_min);
            }
            for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaNum; admixedTaxonIndex++) {

                dxy_p1p2 = dxy_windows_admixed_nativeIntrogressed[windowIndex][admixedTaxonIndex][0];

                // 可能是ILS或introgress population间有基因流
                if (dxy_min < dxy_p1p2) continue;

                for (int introgressedPopIndex = 1; introgressedPopIndex < introgressedPopNum +1; introgressedPopIndex++) {
                    fd = fd_windows_admixed[windowIndex][admixedTaxonIndex][introgressedPopIndex-1];

                    // 保留fd值大于群体混合比例的
                    fd_thresh = fd_threshold[admixedTaxonIndex][introgressedPopIndex-1];
                    if (fd < fd_thresh) continue;

                    dxy_p2p3 = dxy_windows_admixed_nativeIntrogressed[windowIndex][admixedTaxonIndex][introgressedPopIndex];
                    if (dxy_p2p3 < dxy_p1p2){
                        source = Source.getInstanceFromIndex(introgressedPopIndex).get();
                        if (gridSource[admixedTaxonIndex][windowIndex] == 1){
                            gridSource[admixedTaxonIndex][windowIndex] = source.getFeature();
                        }else {
                            gridSource[admixedTaxonIndex][windowIndex] += source.getFeature();
                        }
                    }
                }
            }
        }
        return gridSource;
    }

    /**
     *
     * @param source grid source, dim1 is admixed taxon index, dim2 is window index
     * @param conjunctionNum conjunctionNum
     * @return window start end index of successive window
     * dim1 is taxon
     */
    public static List<int[]>[] getSuccessiveIntrogressionWindow(int[][] source, int conjunctionNum){
        List<int[]>[] successiveWindow_taxon = new List[source.length];
        for (int i = 0; i < successiveWindow_taxon.length; i++) {
            successiveWindow_taxon[i] = new ArrayList<>();
        }
        int windowNum = source[0].length;
        TIntSet introgressionWindowIndexSet;
        for (int taxonIndex = 0; taxonIndex < source.length; taxonIndex++) {
            introgressionWindowIndexSet = new TIntHashSet();
            for (int windowIndex = 0; windowIndex < windowNum; windowIndex++) {
                if (source[taxonIndex][windowIndex] > 1){
                    introgressionWindowIndexSet.add(windowIndex);
                }
            }
            TIntSet resultIndexSet =  new TIntHashSet(introgressionWindowIndexSet);
            TIntIterator intIterator = introgressionWindowIndexSet.iterator();
            while (intIterator.hasNext()){
                int index = intIterator.next();
                int i = 1;
                while (i <= conjunctionNum){
                    resultIndexSet.add(Math.max(index - i, 0));
                    resultIndexSet.add(Math.min(index+i, windowNum-1));
                    i++;
                }
            }

            if (resultIndexSet.size()==0) continue;

            // introgression window index set to sorted arrat
            int[] resultIndexArray = resultIndexSet.toArray();
            Arrays.sort(resultIndexArray);

            // int[2] start end represent successive introgression window index
            List<int[]> indexList = new ArrayList<>();
            int initializeStartIndex=resultIndexArray[0];
            int initializeEndIndex=initializeStartIndex;
            int[] startEnd;
            for (int i = 1; i < resultIndexArray.length; i++) {
                if (resultIndexArray[i]-resultIndexArray[i-1]==1){
                    initializeEndIndex = resultIndexArray[i];
                }else {
                    startEnd = new int[2];
                    startEnd[0] = initializeStartIndex;
                    startEnd[1] = initializeEndIndex;
                    indexList.add(startEnd);
                    initializeStartIndex = resultIndexArray[i];
                    initializeEndIndex = initializeStartIndex ;
                }
            }
            startEnd = new int[2];
            startEnd[0]=initializeStartIndex;
            startEnd[1] = initializeEndIndex;
            indexList.add(startEnd);

            successiveWindow_taxon[taxonIndex] = indexList;

        }
        return successiveWindow_taxon;
    }

    /**
     *
     * @param windowSize The unit of window size is one variant, not base pairs (bp)
     * @param stepSize The unit of step size is one variant, not base pairs (bp)
     * @param taxaGroupFile taxaGroupFile
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @param conjunctionNum conjunctionNum
     * @param switchCostScore switchCostScore
     * @param threadsNum threadsNum
     * @return BitSet[][] local ancestry, dim1 is admixed taxon index, dim2 is n_way admixture, index equal Source.index
     */
    public BitSet[][] calculateLocalAncestry(int windowSize, int stepSize, String taxaGroupFile,
                                             BitSet[] ancestralAlleleBitSet, int conjunctionNum,
                                             double switchCostScore, int maxSolutionCount, int threadsNum){
        int variantsNum = this.getSiteNumber();
        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);
        int n_wayAdmixture = taxaGroup.getIntrogressedPopTaxa().length + 1;

        int[] admixedTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.ADMIXED));
        int[] nativeTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.NATIVE));
        int[][] introgressedPopTaxaIndices = this.getTaxaIndices(taxaGroup.getIntrogressedPopTaxa());
        int[][] native_admixed_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getNative_admixed_introgressed_Taxa());
        int[] native_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getTaxaOf_native_introgressed());

        int[] windowStartIndexArray = GenotypeTable.getWindowStartIndex(windowSize, stepSize, variantsNum);
        double[][][] dxy_windows_admixed = this.calculateAdmixedDxy(admixedTaxaIndices, nativeTaxaIndices,
                introgressedPopTaxaIndices, windowStartIndexArray, windowSize);
        double[][] dxy_pairwise_nativeIntrogressed = this.calculatePairwiseDxy(nativeTaxaIndices,introgressedPopTaxaIndices,
                windowStartIndexArray, windowSize);

        double[][] dafs = this.calculateDaf(threadsNum, native_admixed_introgressed_popIndex, ancestralAlleleBitSet);
        double[] dafs_native = dafs[0];
        double[][] dafs_admixed = new double[taxaGroup.getTaxaOf(Source.ADMIXED).size()][];
        System.arraycopy(dafs, 1, dafs_admixed, 0, taxaGroup.getTaxaOf(Source.ADMIXED).size());
        double[][] dafs_introgressed = new double[taxaGroup.getIntrogressedPopTaxa().length][];
        System.arraycopy(dafs, taxaGroup.getTaxaOf(Source.ADMIXED).size()+1, dafs_introgressed, 0,
                taxaGroup.getIntrogressedPopTaxa().length);

        double[][][] fd = GenotypeTable.calculate_fd(threadsNum, dafs_native, dafs_admixed, dafs_introgressed,
                windowStartIndexArray,
                windowSize, variantsNum);

        double[][] d_f_z_sd = this.get_D_f_z_sd(threadsNum, taxaGroup, ancestralAlleleBitSet, false);
        double[] f_upperLimit = GenotypeTable.get_upperLimit_f(d_f_z_sd);
        int[][] gridSource = GenotypeTable.calculateSource(fd, dxy_pairwise_nativeIntrogressed,
                dxy_windows_admixed, f_upperLimit);

        List<int[]>[] successiveWindow_taxon = GenotypeTable.getSuccessiveIntrogressionWindow(gridSource, conjunctionNum);
        BitSet[] queriesGenotype = this.getTaxaGenotype(admixedTaxaIndices);
        BitSet[] sourcesGenotype = this.getTaxaGenotype(native_introgressed_popIndex);
        List<String> sourceTaxaList = this.getTaxaList(native_introgressed_popIndex);
        Map<String, Source> taxaSourceMap = taxaGroup.getTaxaSourceMap(sourceTaxaList);
        BitSet[][] localAncestry = new BitSet[admixedTaxaIndices.length][];
        for (int i = 0; i < admixedTaxaIndices.length; i++) {
            localAncestry[i] = new BitSet[n_wayAdmixture];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new BitSet();
            }
            // initialize native ancestry to 1
            localAncestry[i][0].set(0, variantsNum);
        }

        ExecutorService executorService = Executors.newFixedThreadPool(threadsNum);
        List<Future<Void>> futures = new ArrayList<>();
        for (int i = 0; i < admixedTaxaIndices.length; i++) {
            for (int j = 0; j < successiveWindow_taxon[i].size(); j++) {
                int admixedTaxonIndex = i;
                int gridIndex = j;
                futures.add(executorService.submit(()->{
                    int gridStart = successiveWindow_taxon[admixedTaxonIndex].get(gridIndex)[0]; // inclusive
                    int gridEnd = successiveWindow_taxon[admixedTaxonIndex].get(gridIndex)[1]; // inclusive
                    int fragmentStartIndex = windowStartIndexArray[gridStart]; // inclusive
                    int fragmentEndIndex = Math.min(windowStartIndexArray[gridEnd]+windowSize, variantsNum); // exclusive
                    BitSet queryFragment = queriesGenotype[admixedTaxonIndex].get(fragmentStartIndex, fragmentEndIndex);
                    BitSet[] sourcesFragment = new BitSet[sourcesGenotype.length];
                    for (int sourceIndex = 0; sourceIndex < sourcesGenotype.length; sourceIndex++) {
                        sourcesFragment[sourceIndex]=sourcesGenotype[sourceIndex].get(fragmentStartIndex,fragmentEndIndex);
                    }
                    int fragmentLen = fragmentEndIndex - fragmentStartIndex;
                    EnumMap<Direction, IntList[]> biDirectionCandidateSolution =
                            Solution.getBiDirectionCandidateSourceSolution(sourcesFragment, queryFragment, fragmentLen,
                                    switchCostScore, sourceTaxaList, taxaSourceMap, maxSolutionCount);
                    IntList solution = Solution.calculateBreakPoint(biDirectionCandidateSolution);
                    EnumSet<Source> sources;
                    int intervalStartIndex, intervalEndIndex;
                    for (int k = 0; k < solution.size(); k=k+3) {
                        sources = Source.getSourcesFrom(solution.getInt(k));
                        intervalStartIndex = solution.getInt(k+1)+fragmentStartIndex; // inclusive
                        intervalEndIndex = solution.getInt(k+2)+fragmentStartIndex; // inclusive
                        for (Source source : sources){
                            if (source.equals(Source.NATIVE)) continue;
                            localAncestry[admixedTaxonIndex][source.getIndex()].set(intervalStartIndex, intervalEndIndex+1);
                            localAncestry[admixedTaxonIndex][0].set(intervalStartIndex, intervalEndIndex+1, false);
                        }
                    }
                    return null;
                }));
            }
        }
        for (Future<Void> future : futures) {
            try {
                future.get();
            } catch (InterruptedException | ExecutionException e) {
                executorService.shutdown();
                throw new RuntimeException(e);
            }
        }
        executorService.shutdown();
        return localAncestry;
    }

    /**
     *
     * @param windowSize The unit of window size is one variant, not base pairs (bp)
     * @param stepSize The unit of step size is one variant, not base pairs (bp)
     * @param taxaGroupFile taxaGroupFile
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @param conjunctionNum conjunctionNum
     * @param switchCostScore switchCostScore
     * @param threadsNum threadsNum
     * @return BitSet[][] local ancestry, dim1 is admixed taxon index, dim2 is n_way admixture, index equal Source.index
     */
    public BitSet[][] calculateLocalAncestry_singleThread(int windowSize, int stepSize, String taxaGroupFile,
                                                          BitSet[] ancestralAlleleBitSet, int conjunctionNum,
                                                          double switchCostScore, int maxSolutionCount, int threadsNum){
        int variantsNum = this.getSiteNumber();
        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);
        int n_wayAdmixture = taxaGroup.getIntrogressedPopTaxa().length + 1;

        int[] admixedTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.ADMIXED));
        int[] nativeTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.NATIVE));
        int[][] introgressedPopTaxaIndices = this.getTaxaIndices(taxaGroup.getIntrogressedPopTaxa());
        int[][] native_admixed_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getNative_admixed_introgressed_Taxa());
        int[] native_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getTaxaOf_native_introgressed());

        int[] windowStartIndexArray = GenotypeTable.getWindowStartIndex(windowSize, stepSize, variantsNum);
        double[][][] dxy_windows_admixed = this.calculateAdmixedDxy(admixedTaxaIndices, nativeTaxaIndices,
                introgressedPopTaxaIndices, windowStartIndexArray, windowSize);
        double[][] dxy_pairwise_nativeIntrogressed = this.calculatePairwiseDxy(nativeTaxaIndices,introgressedPopTaxaIndices,
                windowStartIndexArray, windowSize);

        double[][] dafs = this.calculateDaf(threadsNum, native_admixed_introgressed_popIndex, ancestralAlleleBitSet);
        double[] dafs_native = dafs[0];
        double[][] dafs_admixed = new double[taxaGroup.getTaxaOf(Source.ADMIXED).size()][];
        System.arraycopy(dafs, 1, dafs_admixed, 0, taxaGroup.getTaxaOf(Source.ADMIXED).size());
        double[][] dafs_introgressed = new double[taxaGroup.getIntrogressedPopTaxa().length][];
        System.arraycopy(dafs, taxaGroup.getTaxaOf(Source.ADMIXED).size()+1, dafs_introgressed, 0,
                taxaGroup.getIntrogressedPopTaxa().length);

        double[][][] fd = GenotypeTable.calculate_fd(threadsNum, dafs_native, dafs_admixed, dafs_introgressed,
                windowStartIndexArray,
                windowSize, variantsNum);

        double[][] d_f_z_sd = this.get_D_f_z_sd(threadsNum, taxaGroup, ancestralAlleleBitSet, false);
        double[] f_upperLimit = GenotypeTable.get_upperLimit_f(d_f_z_sd);
        int[][] gridSource = GenotypeTable.calculateSource(fd, dxy_pairwise_nativeIntrogressed,
                dxy_windows_admixed, f_upperLimit);

        List<int[]>[] successiveWindow_taxon = GenotypeTable.getSuccessiveIntrogressionWindow(gridSource, conjunctionNum);
        BitSet[] queriesGenotype = this.getTaxaGenotype(admixedTaxaIndices);
        BitSet[] sourcesGenotype = this.getTaxaGenotype(native_introgressed_popIndex);
        List<String> sourceTaxaList = this.getTaxaList(native_introgressed_popIndex);
        Map<String, Source> taxaSourceMap = taxaGroup.getTaxaSourceMap(sourceTaxaList);
        BitSet[][] localAncestry = new BitSet[admixedTaxaIndices.length][];
        for (int i = 0; i < admixedTaxaIndices.length; i++) {
            localAncestry[i] = new BitSet[n_wayAdmixture];
            for (int j = 0; j < localAncestry[i].length; j++) {
                localAncestry[i][j] = new BitSet();
            }
            // initialize native ancestry to 1
            localAncestry[i][0].set(0, variantsNum);
        }
        int gridStart, gridEnd, fragmentStartIndex, fragmentEndIndex, fragmentLen, intervalStartIndex, intervalEndIndex;
        BitSet queryFragment;
        BitSet[] sourcesFragment;
        EnumMap<Direction, IntList[]> biDirectionCandidateSolution;
        IntList solution;
        EnumSet<Source> sources;
        for (int admixedTaxonIndex = 0; admixedTaxonIndex < admixedTaxaIndices.length; admixedTaxonIndex++) {
            for (int gridIndex = 0; gridIndex < successiveWindow_taxon[admixedTaxonIndex].size(); gridIndex++) {
//                int admixedTaxonIndex = i;
//                int gridIndex = j;
                gridStart = successiveWindow_taxon[admixedTaxonIndex].get(gridIndex)[0]; // inclusive
                gridEnd = successiveWindow_taxon[admixedTaxonIndex].get(gridIndex)[1]; // inclusive
                fragmentStartIndex = windowStartIndexArray[gridStart]; // inclusive
                fragmentEndIndex = Math.min(windowStartIndexArray[gridEnd]+windowSize, variantsNum); // exclusive
                queryFragment = queriesGenotype[admixedTaxonIndex].get(fragmentStartIndex, fragmentEndIndex);
                sourcesFragment = new BitSet[sourcesGenotype.length];
                for (int sourceIndex = 0; sourceIndex < sourcesGenotype.length; sourceIndex++) {
                    sourcesFragment[sourceIndex]=sourcesGenotype[sourceIndex].get(fragmentStartIndex,fragmentEndIndex);
                }
                fragmentLen = fragmentEndIndex - fragmentStartIndex;
                biDirectionCandidateSolution =
                        Solution.getBiDirectionCandidateSourceSolution(sourcesFragment, queryFragment, fragmentLen,
                                switchCostScore, sourceTaxaList, taxaSourceMap, maxSolutionCount);
                solution = Solution.calculateBreakPoint(biDirectionCandidateSolution);
                for (int k = 0; k < solution.size(); k=k+3) {
                    sources = Source.getSourcesFrom(solution.getInt(k));
                    intervalStartIndex = solution.getInt(k+1)+fragmentStartIndex; // inclusive
                    intervalEndIndex = solution.getInt(k+2)+fragmentStartIndex; // inclusive
                    for (Source source : sources){
                        if (source.equals(Source.NATIVE)) continue;
                        localAncestry[admixedTaxonIndex][source.getIndex()].set(intervalStartIndex, intervalEndIndex+1);
                        localAncestry[admixedTaxonIndex][0].set(intervalStartIndex, intervalEndIndex+1, false);
                    }
                }
            }
        }
        return localAncestry;
    }

    public double[][] get_D_f_z_sd(int threadsNum, TaxaGroup taxaGroup, BitSet[] ancestralAlleleBitSet,
                                   boolean ifZsocre){
        int introgressedPopNum = taxaGroup.getIntrogressedPopTaxa().length;
        double[][] d_f_z_sd = new double[introgressedPopNum][6];
        int[][] native_admixed_introgressed_popIndex;
        List<Source> sources = new ArrayList<>();
        sources.add(Source.NATIVE);
        sources.add(Source.ADMIXED);
        StringBuilder sb = new StringBuilder();
        Source source;
        for (int i = 0; i < introgressedPopNum; i++) {
            source = Source.getInstanceFromIndex(i+1).get();
            sources.add(source);
            native_admixed_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getTaxaOf(sources));
            d_f_z_sd[i] = this.calculatePattersonD_f(threadsNum, native_admixed_introgressed_popIndex,
                    ancestralAlleleBitSet, ifZsocre);
            sources.remove(2);
            System.out.println();
            System.out.println(source.name());
            if (ifZsocre){
                System.out.println("PattersonD, f, zscore_PattersonD, zscore_f, sd_PattersonD, sd_f");
            }else {
                System.out.println("PattersonD, f");
            }

            sb.setLength(0);
            sb.append(Doubles.join(", ", d_f_z_sd[i]));
            System.out.println(sb);
            System.out.println();
        }
        return d_f_z_sd;
    }

    public static double[] get_upperLimit_f(double[][] d_f_z_sd){
        int introgressedPopNum = d_f_z_sd.length;
        double[] f_upperLimit = new double[introgressedPopNum];
        for (int i = 0; i < introgressedPopNum; i++) {
            // 2 sd
//            f_upperLimit[i] = d_f_z_sd[i][1]+2*d_f_z_sd[i][5];
            f_upperLimit[i] = d_f_z_sd[i][1];
        }
        return f_upperLimit;
    }

    /**
     *
     * @param windowSize The unit of window size is one variant, not base pairs (bp)
     * @param stepSize The unit of step size is one variant, not base pairs (bp)
     * @param taxaGroupFile taxaGroupFile
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @param threadsNum threadsNum
     * @param fdOutDir fdOutDir
     * @param individualLocalAncestryDir individualLocalAncestryDir
     */
    public void calculateLocalAncestry_check_fd_gridSource(int windowSize, int stepSize, String taxaGroupFile,
                                             BitSet[] ancestralAlleleBitSet, int threadsNum, String fdOutDir,
                                                           String individualLocalAncestryDir){
        int variantsNum = this.getSiteNumber();
        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);

        int[] admixedTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.ADMIXED));
        int[] nativeTaxaIndices = this.getTaxaIndices(taxaGroup.getTaxaOf(Source.NATIVE));
        int[][] introgressedPopTaxaIndices = this.getTaxaIndices(taxaGroup.getIntrogressedPopTaxa());
        int[][] native_admixed_introgressed_popIndex = this.getTaxaIndices(taxaGroup.getNative_admixed_introgressed_Taxa());

        int[] windowStartIndexArray = GenotypeTable.getWindowStartIndex(windowSize, stepSize, variantsNum);
        double[][][] dxy_windows_admixed = this.calculateAdmixedDxy(admixedTaxaIndices, nativeTaxaIndices,
                introgressedPopTaxaIndices, windowStartIndexArray, windowSize);
        double[][] dxy_pairwise_nativeIntrogressed = this.calculatePairwiseDxy(nativeTaxaIndices,introgressedPopTaxaIndices,
                windowStartIndexArray, windowSize);

        double[][] dafs = this.calculateDaf(threadsNum, native_admixed_introgressed_popIndex, ancestralAlleleBitSet);
        double[] dafs_native = dafs[0];
        double[][] dafs_admixed = new double[taxaGroup.getTaxaOf(Source.ADMIXED).size()][];
        System.arraycopy(dafs, 1, dafs_admixed, 0, taxaGroup.getTaxaOf(Source.ADMIXED).size());
        double[][] dafs_introgressed = new double[taxaGroup.getIntrogressedPopTaxa().length][];
        System.arraycopy(dafs, taxaGroup.getTaxaOf(Source.ADMIXED).size()+1, dafs_introgressed, 0,
                taxaGroup.getIntrogressedPopTaxa().length);

        double[][][] fd = GenotypeTable.calculate_fd(threadsNum, dafs_native, dafs_admixed, dafs_introgressed,
                windowStartIndexArray,
                windowSize, variantsNum);
        double[][] d_f_z_sd = this.get_D_f_z_sd(threadsNum, taxaGroup, ancestralAlleleBitSet, true);
        double[] f_upperLimit = GenotypeTable.get_upperLimit_f(d_f_z_sd);
        int[][] gridSource = GenotypeTable.calculateSource(fd, dxy_pairwise_nativeIntrogressed,
                dxy_windows_admixed, f_upperLimit);

        double[][][] fd_byTaxon = new double[fd[0].length][fd.length][fd[0][0].length];
        for (int i = 0; i < fd.length; i++) {
            for (int j = 0; j < fd[i].length; j++) {
                fd_byTaxon[j][i] = fd[i][j];
            }
        }
        String[] outName_fd = IntStream.range(60, 90).boxed().map(e->"tsk_"+e+"_fd.txt").toArray(String[]::new);
        String[] outName_individualLocalAncestry = IntStream.range(60, 90).boxed().map(e->"tsk_"+e+".txt").toArray(String[]::new);

        BufferedWriter bw;
        try {
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < fd_byTaxon.length; i++) {
                bw = IOTool.getWriter(new File(fdOutDir, outName_fd[i]));
                sb.setLength(0);
                sb.append("scaffold\tstart\tend\t");
                for (int j = 0; j < fd_byTaxon[0][0].length; j++) {
                    sb.append(Source.getInstanceFromIndex(j+1).get().name()).append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < fd_byTaxon[i].length; j++) {
                    int windowStartSiteIndex = windowStartIndexArray[j]; // inclusive
                    int windowEndSiteIndex = Math.min(windowStartIndexArray[j]+windowSize, variantsNum); // inclusive
                    sb.setLength(0);
                    sb.append("1").append("\t");
                    sb.append(this.getPosition(windowStartSiteIndex)).append("\t");
                    sb.append(this.getPosition(windowEndSiteIndex-1)).append("\t");
                    sb.append(Doubles.join("\t", fd_byTaxon[i][j]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }

            for (int i = 0; i < gridSource.length; i++) {
                bw = IOTool.getWriter(new File(individualLocalAncestryDir, outName_individualLocalAncestry[i]));
                sb.setLength(0);
                sb.append("scaffold\tstart\tend\tsourceFeature");
                bw.write(sb.toString());
                bw.newLine();
                for (int j = 0; j < gridSource[i].length; j++) {
                    int windowStartSiteIndex = windowStartIndexArray[j]; // inclusive
                    int windowEndSiteIndex = Math.min(windowStartIndexArray[j]+windowSize, variantsNum); // inclusive
                    sb.setLength(0);
                    sb.append("1").append("\t");
                    sb.append(this.getPosition(windowStartSiteIndex)).append("\t");
                    sb.append(this.getPosition(windowEndSiteIndex-1)).append("\t");
                    sb.append(Doubles.join("\t", gridSource[i][j]));
                    bw.write(sb.toString());
                    bw.newLine();
                }
                bw.flush();
                bw.close();
            }
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
    }

    /**
     *
     * @param localAncestry dim1 is admixed taxon index, dim1 is n-way admixture (index seen Source.index)
     * @param localAncestryOutFile localAncestryOutFile
     * @param taxaGroupFile taxaGroupFile
     */
    public void write_localAncestry(BitSet[][] localAncestry, String localAncestryOutFile,
                                           String taxaGroupFile){
        TaxaGroup taxaGroup = TaxaGroup.buildFrom(taxaGroupFile);
        int variantsNum = this.getSiteNumber();
        try (BufferedWriter bw = IOTool.getWriter(localAncestryOutFile)) {
            List<String> admixedTaxaList = taxaGroup.getTaxaOf(Source.ADMIXED);
            StringBuilder sb = new StringBuilder();
//            sb.append("## native ancestry,introgressed ancestry 1, introgressed ancestry 2, et al.");
            sb.append("pos").append("\t");
            sb.append(String.join("\t", admixedTaxaList));
            bw.write(sb.toString());
            bw.newLine();
            int sourceNum = localAncestry[0].length;
            int ancestry;
            int pos;
            for (int variantIndex = 0; variantIndex < variantsNum; variantIndex++) {
                sb.setLength(0);
                pos = this.getSnps()[variantIndex].getPos();
                sb.append(pos).append("\t");
                for (BitSet[] bitSets : localAncestry) {
                    for (int sourceIndex = 0; sourceIndex < sourceNum; sourceIndex++) {
                        ancestry = bitSets[sourceIndex].get(variantIndex) ? 1 : 0;
                        sb.append(ancestry).append(",");
                    }
                    sb.deleteCharAt(sb.length() - 1);
                    sb.append("\t");
                }
                sb.deleteCharAt(sb.length()-1);
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public static void run_LAIDP(String genotypeFile, int windowSize, int stepSize, String taxaGroupFile,
                                 String ancestryAllele, int conjunctionNum,
                                 double switchCostScore, int maxSolutionCount, String localAnceOutFile, int threadsNum){
        GenotypeTable genotypeTable = new GenotypeTable(genotypeFile);
        BitSet[] ancestralAlleleBitSet = genotypeTable.getAncestralAlleleBitSet(ancestryAllele);
        BitSet[][] localAnc = genotypeTable.calculateLocalAncestry(windowSize, stepSize, taxaGroupFile,
                ancestralAlleleBitSet, conjunctionNum, switchCostScore, maxSolutionCount, threadsNum);
//        BitSet[][] localAnc = genotypeTable.calculateLocalAncestry_singleThread(windowSize, stepSize, taxaGroupFile,
//                ancestralAlleleBitSet, conjunctionNum, threadsNum);
        genotypeTable.write_localAncestry(localAnc, localAnceOutFile, taxaGroupFile);
    }

    /**
     *
     * @param ancestryAllele [simulation, file, int[]],
     * @return ancestryAllele bitset, dim1 represent alt allele is ancestral allele, dim2 mean ancestral allele is missing
     */
    public BitSet[] getAncestralAlleleBitSet(String ancestryAllele){
        if (ancestryAllele.equals("simulation")){
            BitSet[] ancestralAlleleBitSet = new BitSet[2];
            ancestralAlleleBitSet[0] = new BitSet();
            ancestralAlleleBitSet[1] = new BitSet();
            return ancestralAlleleBitSet;
        }
        List<String> temp= PStringUtils.fastSplit(ancestryAllele, ",");
        int[] outGroupTaxaIndices = new int[temp.size()];
        if (temp.size() > 1){
            for (int i = 0; i < temp.size(); i++) {
                outGroupTaxaIndices[i] = Integer.parseInt(temp.get(i));
            }
            return this.getAncestralAlleleFromTaxa(outGroupTaxaIndices);
        }else {
            return GenotypeTable.getAncestralAlleleBitSetFromFile(ancestryAllele);
        }
    }



    /**
     *
     * @param ancestryAlleleFile
     *  pos ref alt ancestralState
     *  2 A T 0
     *  3 C A 1
     *  4 G T -9
     *  ancestralState 0 means ancestral allele is reference allele
     *                 1 means ancestral allele is alt allele
     *                 -9 means ancestral allele is missing
     *  total pos number must be equal genotype file
     * @return ancestryAllele bitset, dim1 represent alt allele is ancestral allele,
     * dim2 mean ancestral allele is missing
     */
    public static BitSet[] getAncestralAlleleBitSetFromFile(String ancestryAlleleFile){
        BitSet[] ancestralAlleleBitset = new BitSet[2];
        for (int i = 0; i < ancestralAlleleBitset.length; i++) {
            ancestralAlleleBitset[i] = new BitSet();
        }
        try (BufferedReader br = IOTool.getReader(ancestryAlleleFile)) {
            br.readLine();
            String line;
            List<String> temp;
            int k = 0;
            int ancestralState;
            while ((line = br.readLine())!=null){
                temp =PStringUtils.fastSplit(line);
                ancestralState = Integer.parseInt(temp.get(3));
                if (ancestralState < 0){
                    ancestralAlleleBitset[1].set(k);
                }else if (ancestralState == 1){
                    ancestralAlleleBitset[0].set(k);
                }
                k++;
            }
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
        return ancestralAlleleBitset;
    }

    public BitSet[] getTaxaGenotype(int[] taxaIndices){
        BitSet[] genoTaxa = new BitSet[taxaIndices.length];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxa[i] = (BitSet) this.genoTaxon[taxaIndices[i]][0].clone();
        }
        return genoTaxa;
    }

    public List<String> getTaxaList(int[] taxaIndices){
        List<String> taxaList = new ArrayList<>();
        for (int index : taxaIndices){
            taxaList.add(this.taxa[index]);
        }
        return taxaList;
    }


    public GenotypeTable getSubGenotypeTableByTaxa(int[] taxaIndices){
        BitSet[][] genoTaxa = new BitSet[taxaIndices.length][2];
        for (int i = 0; i < taxaIndices.length; i++) {
            System.arraycopy(this.getGenoTaxon()[taxaIndices[i]], 0, genoTaxa[i], 0, genoTaxa[0].length);
        }
        String[] taxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            taxa[i]=this.getTaxa()[taxaIndices[i]];
        }
        SNP[] snps = new SNP[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            snps[i] = this.getSnps()[i].replicateWithoutFeature();
        }
        return GenotypeTable.newInstanceFromGenoTaxon(taxa, snps, genoTaxa);
    }

}
