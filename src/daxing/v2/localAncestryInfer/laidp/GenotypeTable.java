package daxing.v2.localAncestryInfer.laidp;

import daxing.common.bisnp.SNP;
import daxing.common.chrrange.ChrPos;
import daxing.common.utiles.IOTool;
import org.apache.commons.lang.ArrayUtils;
import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
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
            List<String> vcfAnnotationList = new ArrayList<>();
            BufferedReader br = IOTool.getReader(haploidGenotypeFile);
            String temp;
            while ((temp = br.readLine()).startsWith("##")) {
                vcfAnnotationList.add(temp);
            }
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
            List<String> lines = new ArrayList ();
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
//        this.sortByTaxa();
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

    public static GenotypeTable newInstanceFromGenoSite(String[] taxa, SNP[] snps, BitSet[][] genoSite){
        GenotypeTable genotypeTable = new GenotypeTable(taxa, snps, genoSite, null);
        genotypeTable.transposeSiteToTaxon();
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

//    public void sortByTaxa() {
//        System.out.println("Start sorting genotype table by taxon");
//        long start = System.nanoTime();
//        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
//        Arrays.sort(this.taxa);
//        BitSet[][] nGenoTaxon = new BitSet[this.getTaxaNumber()][];
//        for (int i = 0; i < this.getTaxaNumber(); i++) {
//            nGenoTaxon[i] = this.genoTaxon[indices[i]];
//        }
//        this.genoTaxon = nGenoTaxon;
//        System.out.println("Sorting finished in " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
//        this.transposeTaxonToSite();
//    }

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
        double[][] mafs = new double[taxaIndices.length][numVariants];
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
                        System.arraycopy(bitSets, 0, clonedBitSets, 0, bitSetsLen);
                        clonedBitSets[0].and(pop_bitSet[i]);
                        clonedBitSets[1].and(pop_bitSet[i]);
                        count1 = clonedBitSets[0].cardinality();
                        missing =  clonedBitSets[1].cardinality();
                        total = pop_bitSet[i].cardinality();
                        if (missing == total){
                            // all sample in these pop is missing
                            mafs[i][variantsIndex] = -1;
                        }else {
                            mafs[i][variantsIndex] = (double) count1 / (total - missing);
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
        return mafs;
    }

    /**
     *
     * @param threadsNum threadsNum
     * @return alternative allele frequency, all taxa will be treated as one population
     * -1 means sample in one population is all missing
     */
    public double[] calculateAltAlleleFrequency(int threadsNum) {
        int taxaNum = this.getTaxa().length;
        int[][] taxaIndex = new int[1][];
        for (int i = 0; i < taxaIndex.length; i++) {
            taxaIndex[i] = IntStream.range(0, taxaNum).toArray();
        }
        return this.calculateAltAlleleFrequency(threadsNum, taxaIndex)[0];
    }

    /**
     *
     * @param threadsNum threadsNum
     * @return minor allele frequency, all taxa will be treated as one population.
     * -1 means sample in one population is all missing
     */
    public double[] calculateMaf(int threadsNum){
        double[] alts = this.calculateAltAlleleFrequency(threadsNum);
        double[] mafs = new double[alts.length];
        for (int i = 0; i < alts.length; i++) {
            if (alts[i] > 0.5){
                mafs[i] = 1-alts[i];
            }else {
                // include missing in alts
                mafs[i] = alts[i];
            }
        }
        return mafs;
    }

    /**
     *
     * @param threadsNum threadsNum
     * @param taxaIndices dim1 is different populations, dim2 is different taxon
     * @return alternative allele frequency, dim1 is different populations, dim2 is variants.
     * -1 means sample in one population is all missing
     */
    public double[][] calculateMaf(int threadsNum, int[][] taxaIndices){
        double[][] alts = this.calculateAltAlleleFrequency(threadsNum, taxaIndices);
        double[][] mafs = new double[alts.length][alts[0].length];
        for (int i = 0; i < alts.length; i++) {
            for (int j = 0; j < alts[i].length; j++) {
                if (alts[i][j] > 0.5){
                    mafs[i][j] = 1-alts[i][j];
                }else {
                    // include missing in alts
                    mafs[i][j] = alts[i][j];
                }
            }
        }
        return mafs;
    }

    /**
     *
     * @param threadsNum threadsNum
     * @param taxaIndices dim1 is different populations, dim2 is different taxon
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @return derived allele frequency, -1 mean missing.
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
     * @param threadsNum threadsNum
     * @param ancestralAlleleBitSet dim1 is genotype(haploid), dim2 is missing.
     * 1 in ancestralAlleleBitSet represent alternative allele is ancestral allele.
     * @return derived allele frequency,  all taxa will be treated as one population.
     * -1 mean missing
     * Missing include two cases:
     * case1: sample in one population is all missing
     * case2: ancestral allele is missing
     */
    public double[] calculateDaf(int threadsNum, BitSet[] ancestralAlleleBitSet){
        double[] alts = this.calculateAltAlleleFrequency(threadsNum);
        double[] dafs = new double[alts.length];
        System.arraycopy(alts, 0, dafs, 0, alts.length);
        for (int i = ancestralAlleleBitSet[0].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[0].nextSetBit(i+1)) {
            // missing
            if (alts[i] < 0){
                dafs[i] = -1;
            }else {
                dafs[i] = 1 - alts[i];
            }
        }
        for (int i = ancestralAlleleBitSet[1].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[1].nextSetBit(i+1)) {
            dafs[i] = -1;
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
     *
     * @param stepSize step size of sliding window, the unit is variants, not bp
     * @return windowStartIndex array, site start index of all windows
     */
    public int[] getWindowStartIndex(int stepSize){
        int numVariants = this.getSiteNumber();
        int numberWindow = (numVariants + stepSize - 1) / stepSize;
        int[] window = new int[numberWindow];
        for (int windowIndex = 0; windowIndex < numberWindow; windowIndex++) {
            window[windowIndex] = Math.min(windowIndex * stepSize, numVariants);
        }
        return window;
    }

    /**
     *
     * @param pop_taxonIndex dim1 is different populations, dim2 is different taxa
     * @return pairwise dxy between all populations
     */
    public double[] calculateDxy(int[][] pop_taxonIndex){
        int pop_num = pop_taxonIndex.length;
        int hapCountA, hapCountB;
        int len_haplotype = this.snps[this.getSiteNumber()-1].getPos()-this.snps[0].getPos() + 1;
        double cumDifference;
        BitSet bitSet, missing;
        double[] dxy = new double[pop_num*(pop_num-1)/2];
        Arrays.fill(dxy, -1);
        int k =0 ;
        for (int popIndexA = 0; popIndexA < pop_num -1; popIndexA++) {
            for (int popIndexB = popIndexA+1; popIndexB < pop_num; popIndexB++) {
                hapCountA = pop_taxonIndex[popIndexA].length;
                hapCountB = pop_taxonIndex[popIndexB].length;
                cumDifference = 0;
                for (int taxonA : pop_taxonIndex[popIndexA]){
                    for (int taxonB : pop_taxonIndex[popIndexB]){
                        bitSet = (BitSet) this.genoTaxon[taxonA][0].clone();
                        missing = (BitSet)this.genoTaxon[taxonA][1].clone();
                        missing.or(this.genoTaxon[taxonB][1]);
                        bitSet.xor(this.genoTaxon[taxonB][0]);
                        bitSet.andNot(missing);
                        cumDifference+=bitSet.cardinality();
                    }
                }
                dxy[k++] = cumDifference/(hapCountA*hapCountB)/len_haplotype;
            }
        }
        return dxy;
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
        SNP snp_start, snp_end;
        int[] windowLen = new int[windowStartIndexArray.length];
        for (int i = 0; i < windowLen.length; i++) {
            snp_start = this.snps[windowStartIndexArray[i]];
            snp_end = this.snps[Math.min(windowStartIndexArray[i]+windowSize, this.getSiteNumber())-1];
            windowLen[i] = snp_end.getPos() - snp_start.getPos() + 1;
        }
        double[][] dxyArrays = new double[windowStartIndexArray.length][pop_num*(pop_num-1)/2];
        int hapCountA, hapCountB, hapCountAB;
        double[] cumDifference;
        BitSet bitSet, missing, subBitset;
        int k =0;
        for (int popIndexA = 0; popIndexA < pop_num -1; popIndexA++) {
            for (int popIndexB = popIndexA+1; popIndexB < pop_num; popIndexB++) {
                hapCountA = pop_taxonIndex[popIndexA].length;
                hapCountB = pop_taxonIndex[popIndexB].length;
                cumDifference = new double[windowStartIndexArray.length];
                for (int taxonA : pop_taxonIndex[popIndexA]){
                    for (int taxonB : pop_taxonIndex[popIndexB]){
                        bitSet = (BitSet) this.genoTaxon[taxonA][0].clone();
                        missing = (BitSet)this.genoTaxon[taxonA][1].clone();
                        missing.or(this.genoTaxon[taxonB][1]);
                        bitSet.xor(this.genoTaxon[taxonB][0]);
                        bitSet.andNot(missing);
                        for (int i = 0; i < windowStartIndexArray.length; i++) {
                            subBitset = bitSet.get(windowStartIndexArray[i], windowStartIndexArray[i]+windowSize);
                            cumDifference[i] += subBitset.cardinality();
                        }
                    }
                }
                hapCountAB = hapCountA*hapCountB;
                for (int i = 0; i < cumDifference.length; i++) {
                    dxyArrays[i][k] = cumDifference[i]/hapCountAB/windowLen[i];
                }
                k++;
            }
        }
        return dxyArrays;
    }


    /**
     *
     * @param taxaIndices taxon可能是不连续的, 也可能是未排序的
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return double[][], the first dimension is taxa, the second dimension is genotype (missing = 0.5, alt=1, ref=0)
     */
    public double[][] getSrcGenotypeFrom(int[] taxaIndices, int[] siteIndices){
        BitSet[][] genoTaxa = new BitSet[taxaIndices.length][];
        for (int i = 0; i < taxaIndices.length; i++) {
            genoTaxa[i] = new BitSet[2];
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, genoTaxa[i], 0, genoTaxa[0].length);
        }
        BitSet[][] genoTaxaSite=new BitSet[taxaIndices.length][2];
        for (int i = 0; i < genoTaxa.length; i++) {
            for (int j = 0; j < genoTaxa[i].length; j++) {
                genoTaxaSite[i][j]=genoTaxa[i][j].get(siteIndices[0],siteIndices[1]);
            }
        }
        double[][] genotype=new double[taxaIndices.length][];
        for (int i = 0; i < genotype.length; i++) {
            genotype[i]= new double[siteIndices[1]-siteIndices[0]+1];
            Arrays.fill(genotype[i], 0);
            for (int j = 0; j < genoTaxaSite[i][0].length(); j++) {
                if (genoTaxaSite[i][0].get(j)){
                    genotype[i][j]=1;
                }
            }
            for (int j = 0; j < genoTaxaSite[i][1].length(); j++) {
                if (genoTaxaSite[i][1].get(j)){
                    genotype[i][j]=0.5;
                }
            }
        }
        return genotype;
    }

    /**
     *
     * @param taxaIndices taxon可能是不连续的, 也可能是未排序的
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return BitSet[], the first dimension is taxa, the bitset dimension is genotype
     */
    public BitSet[] getSrcGenotypeBitSetFrom(int[] taxaIndices, int[] siteIndices){
        BitSet[][] genoTaxa = new BitSet[taxaIndices.length][];
        for (int i = 0; i < taxaIndices.length; i++) {
            genoTaxa[i] = new BitSet[2];
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, genoTaxa[i], 0, genoTaxa[0].length);
        }
        BitSet[] genoTaxaSite=new BitSet[taxaIndices.length];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxaSite[i] = genoTaxa[i][0].get(siteIndices[0],siteIndices[1]);
        }
        return genoTaxaSite;
    }

    /**
     *
     * @param taxonIndex taxon index
     * @param siteIndices 是连续的range, int[2] [start,end] [include,include]
     * @return double[], genotype (missing = 0.5, alt=1, ref=0)
     */
    public double[] getQueryGenotypeFrom(int taxonIndex, int[] siteIndices){
        BitSet[] genoTaxa = new BitSet[2];
        System.arraycopy(this.genoTaxon[taxonIndex], 0, genoTaxa, 0, genoTaxa.length);
        BitSet[] genoTaxaSite=new BitSet[2];
        for (int i = 0; i < genoTaxa.length; i++) {
            genoTaxaSite[i]=genoTaxa[i].get(siteIndices[0],siteIndices[1]);
        }
        double[] query=new double[siteIndices[1]-siteIndices[0]+1];
        Arrays.fill(query, 0);
        for (int i = 0; i < genoTaxaSite[0].length(); i++) {
            if (genoTaxaSite[0].get(i)){
                query[i] = 1;
            }
        }
        for (int i = 0; i < genoTaxaSite[1].length(); i++) {
            if (genoTaxaSite[1].get(i)){
                query[i]=0.5;
            }
        }
        return query;
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
