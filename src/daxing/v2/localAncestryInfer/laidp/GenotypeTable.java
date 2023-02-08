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

    private static final int BLOCK_SIZE = 4000;

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
     * @param taxaIndices dim1 is different pop, dim2 is different taxon
     * @return mafs, dim1 is different pop, dim2 is variants
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
                int count1, total;
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
                        total = pop_bitSet[i].cardinality() - clonedBitSets[1].cardinality();
                        mafs[i][variantsIndex] = (double) count1 / total;
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
     * @return alt frequency of all sites
     */
    public double[] calculateAltAlleleFrequency(int threadsNum) {
        int taxaNum = this.getTaxa().length;
        int[][] taxaIndex = new int[1][];
        for (int i = 0; i < taxaIndex.length; i++) {
            taxaIndex[i] = IntStream.range(0, taxaNum).toArray();
        }
        return this.calculateAltAlleleFrequency(threadsNum, taxaIndex)[0];
    }

    public double[] calculateMaf(int threadsNum){
        double[] alts = this.calculateAltAlleleFrequency(threadsNum);
        double[] mafs = new double[alts.length];
        for (int i = 0; i < alts.length; i++) {
            if (alts[i] > 0.5){
                mafs[i] = 1-alts[i];
            }else {
                mafs[i] = alts[i];
            }
        }
        return mafs;
    }

    public double[][] calculateMaf(int threadsNum, int[][] taxaIndices){
        double[][] alts = this.calculateAltAlleleFrequency(threadsNum, taxaIndices);
        double[][] mafs = new double[alts.length][alts[0].length];
        for (int i = 0; i < alts.length; i++) {
            for (int j = 0; j < alts[i].length; j++) {
                if (alts[i][j] > 0.5){
                    mafs[i][j] = 1-alts[i][j];
                }else {
                    mafs[i][j] = alts[i][j];
                }
            }
        }
        return mafs;
    }

    /**
     *
     * @param threadsNum
     * @param taxaIndices
     * @param ancestralAlleleBitSet 1 in ancestralAlleleBitSet represent allele 1 is ancestral allele, dim1 is
     *                              genotype(haploid), dim2 is missing
     * @return derived allele frequency, -1 mean missing
     */
    public double[][] calculateDaf(int threadsNum, int[][] taxaIndices, BitSet[] ancestralAlleleBitSet){
        double[][] alts = this.calculateAltAlleleFrequency(threadsNum, taxaIndices);
        double[][] dafs = new double[alts.length][alts[0].length];
        for (int i = 0; i < alts.length; i++) {
            System.arraycopy(alts[i], 0, dafs[i], 0, alts[i].length);
        }
        for (int i = ancestralAlleleBitSet[0].nextSetBit(0); i >= 0; i = ancestralAlleleBitSet[0].nextSetBit(i+1)) {
            for (int j = 0; j < dafs.length; j++) {
                dafs[j][i] = 1 - alts[j][i];
            }
        }

        // missing will be filling with -1
        for (int i = ancestralAlleleBitSet[1].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[1].nextSetBit(i+1)) {
            for (int j = 0; j < dafs.length; j++) {
                dafs[j][i] = -1;
            }
        }
        return dafs;
    }

    /**
     *
     * @param threadsNum
     * @param ancestralAlleleBitSet
     * @return derived allele frequency, -1 mean missing
     */
    public double[] calculateDaf(int threadsNum, BitSet[] ancestralAlleleBitSet){
        double[] alts = this.calculateAltAlleleFrequency(threadsNum);
        double[] dafs = new double[alts.length];
        System.arraycopy(alts, 0, dafs, 0, alts.length);
        for (int i = ancestralAlleleBitSet[0].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[0].nextSetBit(i+1)) {
            dafs[i] = 1 - alts[i];
        }
        for (int i = ancestralAlleleBitSet[1].nextSetBit(0); i >=0; i = ancestralAlleleBitSet[1].nextSetBit(i+1)) {
            dafs[i] = -1;
        }
        return dafs;
    }

    /**
     *
     * @param taxonIndex outGroup population
     * @return ancestralAlleleBitSet
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
