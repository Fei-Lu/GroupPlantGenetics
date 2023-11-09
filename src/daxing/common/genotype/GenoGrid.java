package daxing.common.genotype;

import cern.colt.GenericSorting;
import cern.colt.Swapper;
import cern.colt.function.IntComparator;
import daxing.common.genotype.facilities.BlockBinary;
import daxing.common.genotype.facilities.BlockVCF;
import daxing.common.genotype.facilities.SiteGeno;
import daxing.common.utiles.IOTool;
import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeExport;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.Benchmark;
import pgl.infra.utils.PArrayUtils;
import pgl.infra.utils.PStringUtils;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

/**
 * Modify from GenotypeGrid, performance improvement, code simplification, easy customization
 * @author Daxing Xu
 */
public class GenoGrid extends AbstractGenotypeTable implements Swapper, IntComparator {

    String[] taxa = null;
    BiSNP[] snps = null;
    float[] mafs = null;

    public enum Direction {
        BySite, ByTaxon
    }
    /**
     * Bit genotype by site, the first dimension is site; the second dimension is phase 1, phase 2, and missing
     */
    BitSet[][] genoSite = null;
    /**
     * Bit genotype by taxon, the first dimension is taxon; the second dimension is phase 1, phase 2, and missing
     */
    BitSet[][] genoTaxon = null;

    /**
     * Construct an object by reading a file
     * @param infileS inputFile
     * @param format format
     */
    public GenoGrid (String infileS, GenoIOFormat format) {
        switch (format) {
            case VCF:
            case VCF_GZ:
                this.buildFromVCF(infileS);
                break;
            case Binary:
            case Binary_GZ:
                this.buildFromBinary(infileS);
                break;
            case HDF5:
                throw new UnsupportedOperationException("Not supported yet.");
        }
        this.sortByTaxa();
    }

    public GenoGrid (BitSet[][] bArray, Direction gd, String[] taxa, BiSNP[] snps) {
        this.taxa = taxa;
        this.snps = snps;
        if (gd == Direction.BySite) {
            genoSite = bArray;
            this.transposeSiteToTaxon();
        }
        else if (gd == Direction.ByTaxon) {
            genoTaxon = bArray;
            this.transposeTaxonToSite();
        }
        this.mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    @Override
    public int getTaxaNumber() {
        return this.taxa.length;
    }

    @Override
    public int getSiteNumber() {
        return this.snps.length;
    }

    @Override
    public String getTaxonName(int taxonIndex) {
        return this.taxa[taxonIndex];
    }

    @Override
    public String[] getTaxaNames() {
        return this.taxa;
    }

    @Override
    public short getChromosome(int siteIndex) {
        return this.snps[siteIndex].getChromosome();
    }

    @Override
    public int getPosition(int siteIndex) {
        return this.snps[siteIndex].getPosition();
    }

    @Override
    public void sortBySite() {
        System.out.println("Start sorting genotype table by site");
        long start = System.nanoTime();
        GenericSorting.quickSort(0, this.getSiteNumber(), this, this);
        System.out.println("Sorting finished in " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
        this.transposeSiteToTaxon();
    }

    @Override
    public void sortByTaxa() {
        System.out.println("Start sorting genotype table by taxon");
        long start = System.nanoTime();
        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
        Arrays.sort(this.taxa);
        BitSet[][] nGenoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            nGenoTaxon[i] = this.genoTaxon[indices[i]];
        }
        this.genoTaxon = nGenoTaxon;
        System.out.println("Sorting finished in " + Benchmark.getTimeSpanSeconds(start) + " seconds.");
        this.transposeTaxonToSite();
    }

    @Override
    public int getTaxonIndex(String taxon) {
        return Arrays.binarySearch(taxa, taxon);
    }

    @Override
    public int getSiteIndex(short chromosome, int position) {
        ChrPos query = new ChrPos (chromosome, position);
        return Arrays.binarySearch(this.snps, query);
    }

    @Override
    public byte getGenotypeByte(int siteIndex, int taxonIndex) {
        if (isMissing(siteIndex, taxonIndex)) return AlleleEncoder.genotypeMissingCoding;
        byte ref = this.getReferenceAlleleByte(siteIndex);
        byte alt = this.getAlternativeAlleleByte(siteIndex);
        byte b1;
        byte b2;
        if (isPhase1Alternative(siteIndex, taxonIndex)) b1 = alt;
        else b1 = ref;
        if (isPhase2Alternative(siteIndex, taxonIndex)) b2 = alt;
        else b2 = ref;
        return AlleleEncoder.getGenotypeCoding(b1, b2);
    }

    @Override
    public boolean isHeterozygous(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex,taxonIndex)) return false;
        return this.isPhase1Alternative(siteIndex, taxonIndex) != this.isPhase2Alternative(siteIndex, taxonIndex);
    }

    @Override
    public boolean isHomozygous(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex,taxonIndex)) return false;
        return this.isPhase1Alternative(siteIndex, taxonIndex) == this.isPhase2Alternative(siteIndex, taxonIndex);
    }

    @Override
    public boolean isMissing(int siteIndex, int taxonIndex) {
        return this.genoSite[siteIndex][2].get(taxonIndex);
    }

    @Override
    public boolean isPhase1Alternative(int siteIndex, int taxonIndex) {
        return this.genoSite[siteIndex][0].get(taxonIndex);
    }

    @Override
    public boolean isPhase2Alternative(int siteIndex, int taxonIndex) {
        return this.genoSite[siteIndex][1].get(taxonIndex);
    }

    @Override
    public boolean isPhase1Reference(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex, taxonIndex)) return false;
        return !this.isPhase1Alternative(siteIndex, taxonIndex);
    }

    @Override
    public boolean isPhase2Reference(int siteIndex, int taxonIndex) {
        if (this.isMissing(siteIndex, taxonIndex)) return false;
        return !this.isPhase2Alternative(siteIndex, taxonIndex);
    }

    @Override
    public boolean isAlternativeAlleleTypeOf(AlleleType at, int siteIndex) {
        return snps[siteIndex].isAlternativeAlleleTypeOf(at);
    }

    @Override
    public boolean isReferenceAlleleTypeOf(AlleleType at, int siteIndex) {
        return snps[siteIndex].isReferenceAlleleTypeOf(at);
    }

    @Override
    public int getMissingNumberBySite(int siteIndex) {
        return this.genoSite[siteIndex][2].cardinality();
    }

    @Override
    public int getMissingNumberByTaxon(int taxonIndex) {
        return this.genoTaxon[taxonIndex][2].cardinality();
    }

    @Override
    public int getNonMissingNumberBySite(int siteIndex) {
        return this.getTaxaNumber() - this.getMissingNumberBySite(siteIndex);
    }

    @Override
    public int getNonMissingNumberByTaxon(int taxonIndex) {
        return this.getSiteNumber() - this.getMissingNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHomozygoteNumberBySite(int siteIndex) {
        return this.getNonMissingNumberBySite(siteIndex) - this.getHeterozygoteNumberBySite(siteIndex);
    }

    @Override
    public int getHomozygoteNumberByTaxon(int taxonIndex) {
        return this.getNonMissingNumberByTaxon(taxonIndex) - this.getHeterozygoteNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHeterozygoteNumberBySite(int siteIndex) {
        BitSet phase1C = this.genoSite[siteIndex][0].get(0, this.getTaxaNumber());
        phase1C.xor(this.genoSite[siteIndex][1]);
        return phase1C.cardinality();
    }

    @Override
    public int getHeterozygoteNumberByTaxon(int taxonIndex) {
        BitSet phase1C = this.genoTaxon[taxonIndex][0].get(0, this.getSiteNumber());
        phase1C.xor(this.genoTaxon[taxonIndex][1]);
        return phase1C.cardinality();
    }

    @Override
    public int getAlternativeAlleleNumberBySite(int siteIndex) {
        int cnt = 0;
        for (int i = 0; i < 2; i++) {
            cnt+=this.genoSite[siteIndex][i].cardinality();
        }
        return cnt;
    }

    @Override
    public int getAlternativeAlleleOccurrenceBySite (int siteIndex) {
        BitSet bs = (BitSet)this.genoSite[siteIndex][0].clone();
        bs.or(this.genoSite[siteIndex][1]);
        return bs.cardinality();
    }

    @Override
    public float getHeterozygousProportionByTaxon(int taxonIndex) {
        return (float)((double)this.getHeterozygoteNumberByTaxon(taxonIndex)/this.getNonMissingNumberByTaxon(taxonIndex));
    }

    @Override
    public float getHeterozygousProportionBySite(int siteIndex) {
        return (float)((double)this.getHeterozygoteNumberBySite(siteIndex)/this.getNonMissingNumberBySite(siteIndex));
    }

    @Override
    public byte getMinorAlleleByte(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleByte(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleByte(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMinorAlleleByte(siteIndex);
    }

    @Override
    public char getMinorAlleleBase(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Minor)) return this.getReferenceAlleleBase(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Minor)) return this.getAlternativeAlleleBase(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMinorAlleleBase(siteIndex);
    }

    @Override
    public float getMinorAlleleFrequency(int siteIndex) {
        if (this.mafs[siteIndex] != Float.MIN_VALUE) return this.mafs[siteIndex];
        float altFre = this.getAlternativeAlleleFrequency(siteIndex);
        if (altFre < 0.5) {
            this.setAlternativeAlleleType(AlleleType.Minor, siteIndex);
            this.setReferenceAlleleType(AlleleType.Major, siteIndex);
        }
        else {
            this.setReferenceAlleleType(AlleleType.Minor, siteIndex);
            this.setAlternativeAlleleType(AlleleType.Major, siteIndex);
            altFre = 1 - altFre;
        }
        mafs[siteIndex] = altFre;
        return altFre;
    }

    @Override
    public byte getMajorAlleleByte(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleByte(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleByte(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMajorAlleleByte(siteIndex);
    }

    @Override
    public char getMajorAlleleBase(int siteIndex) {
        if (this.snps[siteIndex].reference.isAlleleTypeOf(AlleleType.Major)) return this.getReferenceAlleleBase(siteIndex);
        if (this.snps[siteIndex].alternative.isAlleleTypeOf(AlleleType.Major)) return this.getAlternativeAlleleBase(siteIndex);
        this.getMinorAlleleFrequency(siteIndex);
        return this.getMajorAlleleBase(siteIndex);
    }

    @Override
    public float getMajorAlleleFrequency(int siteIndex) {
        return 1-this.getMinorAlleleFrequency(siteIndex);
    }

    @Override
    public byte getReferenceAlleleByte(int siteIndex) {
        return snps[siteIndex].getReferenceAlleleByte();
    }

    @Override
    public char getReferenceAlleleBase(int siteIndex) {
        return snps[siteIndex].getReferenceAlleleBase();
    }

    @Override
    public float getReferenceAlleleFrequency(int siteIndex) {
        if (this.isReferenceAlleleTypeOf(AlleleType.Minor, siteIndex)) return this.mafs[siteIndex];
        if (this.isReferenceAlleleTypeOf(AlleleType.Major, siteIndex)) return 1-this.mafs[siteIndex];
        this.getMinorAlleleFrequency(siteIndex);
        return this.getReferenceAlleleFrequency(siteIndex);
    }

    @Override
    public byte getAlternativeAlleleByte(int siteIndex) {
        return snps[siteIndex].getAlternativeAlleleByte();
    }

    @Override
    public char getAlternativeAlleleBase(int siteIndex) {
        return snps[siteIndex].getAlternativeAlleleBase();
    }

    @Override
    public float getAlternativeAlleleFrequency(int siteIndex) {
        if (this.isAlternativeAlleleTypeOf(AlleleType.Minor, siteIndex)) return mafs[siteIndex];
        if (this.isAlternativeAlleleTypeOf(AlleleType.Major, siteIndex)) return 1-mafs[siteIndex];
        return (float)((double)this.getAlternativeAlleleNumberBySite(siteIndex)/(this.getNonMissingNumberBySite(siteIndex)*2));
    }

    @Override
    public int getStartIndexOfChromosome(short chromosome) {
        int index = this.getSiteIndex(chromosome, Integer.MIN_VALUE);
        if (index < 0) {
            index = -index - 1;
            if (index < this.getSiteNumber() && this.getChromosome(index) == chromosome) return index;
            return -1;
        }
        else {
            while (index > 0 && this.getChromosome(index-1) == chromosome) {
                index--;
            }
            return index;
        }
    }

    @Override
    public int getEndIndexOfChromosome(short chromosome) {
        int index = this.getSiteIndex(chromosome, Integer.MAX_VALUE);
        if (index < 0) {
            index = -index - 2;
            if (this.getChromosome(index) == chromosome) return index+1;
            else return -1;
        }
        else {
            while ((index+1) < this.getSiteNumber() && this.getChromosome(index+1) == chromosome) {
                index++;
            }
            return index+1;
        }
    }

    public int getIBSDistanceNonMissingSite(int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex){
        BitSet bothMissing = this.genoTaxon[taxonIndex1][2].get(startSiteIndex, endSiteIndex);
        bothMissing.or(this.genoTaxon[taxonIndex2][2].get(startSiteIndex, endSiteIndex));
        return (endSiteIndex-startSiteIndex-bothMissing.cardinality());
    }

    public int getIBSDistanceNonMissingSite(int taxonIndex1, int taxonIndex2){
        return this.getIBSDistanceNonMissingSite(taxonIndex1, taxonIndex2, 0, this.getSiteNumber());
    }

    /**
     * parent1和parent2分离的site总数
     * @param taxonIndex1
     * @param taxonIndex2
     * @param startSiteIndex
     * @param endSiteIndex
     * @return
     */
    public int[] getSegregatingSitesIndexBetweenParents(int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex){
        BitSet bothHetero = this.genoTaxon[taxonIndex1][0].get(startSiteIndex, endSiteIndex);
        bothHetero.xor(this.genoTaxon[taxonIndex1][1].get(startSiteIndex, endSiteIndex));
        BitSet taxon2Hetero = this.genoTaxon[taxonIndex2][0].get(startSiteIndex, endSiteIndex);
        taxon2Hetero.xor(this.genoTaxon[taxonIndex2][1].get(startSiteIndex, endSiteIndex));

        // both hetero sites
        bothHetero.or(taxon2Hetero);

        // bothMissing sites
        BitSet bothMissing = this.genoTaxon[taxonIndex1][2].get(startSiteIndex, endSiteIndex);
        bothMissing.or(this.genoTaxon[taxonIndex2][2].get(startSiteIndex, endSiteIndex));

        // segregation sites
        BitSet segregating = this.genoTaxon[taxonIndex1][0].get(startSiteIndex, endSiteIndex);
        segregating.xor(this.genoTaxon[taxonIndex2][0].get(startSiteIndex, endSiteIndex));
        segregating.andNot(bothMissing);
        segregating.andNot(bothHetero);

        return segregating.stream().toArray();
    }

    /**
     * parent1和parent2分离的site总数
     * @param taxonIndex1
     * @param taxonIndex2
     * @return
     */
    public int[] getSegregatingSitesIndexBetweenParents(int taxonIndex1, int taxonIndex2){
        return this.getSegregatingSitesIndexBetweenParents(taxonIndex1, taxonIndex2, 0, this.getSiteNumber());
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2) {
        return this.getIBSDistance(taxonIndex1, taxonIndex2, 0, this.getSiteNumber());
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex) {
        BitSet[] t1 = new BitSet[2];
        BitSet[] t2 = new BitSet[2];
        BitSet bothMissing = this.genoTaxon[taxonIndex1][2].get(startSiteIndex, endSiteIndex);
        bothMissing.or(this.genoTaxon[taxonIndex2][2].get(startSiteIndex, endSiteIndex));
        for (int i = 0; i < t1.length; i++) {
            t1[i] = this.genoTaxon[taxonIndex1][i].get(startSiteIndex, endSiteIndex);
            t1[i].or(bothMissing);
            t2[i] = this.genoTaxon[taxonIndex2][i].get(startSiteIndex, endSiteIndex);
            t2[i].or(bothMissing);
        }
        int cnt = 0;
        for (BitSet set : t1) {
            for (BitSet bitSet : t2) {
                BitSet bs = (BitSet) set.clone();
                bs.xor(bitSet);
                cnt += bs.cardinality();
            }
        }
        return (float) ((double)cnt/4/(endSiteIndex-startSiteIndex-bothMissing.cardinality()));
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int[] siteIndices) {
        BitSet selected = new BitSet(this.getSiteNumber());
        for (int siteIndex : siteIndices) {
            selected.set(siteIndex);
        }
        int cnt = 0;
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                BitSet bs1 = this.genoTaxon[taxonIndex1][i].get(0, this.getSiteNumber());
                bs1.and(selected);
                BitSet bs2 = this.genoTaxon[taxonIndex2][i].get(0, this.getSiteNumber());
                bs2.and(selected);
                bs1.xor(bs2);
                cnt+=bs1.cardinality();
            }
        }
        BitSet bs1 = this.genoTaxon[taxonIndex1][2].get(0, this.getSiteNumber());
        bs1.and(selected);
        BitSet bs2 = this.genoTaxon[taxonIndex2][2].get(0, this.getSiteNumber());
        bs2.and(selected);
        bs1.or(bs2);
        return (float) ((double)cnt/4/(siteIndices.length-bs1.cardinality()));
    }

    @Override
    public float[][] getIBSDistanceMatrix() {
        return this.getIBSDistanceMatrix(0, this.getSiteNumber());
    }

    @Override
    public float[][] getIBSDistanceMatrix(int startIndex, int endIndex) {
        float[][] matrix = new float[this.getTaxaNumber()][this.getTaxaNumber()];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < matrix.length-1; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            for (int j = i + 1; j < this.getTaxaNumber(); j++) {
                matrix[i][j] = this.getIBSDistance(i, j, startIndex, endIndex);
                matrix[j][i] = matrix[i][j];
            }
        });
        return matrix;
    }

    @Override
    public float[][] getIBSDistanceMatrix(int[] siteIndices) {
        float[][] matrix = new float[this.getTaxaNumber()][this.getTaxaNumber()];
        List<Integer> indexList = new ArrayList<>();
        for (int i = 0; i < matrix.length-1; i++) {
            indexList.add(i);
        }
        indexList.parallelStream().forEach(i -> {
            for (int j = i + 1; j < this.getTaxaNumber(); j++) {
                matrix[i][j] = this.getIBSDistance(i, j, siteIndices);
                matrix[j][i] = matrix[i][j];
            }
        });
        return matrix;
    }

    @Override
    public GenotypeTable getSubGenotypeTableBySite(int[] siteIndices) {
        BitSet[][] bArray = new BitSet[siteIndices.length][3];
        for (int i = 0; i < siteIndices.length; i++) {
            System.arraycopy(this.genoSite[siteIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        BiSNP[] nsnps = new BiSNP[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            nsnps[i] = this.snps[siteIndices[i]];
        }
        return new GenoGrid(bArray, Direction.BySite, this.taxa, nsnps);
    }

    @Override
    public GenotypeTable getSubGenotypeTableByTaxa(int[] taxaIndices) {
        BitSet[][] bArray = new BitSet[taxaIndices.length][3];
        for (int i = 0; i < taxaIndices.length; i++) {
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            nTaxa[i] = this.getTaxonName(taxaIndices[i]);
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            nsnps[i] = this.snps[i].replicateWithoutFeature();
        }
        return new GenoGrid(bArray, Direction.ByTaxon, nTaxa, nsnps);
    }

    @Override
    public String getUnphasedVCFRecord(int siteIndex) {
        StringBuilder vsb = new StringBuilder();
        char delimiter = '/';
        vsb.setLength(0);
        vsb.append(this.getChromosome(siteIndex)).append("\t").append(this.getPosition(siteIndex)).append("\t").append(this.getChromosome(siteIndex)).append("-").append(this.getPosition(siteIndex)).append("\t");
        vsb.append(this.getReferenceAlleleBase(siteIndex)).append("\t").append(this.getAlternativeAlleleBase(siteIndex)).append("\t.\t.\t");
        if (this.snps[siteIndex].getSNPInfo() == null) vsb.append(".");
        else vsb.append(this.snps[siteIndex].getSNPInfo());
        vsb.append("\t").append("GT");
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            if (isMissing(siteIndex, i)) vsb.append("\t").append(".").append(delimiter).append(".");
            else {
                vsb.append("\t");
                if (isPhase1Alternative(siteIndex, i)) vsb.append("1");
                else vsb.append("0");
                vsb.append(delimiter);
                if (isPhase2Alternative(siteIndex, i)) vsb.append("1");
                else vsb.append("0");
            }
        }
        return vsb.toString();
    }

    @Override
    public ByteBuffer getBinaryOutput(int siteIndex, ByteBuffer bb) {
        bb.putShort(this.getChromosome(siteIndex));
        bb.putInt(this.getPosition(siteIndex));
        bb.put(AlleleEncoder.getGenotypeCoding(this.getReferenceAlleleByte(siteIndex), this.getAlternativeAlleleByte(siteIndex)));
        bb.put(this.snps[siteIndex].getReferenceAlleleFeature());
        bb.put(this.snps[siteIndex].getAlternativeAlleleFeature());
        int size = (bb.capacity()-GenotypeExport.getByteSizeOfSNPInBinary())/3;
        for (int i = 0; i < this.genoSite[siteIndex].length; i++) {
            byte[] b = new byte[size];
            byte[] ob = this.genoSite[siteIndex][i].toByteArray();
            System.arraycopy(ob,0,b,0,ob.length);
            bb.put(b);
        }
        return bb;
    }

    @Override
    public void setAlternativeAlleleType(AlleleType at, int siteIndex) {
        this.snps[siteIndex].setAlternativeAlleleType(at);
    }

    @Override
    public void setReferenceAlleleType(AlleleType at, int siteIndex) {
        this.snps[siteIndex].setReferenceAlleleType(at);
    }

    /**
     * Reader of a binary genotype file
     * @param infileS vcf binaryFile
     */
    private void buildFromBinary (String infileS) {
        try{
            DataInputStream dis=IOTool.getBinaryReader(infileS);
            int siteNumber = dis.readInt();
            short taxaNumber = (short)dis.readInt();
            this.taxa = new String[taxaNumber];
            for (int i = 0; i < taxaNumber; i++) {
                this.taxa[i] = dis.readUTF();
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<BlockBinary>> resultList = new ArrayList<>();
            byte[][] input = new byte[BlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
            int siteCount = 0;
            int startIndex = 0;
            int actBlockSize = 0;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < siteNumber; i++) {
                dis.read(input[actBlockSize]);
                siteCount++;
                actBlockSize++;
                if (siteCount% BlockBinary.blockSize == 0) {
                    BlockBinary gsb = new BlockBinary(input, startIndex, actBlockSize);
                    Future<BlockBinary> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=actBlockSize;
                    actBlockSize = 0;
                    input = new byte[BlockBinary.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
                }
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb);
                }
            }
            dis.close();
            if (actBlockSize != 0) {
                BlockBinary gsb = new BlockBinary(input, startIndex, actBlockSize);
                Future<BlockBinary> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new BiSNP[siteCount];
            for (Future<BlockBinary> blockBinaryFuture : resultList) {
                BlockBinary block = blockBinaryFuture.get();
                for (int j = 0; j < block.getActBlockSize(); j++) {
                    this.snps[block.getStartIndex() + j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex() + j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(this.getSiteNumber()).append(" SNPs are in ").append(infileS).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
        mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    private void buildFromVCF (String infileS) {
        try {
            List<String> vcfAnnotationList = new ArrayList<>();
            BufferedReader br = IOTool.getReader(infileS);
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
            List<Future<BlockVCF>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            List<String> lines = new ArrayList ();
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()%BlockVCF.blockSize == 0) {
                    BlockVCF gsb = new BlockVCF(lines, startIndex);
                    Future<BlockVCF> result = pool.submit(gsb);
                    resultList.add(result);
                    startIndex+=lines.size();
                    lines = new ArrayList<>();
                }
                siteCount++;
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb);
                }
            }
            br.close();
            if (lines.size() != 0) {
                BlockVCF gsb = new BlockVCF(lines, startIndex);
                Future<BlockVCF> result = pool.submit(gsb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.genoSite = new BitSet[siteCount][];
            this.snps = new BiSNP[siteCount];
            for (Future<BlockVCF> blockVCFFuture : resultList) {
                BlockVCF block = blockVCFFuture.get();
                for (int j = 0; j < block.getActBlockSize(); j++) {
                    this.snps[block.getStartIndex() + j] = block.getSNPBlock()[j];
                    this.genoSite[block.getStartIndex() + j] = block.getGenoSiteBlock()[j];
                }
            }
            sb.setLength(0);
            sb.append("A total of ").append(this.getSiteNumber()).append(" SNPs are in ").append(infileS).append("\n");
            sb.append("Genotype table is successfully built");
            System.out.println(sb);
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        this.transposeSiteToTaxon();
        mafs = new float[this.getSiteNumber()];
        Arrays.fill(mafs, Float.MIN_VALUE);
    }

    private void transposeSiteToTaxon () {
        long start = System.nanoTime();
        genoTaxon = new BitSet[this.getTaxaNumber()][];
        for (int i = 0; i < genoTaxon.length; i++) {
            genoTaxon[i]=new BitSet[3];
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
            genoSite[i]=new BitSet[3];
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

    @Override
    public void swap(int index1, int index2) {
        BiSNP tempB;
        tempB = snps[index1];
        snps[index1] = snps[index2];
        snps[index2] = tempB;
        float tempM;
        tempM = mafs[index1];
        mafs[index1] = mafs[index2];
        mafs[index2] = tempM;
        BitSet[] tempS;
        tempS = genoSite[index1];
        genoSite[index1] = genoSite[index2];
        genoSite[index2] = tempS;
    }

    @Override
    public int compare(int index1, int index2) {
        return snps[index1].compareTo(snps[index2]);
    }



    /**
     * Merge the specified GenoGrid with this GenoGrid
     * Two genotypes should have the same number of taxa in the same order.
     * @param gt gt
     * @return Return a new GenoGrid
     */
    public GenoGrid mergeGenotypesBySite(GenoGrid gt) {
        assert this.getTaxaNumber()==gt.getTaxaNumber() : "check your GenoGrid";
        int snpCount = this.getSiteNumber()+gt.getSiteNumber();
        BitSet[][] bArray = new BitSet[snpCount][3];
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            System.arraycopy(this.genoSite[i], 0, bArray[cnt], 0, this.genoSite[0].length);
            cnt++;
        }
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            System.arraycopy(gt.genoSite[i], 0, bArray[cnt], 0, gt.genoSite[0].length);
            cnt++;
        }
        BiSNP[] nsnps = new BiSNP[snpCount];
        cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            nsnps[cnt] = this.snps[i];
            cnt++;
        }
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            nsnps[cnt] = gt.snps[i];
            cnt++;
        }
        return new GenoGrid(bArray, GenoGrid.Direction.BySite, this.taxa, nsnps);
    }

    /**
     * Merge the specified GenoGrid with this GenoGrid
     * <p> Two genotypes should have the same number of sites in the same order
     * @param gt gt
     * @return Return a new GenoGrid
     */
    public GenoGrid mergeGenotypesByTaxon(GenoGrid gt) {
        assert this.getSiteNumber()==gt.getSiteNumber() : "check your GenoGrid";
        if (this.getSiteNumber() != gt.getSiteNumber()) return null;
        int taxaCount = this.getTaxaNumber()+gt.getTaxaNumber();
        BitSet[][] bArray = new BitSet[taxaCount][3];
        String[] taxa = new String[taxaCount];
        int cnt = 0;
        for (int i = 0; i < this.getTaxaNumber(); i++) {
            System.arraycopy(this.genoTaxon[i], 0, bArray[cnt], 0, this.genoTaxon[0].length);
            taxa[cnt] = this.getTaxonName(i);
            cnt++;
        }
        for (int i = 0; i < gt.getTaxaNumber(); i++) {
            System.arraycopy(gt.genoTaxon[i], 0, bArray[cnt], 0, gt.genoTaxon[0].length);
            taxa[cnt] = gt.getTaxonName(i);
            cnt++;
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        for (int i = 0; i < nsnps.length; i++) {
            nsnps[i] = this.snps[i].replicateWithoutFeature();
        }
        return new GenoGrid(bArray, GenoGrid.Direction.ByTaxon, taxa, nsnps);
    }

    /**
     * Convert {@link GenoGrid} to {@link GenoRows}
     * @return
     */
    public GenoRows getConvertedGenotype() {
        SiteGeno[] geno = new SiteGeno[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            short chr = this.getChromosome(i);
            int pos = this.getPosition(i);
            char refBase = this.getReferenceAlleleBase(i);
            char altBase = this.getAlternativeAlleleBase(i);
            BitSet phase1 = this.genoSite[i][0];
            BitSet phase2 = this.genoSite[i][1];
            BitSet missingP = this.genoSite[i][2];
            geno[i] = new SiteGeno(chr, pos, refBase, altBase, null, phase1, phase2, missingP, this.getTaxaNumber());
        }
        return new GenoRows(geno, this.taxa);
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param siteIndices siteIndices
     * @return
     */
    public GenoGrid getSubsetGenotypeBySite(int[] siteIndices) {
        BitSet[][] bArray = new BitSet[siteIndices.length][3];
        for (int i = 0; i < siteIndices.length; i++) {
            System.arraycopy(this.genoSite[siteIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        BiSNP[] nsnps = new BiSNP[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            nsnps[i] = this.snps[siteIndices[i]];
        }
        return new GenoGrid(bArray, Direction.BySite, this.taxa, nsnps);
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param taxaIndices
     * @return
     */
    public GenoGrid getSubsetGenotypeByTaxon(int[] taxaIndices) {
        BitSet[][] bArray = new BitSet[taxaIndices.length][3];
        for (int i = 0; i < taxaIndices.length; i++) {
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            nTaxa[i] = this.getTaxonName(taxaIndices[i]);
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            nsnps[i] = this.snps[i].replicateWithoutFeature();
        }
        return new GenoGrid(bArray, Direction.ByTaxon, nTaxa, nsnps);
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param siteIndices
     */
    public void subsetsGenotypeBySite(int[] siteIndices) {
        BitSet[][] bArray = new BitSet[siteIndices.length][3];
        for (int i = 0; i < siteIndices.length; i++) {
            System.arraycopy(this.genoSite[siteIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        BiSNP[] nsnps = new BiSNP[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            nsnps[i] = this.snps[siteIndices[i]];
        }
        this.genoSite=bArray;
        this.snps=nsnps;
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param taxaIndices
     */
    public void subsetsGenotypeByTaxon(int[] taxaIndices) {
        BitSet[][] bArray = new BitSet[taxaIndices.length][3];
        for (int i = 0; i < taxaIndices.length; i++) {
            System.arraycopy(this.genoTaxon[taxaIndices[i]], 0, bArray[i], 0, bArray[0].length);
        }
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < taxaIndices.length; i++) {
            nTaxa[i] = this.getTaxonName(taxaIndices[i]);
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            nsnps[i] = this.snps[i].replicateWithoutFeature();
        }
        this.genoTaxon=bArray;
        this.taxa=nTaxa;
        this.snps=nsnps;
    }
}
