package daxing.common.genotype;

import daxing.common.utiles.IOTool;
import daxing.common.genotype.facilities.SiteGeno;
import pgl.PGLConstraints;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.genot.*;
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
import java.util.concurrent.*;
import java.util.stream.IntStream;

/**
 * Modify from GenotypeRows, performance improvement, code simplification, easy customization
 * The bit version implementation of a genotype table, see {@link GenotypeTable}
 * <p>
 * The class implements bitset genotype by site, which means it is smaller than {@link GenoGrid} in memory,
 * but slower in taxon based calculation, e.g. matrix of genetic divergence between taxa
 * <p>
 *  * Supports only bi-allelic SNPs, 3rd+ allele will be ignored. Allele depth is ignored.
 *
 * @author feilu
 */
public class GenoRows extends AbstractGenotypeTable {

    /**
     * The taxa in the genotype table
     */
    String[] taxa = null;
    /**
     * SNP, alleles and their presentation in individuals
     */
    SiteGeno[] geno = null;

    /**
     * Construct an object
     */
    public GenoRows() {

    }

    /**
     * Construct an object by reading a file
     * @param infileS file
     * @param format format
     */
    public GenoRows(String infileS, GenoIOFormat format) {
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
    }

    /**
     * Construct an object
     * @param geno siteGeno
     * @param taxa taxa
     */
    public GenoRows (SiteGeno[] geno, String[] taxa) {
        this.taxa = taxa;
        this.geno = geno;
    }

    @Override
    public int getTaxaNumber() {
        return taxa.length;
    }

    @Override
    public String getTaxonName(int taxonIndex) {
        return taxa[taxonIndex];
    }

    @Override
    public String[] getTaxaNames () {
        return this.taxa;
    }

    @Override
    public int getTaxonIndex(String taxon) {
        return Arrays.binarySearch(taxa, taxon);
    }

    @Override
    public int getSiteNumber () {
        return geno.length;
    }

    @Override
    public short getChromosome(int siteIndex) {
        return geno[siteIndex].getChromosome();
    }

    @Override
    public int getPosition(int siteIndex) {
        return geno[siteIndex].getPosition();
    }

    @Override
    public void sortBySite() {
        Arrays.sort(geno);
    }

    @Override
    public void sortByTaxa() {
        int[] indices = PArrayUtils.getIndicesByAscendingValue(this.taxa);
        Arrays.sort(this.taxa);
        Arrays.stream(this.geno).parallel().forEach(geno->geno.sortByTaxa(indices));
    }

    @Override
    public byte getGenotypeByte (int siteIndex, int taxonIndex) {
        return geno[siteIndex].getGenotypeByte(taxonIndex);
    }

    @Override
    public boolean isHeterozygous(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isHeterozygous(taxonIndex);
    }

    @Override
    public boolean isHomozygous(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isHomozygous(taxonIndex);
    }

    @Override
    public boolean isMissing(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isMissing(taxonIndex);
    }

    @Override
    public boolean isPhase1Alternative(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase1Alternative(taxonIndex);
    }

    @Override
    public boolean isPhase2Alternative(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase2Alternative(taxonIndex);
    }

    @Override
    public boolean isPhase1Reference(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase1Reference(taxonIndex);
    }

    @Override
    public boolean isPhase2Reference(int siteIndex, int taxonIndex) {
        return geno[siteIndex].isPhase2Reference(taxonIndex);
    }

    @Override
    public boolean isAlternativeAlleleTypeOf(AlleleType at, int siteIndex) {
        return geno[siteIndex].isAlternativeAlleleTypeOf(at);
    }

    @Override
    public boolean isReferenceAlleleTypeOf(AlleleType at, int siteIndex) {
        return geno[siteIndex].isReferenceAlleleTypeOf(at);
    }

    @Override
    public int getSiteIndex(short chromosome, int position) {
        ChrPos query = new ChrPos (chromosome, position);
        return Arrays.binarySearch(geno, query);
    }

    @Override
    public int getMissingNumberBySite(int siteIndex) {
        return geno[siteIndex].getMissingNumber();
    }

    @Override
    public int getMissingNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isMissing(taxonIndex)) cnt++;
        }
        return cnt;
    }

    @Override
    public int getNonMissingNumberBySite(int siteIndex) {
        return this.getTaxaNumber()-this.getMissingNumberBySite(siteIndex);
    }

    @Override
    public int getNonMissingNumberByTaxon(int taxonIndex) {
        return this.getSiteNumber()-this.getMissingNumberByTaxon(taxonIndex);
    }

    @Override
    public int getHomozygoteNumberBySite(int siteIndex) {
        return geno[siteIndex].getHomozygoteNumber();
    }

    @Override
    public int getHomozygoteNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isHomozygous(taxonIndex)) cnt++;
        }
        return cnt;
    }

    @Override
    public int getHeterozygoteNumberBySite(int siteIndex) {
        return geno[siteIndex].getHeterozygoteNumber();
    }

    @Override
    public int getHeterozygoteNumberByTaxon(int taxonIndex) {
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            if (geno[i].isHeterozygous(taxonIndex)) cnt++;
        }
        return cnt;
    }

    @Override
    public int getAlternativeAlleleNumberBySite(int siteIndex) {
        return geno[siteIndex].getAlternativeAlleleNumber();
    }

    @Override
    public int getAlternativeAlleleOccurrenceBySite(int siteIndex) {
        BitSet bs = geno[siteIndex].getPhase1();
        bs.or(geno[siteIndex].getPhase2());
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
        return geno[siteIndex].getMinorAlleleByte();
    }

    @Override
    public char getMinorAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromByte(this.getMinorAlleleByte(siteIndex));
    }

    @Override
    public float getMinorAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getMinorAlleleFrequency();
    }

    @Override
    public byte getMajorAlleleByte(int siteIndex) {
        return geno[siteIndex].getMinorAlleleByte();
    }

    @Override
    public char getMajorAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromByte(this.getMinorAlleleByte(siteIndex));
    }

    @Override
    public float getMajorAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getMajorAlleleFrequency();
    }

    @Override
    public byte getReferenceAlleleByte(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleByte();
    }

    @Override
    public char getReferenceAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromByte(this.getReferenceAlleleByte(siteIndex));
    }

    @Override
    public float getReferenceAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleFrequency();
    }

    @Override
    public byte getAlternativeAlleleByte(int siteIndex) {
        return geno[siteIndex].getAlternativeAlleleByte();
    }

    @Override
    public char getAlternativeAlleleBase(int siteIndex) {
        return AlleleEncoder.getAlleleBaseFromByte(this.getAlternativeAlleleByte(siteIndex));
    }

    @Override
    public float getAlternativeAlleleFrequency(int siteIndex) {
        return geno[siteIndex].getReferenceAlleleFrequency();
    }

    @Override
    public String getUnphasedVCFRecord(int siteIndex) {
        return geno[siteIndex].getUnphasedVCFOutput();
    }

    @Override
    public ByteBuffer getBinaryOutput (int siteIndex, ByteBuffer bb) {
        return geno[siteIndex].getBinaryOutput(bb);
    }

    @Override
    public void setAlternativeAlleleType(AlleleType at, int siteIndex) {
        geno[siteIndex].setAlternativeAlleleType(at);
    }

    @Override
    public void setReferenceAlleleType(AlleleType at, int siteIndex) {
        geno[siteIndex].setReferenceAlleleType(at);
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

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2) {
        return this.getIBSDistance(taxonIndex1, taxonIndex2, 0, this.getSiteNumber());
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int startSiteIndex, int endSiteIndex) {
        int cnt = 0;
        double dxy = 0;
        for (int i = startSiteIndex; i < endSiteIndex; i++) {
            float dxySite = this.getDxySite(taxonIndex1, taxonIndex2, i);
            if (Float.isNaN(dxySite)) continue;
            dxy+=dxySite;
            cnt++;
        }
        return (float)(dxy/cnt);
    }

    @Override
    public float getIBSDistance(int taxonIndex1, int taxonIndex2, int[] siteIndices) {
        int cnt = 0;
        double dxy = 0;
        for (int siteIndex : siteIndices) {
            float dxySite = this.getDxySite(taxonIndex1, taxonIndex2, siteIndex);
            if (Float.isNaN(dxySite)) continue;
            dxy += dxySite;
            cnt++;
        }
        return (float)(dxy/cnt);
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

    /**
     * Return the dxy bwteen two taxa at a specific site
     * @param taxonIndex1 taxonIndex1
     * @param taxonIndex2 taxonIndex2
     * @param siteIndex siteIndex
     * @return Float.NaN if no shared non-missing site available
     */
    private float getDxySite (int taxonIndex1, int taxonIndex2, int siteIndex) {
        if (this.isMissing(siteIndex, taxonIndex1)) return Float.NaN;
        if (this.isMissing(siteIndex, taxonIndex2)) return Float.NaN;
        boolean isTaxon1Phase1Alt=this.geno[siteIndex].isPhase1Alternative(taxonIndex1);
        boolean isTaxon1Phase2Alt=this.geno[siteIndex].isPhase2Alternative(taxonIndex1);
        boolean isTaxon2Phase1Alt=this.geno[siteIndex].isPhase1Alternative(taxonIndex2);
        boolean isTaxon2Phase2Alt=this.geno[siteIndex].isPhase2Alternative(taxonIndex2);
        if (isTaxon1Phase1Alt && isTaxon1Phase2Alt && isTaxon2Phase1Alt && isTaxon2Phase2Alt){
            return 0f;
        }else if (!isTaxon1Phase1Alt && !isTaxon1Phase2Alt && !isTaxon2Phase1Alt && !isTaxon2Phase2Alt){
            return 0f;
        }else if (isTaxon1Phase1Alt!=isTaxon1Phase2Alt || isTaxon2Phase1Alt!=isTaxon2Phase2Alt){
            return 0.5f;
        }else {
            return 1f;
        }
    }

    @Override
    public GenotypeTable getSubGenotypeTableBySite(int[] siteIndices) {
        SiteGeno[] ge = new SiteGeno[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            ge[i] = geno[siteIndices[i]];
        }
        return new GenoRows(ge, taxa);
    }

    @Override
    public GenotypeTable getSubGenotypeTableByTaxa(int[] taxaIndices) {
        String[] nTaxa = new String[taxaIndices.length];
        for (int i = 0; i < nTaxa.length; i++) {
            nTaxa[i] = taxa[taxaIndices[i]];
        }
        SiteGeno[] nGeno = new SiteGeno[this.getSiteNumber()];
        IntStream.range(0, this.getSiteNumber()).parallel().forEach(index -> nGeno[index] = geno[index].getSubGenotypeByTaxa(taxaIndices));
        return new GenoRows(nGeno, nTaxa);
    }

    /**
     * Build an object from a binary genotype file
     * @param infileS infileS
     */
    private void buildFromBinary (String infileS) {
        try{
            DataInputStream dis = IOTool.getBinaryReader(infileS);
            int siteNumber = dis.readInt();
            short taxaNumber = (short)dis.readInt();
            this.taxa = new String[taxaNumber];
            for (int i = 0; i < taxaNumber; i++) {
                this.taxa[i] = dis.readUTF();
            }
            ExecutorService pool = Executors.newFixedThreadPool(PGLConstraints.parallelLevel);
            List<Future<BlockBinarySiteGeno>> resultList = new ArrayList<>();
            byte[][] input = new byte[BlockBinarySiteGeno.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
            int siteCount = 0;
            int startIndex = 0;
            int actBlockSize = 0;
            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < siteNumber; i++) {
                dis.read(input[actBlockSize]);
                siteCount++;
                actBlockSize++;
                if (siteCount% BlockBinarySiteGeno.blockSize == 0) {
                    BlockBinarySiteGeno sgb = new BlockBinarySiteGeno(input, startIndex, actBlockSize, taxaNumber);
                    Future<BlockBinarySiteGeno> result = pool.submit(sgb);
                    resultList.add(result);
                    startIndex+=actBlockSize;
                    actBlockSize = 0;
                    input = new byte[BlockBinarySiteGeno.blockSize][GenotypeExport.getByteSizeOfSiteInBinary(taxaNumber)];
                }
                if (siteCount%1000000 == 0) {
                    sb.setLength(0);
                    sb.append("Read in ").append(siteCount).append(" SNPs from ").append(infileS);
                    System.out.println(sb);
                }
            }
            dis.close();
            if (actBlockSize != 0) {
                BlockBinarySiteGeno sgb = new BlockBinarySiteGeno(input, startIndex, actBlockSize, taxaNumber);
                Future<BlockBinarySiteGeno> result = pool.submit(sgb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.geno = new SiteGeno[siteCount];
            for (Future<BlockBinarySiteGeno> blockBinarySiteGenoFuture : resultList) {
                BlockBinarySiteGeno block = blockBinarySiteGenoFuture.get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    geno[block.getStartIndex() + j] = block.getSiteGenotypes()[j];
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
    }

    /**
     * Class for parallel reading in a binary genotype file
     */
    static class BlockBinarySiteGeno implements Callable<BlockBinarySiteGeno> {
        public static final int blockSize = 4096;
        byte[][] lines;
        int startIndex;
        int actBlockSize;
        short taxaNumber;
        SiteGeno[] sgbArray = null;

        public BlockBinarySiteGeno(byte[][] lines, int startIndex, int actBlockSize, short taxaNumber) {
            this.lines = lines;
            this.startIndex = startIndex;
            this.actBlockSize = actBlockSize;
            this.taxaNumber = taxaNumber;
        }

        public int getStartIndex () {
            return this.startIndex;
        }

        public SiteGeno[] getSiteGenotypes () {
            return this.sgbArray;
        }

        @Override
        public BlockBinarySiteGeno call() throws Exception {
            this.sgbArray = new SiteGeno[this.actBlockSize];
            ByteBuffer bb = ByteBuffer.allocate(lines[0].length);
            for (int i = 0; i < this.actBlockSize; i++) {
                bb.put(lines[i]);
                sgbArray[i] = SiteGeno.buildFromBinaryLine(bb, taxaNumber);
            }
            lines = null;
            return this;
        }
    }

    /**
     * Build an object from a VCF file
     * @param infileS infileS
     */
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
            List<Future<BlockSiteGeno>> resultList = new ArrayList<>();
            int siteCount = 0;
            int startIndex = 0;
            ArrayList lines = new ArrayList ();
            StringBuilder sb = new StringBuilder();
            while ((temp = br.readLine()) != null) {
                lines.add(temp);
                if (lines.size()% BlockSiteGeno.blockSize == 0) {
                    BlockSiteGeno sgb = new BlockSiteGeno(lines, startIndex);
                    Future<BlockSiteGeno> result = pool.submit(sgb);
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
                BlockSiteGeno sgb = new BlockSiteGeno(lines, startIndex);
                Future<BlockSiteGeno> result = pool.submit(sgb);
                resultList.add(result);
            }
            pool.shutdown();
            pool.awaitTermination(Long.MAX_VALUE, TimeUnit.MICROSECONDS);
            this.geno = new SiteGeno[siteCount];
            for (Future<BlockSiteGeno> blockSiteGenoFuture : resultList) {
                BlockSiteGeno block = blockSiteGenoFuture.get();
                for (int j = 0; j < block.actBlockSize; j++) {
                    geno[block.getStartIndex() + j] = block.getSiteGenotypes()[j];
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
    }

    /**
     * Class for parallel reading in VCF
     */
    static class BlockSiteGeno implements Callable<BlockSiteGeno> {
        public static final int blockSize = 4096;
        List<String> lines;
        int startIndex;
        int actBlockSize;
        SiteGeno[] sgbArray = null;

        public BlockSiteGeno(List<String> lines, int startIndex) {
            this.lines = lines;
            this.startIndex = startIndex;
            this.actBlockSize = lines.size();
        }

        public int getStartIndex () {
            return this.startIndex;
        }

        public SiteGeno[] getSiteGenotypes () {
            return this.sgbArray;
        }

        @Override
        public BlockSiteGeno call() throws Exception {
            this.sgbArray = new SiteGeno[this.actBlockSize];
            for (int i = 0; i < this.actBlockSize; i++) {
                sgbArray[i] = SiteGeno.buildFromVCFLine(lines.get(i));
            }
            lines = null;
            return this;
        }
    }

    /**
     * Merge the specified GenoRows with this GenoRows
     * <p> Two genotypes should have the same number of taxa in the same order.
     * @param gt gt
     * @return Return a new GenoGrid
     */
    public GenoRows mergeGenotypesBySite(GenoRows gt) {
        long start=System.nanoTime();
        assert this.getTaxaNumber()== gt.getTaxaNumber() : "check your GenoRows";
        int snpCount = this.getSiteNumber()+gt.getSiteNumber();
        SiteGeno[] geno = new SiteGeno[snpCount];
        int cnt = 0;
        for (int i = 0; i < this.getSiteNumber(); i++) {
            geno[cnt] = this.geno[i];
            cnt++;
        }
        for (int i = 0; i < gt.getSiteNumber(); i++) {
            geno[cnt] = gt.geno[i];
            cnt++;
        }
        System.out.println("merge spend "+ Benchmark.getTimeSpanSeconds(start)+ " seconds");
        return new GenoRows(geno, this.taxa);
    }

    /**
     * Convert {@link GenoGrid} to {@link GenoRows}
     * @return GenoGrid
     */
    public GenoGrid getConvertedGenotype() {
        BitSet[][] bArray = new BitSet[this.getSiteNumber()][3];
        for (int i = 0; i < this.getSiteNumber(); i++) {
            bArray[i][0] = this.geno[i].getPhase1();
            bArray[i][1] = this.geno[i].getPhase2();
            bArray[i][2] = this.geno[i].getMissing();
        }
        BiSNP[] nsnps = new BiSNP[this.getSiteNumber()];
        System.arraycopy(this.geno, 0, nsnps, 0, nsnps.length);
        return new GenoGrid(bArray, GenoGrid.Direction.BySite, this.taxa, nsnps);
    }

    /**
     * Return a new subset of original genotype. The original genotype is unchanged.
     * @param siteIndices siteIndices
     * @return GenoRows
     */
    public GenoRows getSubsetGenotypeBySite(int[] siteIndices) {
        SiteGeno[] ge = new SiteGeno[siteIndices.length];
        for (int i = 0; i < siteIndices.length; i++) {
            ge[i] = this.geno[siteIndices[i]];
        }
        return new GenoRows(ge, this.taxa);
    }

    /**
     * Subsets the original genotype. The original genotype is changed.
     * @param siteIndices siteIndices
     */
    public void subsetsGenotypeBySite(int[] siteIndices) {
        getSubsetGenotypeBySite(siteIndices);
    }
}
