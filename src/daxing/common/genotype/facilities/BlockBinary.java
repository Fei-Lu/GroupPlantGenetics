package daxing.common.genotype.facilities;

import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.genot.GenotypeExport;
import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.Dyad;
import java.nio.ByteBuffer;
import java.util.BitSet;
import java.util.concurrent.Callable;

/**
 * Modify from GenoSiteBlockBinary, performance improvement, code simplification, easy customization
 * Multithreading reading facilities of GenoGrid
 */
public class BlockBinary implements Callable<BlockBinary> {
    public static final int blockSize = 4096;
    byte[][] lines;
    int startIndex;
    int actBlockSize;
    BitSet[][] genoSiteBlock = null;
    BiSNP[] snpBlock = null;

    public BlockBinary (byte[][] lines, int startIndex, int actBlockSize) {
        this.lines = lines;
        this.startIndex = startIndex;
        this.actBlockSize = actBlockSize;
    }

    public BiSNP[] getSNPBlock () {
        return this.snpBlock;
    }

    public int getActBlockSize() {
        return actBlockSize;
    }

    public BitSet[][] getGenoSiteBlock () {
        return this.genoSiteBlock;
    }

    public int getStartIndex () {
        return this.startIndex;
    }

    @Override
    public BlockBinary call() throws Exception {
        this.genoSiteBlock = new BitSet[this.actBlockSize][3];
        this.snpBlock = new BiSNP[this.actBlockSize];
        ByteBuffer bb = ByteBuffer.allocate(lines[0].length);
        for (int i = 0; i < this.actBlockSize; i++) {
            bb.put(lines[i]);
            Dyad<BiSNP, BitSet[]> d = buildFromBinaryLine(bb);
            snpBlock[i] = d.getFirstElement();
            genoSiteBlock[i] = d.getSecondElement();
        }
        lines = null;
        return this;
    }

    private Dyad<BiSNP, BitSet[]> buildFromBinaryLine (ByteBuffer bb) {
        //short chr, int pos, char refBase, char altBase, String info, BitSet phase1, BitSet phase2, BitSet missing, int taxaNumber
        bb.flip();
        short chr = bb.getShort();
        int pos = bb.getInt();
        byte geno = bb.get();
        char refBase = AlleleEncoder.getAlleleBase1FromGenotypeByte(geno);
        char altBase = AlleleEncoder.getAlleleBase2FromGenotypeByte(geno);
        byte refFeature = bb.get();
        byte altFeature = bb.get();
        int size = (bb.capacity()- GenotypeExport.getByteSizeOfSNPInBinary())/3;
        BitSet[] genoSite = new BitSet[3];
        for (int i = 0; i < genoSite.length; i++) {
            byte[] ba = new byte[size];
            bb.get(ba);
            genoSite[i] = BitSet.valueOf(ba);
        }
        BiSNP snp = new BiSNP(chr, pos, refBase, altBase, null);
        snp.setReferenceAlleleFeature(refFeature);
        snp.setAlternativeAlleleFeature(altFeature);
        bb.clear();
        return new Dyad<>(snp, genoSite);
    }
}
