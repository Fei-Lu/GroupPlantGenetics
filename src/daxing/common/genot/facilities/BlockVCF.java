package daxing.common.genot.facilities;

import pgl.infra.dna.snp.BiSNP;
import pgl.infra.utils.Dyad;
import pgl.infra.utils.PStringUtils;
import java.util.BitSet;
import java.util.List;
import java.util.concurrent.Callable;

/**
 * Modify from GenoSiteBlockVCF, performance improvement, code simplification, easy customization
 * Multithreading reading facilities of GenoGrid
 */
public class BlockVCF implements Callable<BlockVCF> {
    public static final int blockSize = 4096;
    List<String> lines;
    int startIndex;
    int actBlockSize;
    BitSet[][] genoSiteBlock = null;
    BiSNP[] snpBlock = null;

    public BlockVCF (List<String> lines, int startIndex) {
        this.lines = lines;
        this.startIndex = startIndex;
        this.actBlockSize = lines.size();
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
    public BlockVCF call() throws Exception {
        this.genoSiteBlock = new BitSet[this.actBlockSize][3];
        this.snpBlock = new BiSNP[this.actBlockSize];
        for (int i = 0; i < this.actBlockSize; i++) {
            Dyad<BiSNP, BitSet[]> d = buildFromVCFLine(lines.get(i));
            snpBlock[i] = d.getFirstElement();
            genoSiteBlock[i] = d.getSecondElement();
        }
        lines = null;
        return this;
    }

    public static Dyad<BiSNP, BitSet[]> buildFromVCFLine (String line) {
        List<String> l = PStringUtils.fastSplit(line);
        List<String> ll;
        String current;
        short chr = Short.parseShort(l.get(0));
        int pos = Integer.parseInt(l.get(1));
        char refBase = l.get(3).charAt(0);
        char altBase = l.get(4).charAt(0);
        BiSNP snp = new BiSNP(chr, pos, refBase, altBase, null);
        int taxaNumber = l.size()-9;
        BitSet[] genoSite = new BitSet[3];
        genoSite[0] = new BitSet(taxaNumber);
        genoSite[1] = new BitSet(taxaNumber);
        genoSite[2] = new BitSet(taxaNumber);
        byte[] values;
        for (int i = 0; i < taxaNumber; i++) {
            current = l.get(i+9);
            if (current.startsWith(".")) {
                genoSite[2].set(i);
                continue;
            }
            ll = PStringUtils.fastSplit(current, ":");
            values = ll.get(0).getBytes();
            if (values[0] == 49) {
                genoSite[0].set(i);
            }
            if (values[2] == 49) {
                genoSite[1].set(i);
            }
        }
        Dyad<BiSNP, BitSet[]> d = new Dyad<>(snp, genoSite);
        return d;
    }
}
