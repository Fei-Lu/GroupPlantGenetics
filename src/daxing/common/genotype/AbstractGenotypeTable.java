package daxing.common.genotype;

import daxing.common.utiles.IOTool;
import pgl.infra.dna.genot.GenoIOFormat;
import pgl.infra.dna.genot.GenotypeTable;
import pgl.infra.dna.genot.VCFUtils;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.nio.ByteBuffer;

/**
 * Partial method implementation of {@link GenotypeTable}, providing functions of write, getByteSizeOfSiteInBinary etc.
 * @author Daxing Xu
 */
public abstract class AbstractGenotypeTable implements GenotypeTable {

    public void write(String outFile, GenoIOFormat genoIOFormat){
        switch (genoIOFormat){
            case VCF:
            case VCF_GZ:
                toVCF(outFile);
                break;
            case Binary:
            case Binary_GZ:
                toBinary(outFile);
                break;
            case HDF5:
                throw new UnsupportedOperationException("Not supported yet.");
        }
        System.out.println("Genotype table exported to "+ outFile);
    }

    protected void toVCF (String outFile) {
        BufferedWriter bw = IOTool.getWriter(outFile);
        try {
            bw.write(VCFUtils.getVCFAnnotation());
            bw.newLine();
            bw.write(VCFUtils.getVCFHeader(this.getTaxaNames()));
            bw.newLine();
            for (int i = 0; i < this.getSiteNumber(); i++) {
                bw.write(this.getUnphasedVCFRecord(i));
                bw.newLine();
            }
            bw.flush();
            bw.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    protected void toBinary (String outFile) {
        try {
            DataOutputStream dos= IOTool.getBinaryWriter(outFile);
            int siteNumber = this.getSiteNumber();
            int taxaNumber = this.getTaxaNumber();
            dos.writeInt(siteNumber);
            dos.writeInt(taxaNumber);
            for (int i = 0; i < taxaNumber; i++) {
                dos.writeUTF(this.getTaxonName(i));
            }
            int byteSize = this.getByteSizeOfSiteInBinary();
            ByteBuffer bb = ByteBuffer.allocate(byteSize);
            for (int i = 0; i < siteNumber; i++) {
                bb = this.getBinaryOutput(i, bb);
                bb.flip();
                dos.write(bb.array());
                bb.clear();
            }
            dos.close();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     * Return the number of byte for a {@link pgl.infra.dna.snp.BiSNP}
     * @return
     */
    protected static int getByteSizeOfSNPInBinary () {
        //short chr, int pos, genotype (ref+alt), byte feature *2
        return 9;
    }

    /**
     * Return the number of byte for genotype of a site (chr, pos, ref+alt, alternative at phase 1 and phase 2, and missing)
     * @return
     */
    protected int getByteSizeOfSiteInBinary () {
        int taxaNumber= this.getTaxaNumber();
        //short chr, int pos, genotype (ref+alt), byte feature *2, BitSet phase1, BitSet phase2, BitSet missing
        int n = getByteSizeOfSNPInBinary();
        if (taxaNumber%64 == 0) {
            n += (taxaNumber/64)*8*3;
        }
        else  {
            n += (taxaNumber/64+1)*8*3;
        }
        return n;
    }
}
