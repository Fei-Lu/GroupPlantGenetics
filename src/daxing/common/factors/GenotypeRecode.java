package daxing.common.factors;

public enum GenotypeRecode {
    REF("0/0",(byte)0),
    ALT("1/1", (byte)1),
    HETERO("0/1", (byte)2),
    MISSING("./.", (byte)9);

    String genotype;
    byte recodeByte;
    GenotypeRecode(String genotype, byte genotypeByte) {
        this.genotype=genotype;
        this.recodeByte =genotypeByte;
    }

    public byte getRecodeByte() {
        return recodeByte;
    }

    public String getGenotype() {
        return genotype;
    }

    /**
     * Allele byte from AlleleEncoder
     * missing ./. is N in AlleleEncoder
     * @param refAllele
     * @param altAllele
     * @param allele1
     * @param allele2
     * @return
     */
    public static GenotypeRecode newInstanceFromChar(char refAllele, char altAllele, char allele1, char allele2){
        if (allele1==allele2 && allele1 == refAllele) return GenotypeRecode.REF;
        if (allele1==allele2 && allele1 == altAllele) return GenotypeRecode.ALT;
        if (allele1==refAllele && allele2==altAllele) return GenotypeRecode.HETERO;
        if (allele2==refAllele && allele1==altAllele) return GenotypeRecode.HETERO;
        if (allele1==allele2 && allele1=='N') return GenotypeRecode.MISSING;
        return null;
    }

}
