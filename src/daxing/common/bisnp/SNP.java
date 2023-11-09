package daxing.common.bisnp;

import daxing.common.chrrange.ChrPos;
import pgl.infra.dna.allele.Allele;
import pgl.infra.dna.allele.AlleleEncoder;
import pgl.infra.dna.allele.AlleleType;

/**
 * modified from BiSNP, extend from RefChrPos
 */
public class SNP extends ChrPos {

    /**
     * The reference allele
     */
    public Allele reference = null;
    /**
     * The alternative allele
     */
    public Allele alternative = null;

    public String snpID = null;
    /**
     * Annotation information of SNP
     */
    public String info = null;

    /**
     * Construct an object with all fields initialized
     * @param chr
     * @param pos
     * @param refBase
     * @param altBase
     * @param info
     */
    public SNP(String chr, int pos, char refBase, char altBase, String snpID, String info) {
        super(chr, pos);
        this.initializeRefAllele(refBase);
        this.initializeAltAllele(altBase);
        this.SetSNP_ID(snpID);
        this.setSNPInfo(info);

    }

    /**
     * Initialize the reference allele
     * @param refBase
     */
    private void initializeRefAllele(char refBase) {
        this.reference = new Allele (refBase);
        reference.setAlleleType(AlleleType.Reference);
    }

    /**
     * Initialize the alternative allele
     * @param altBase
     */
    private void initializeAltAllele(char altBase) {
        this.alternative = new Allele (altBase);
        alternative.setAlleleType(AlleleType.Alternative);
    }

    /**
     * Set the allele type of the reference allele
     * see {@link AlleleType}
     * @param at
     */
    public void setReferenceAlleleType (AlleleType at) {
        reference.setAlleleType(at);
    }

    /**
     * Set the allele feature of the reference allele
     * @param feature
     */
    public void setReferenceAlleleFeature (byte feature) {
        reference.setAlleleFeature(feature);
    }

    /**
     * Set the allele type of the alternative allele
     * see {@link AlleleType}
     * @param at
     */
    public void setAlternativeAlleleType (AlleleType at) {
        alternative.setAlleleType(at);
    }

    /**
     * Set the allele feature of the alternative allele
     * @param feature
     */
    public void setAlternativeAlleleFeature (byte feature) {
        this.alternative.setAlleleFeature(feature);
    }

    /**
     * Remove the allele type of the reference allele
     * see {@link AlleleType}
     * @param at
     */
    public void removeReferenceAlleleType (AlleleType at) {
        reference.removeAlleleType(at);
    }

    /**
     * Remove the allele type of the alternative allele
     * see {@link AlleleType}
     * @param at
     */
    public void removeAlternativeAlleleType (AlleleType at) {
        alternative.removeAlleleType(at);
    }

    /**
     * Return if the reference allele is a certain allele type
     * see {@link AlleleType}
     * @param at
     * @return
     */
    public boolean isReferenceAlleleTypeOf (AlleleType at) {
        return reference.isAlleleTypeOf(at);
    }

    /**
     * Return if the alternative allele is a certain allele type
     * see {@link AlleleType}
     * @param at
     * @return
     */
    public boolean isAlternativeAlleleTypeOf (AlleleType at) {
        return alternative.isAlleleTypeOf(at);
    }

    /**
     * Return the byte value of the reference allele
     * @return
     */
    public byte getReferenceAlleleByte() {
        return reference.getAlleleCoding();
    }

    /**
     * Return the base of the reference allele
     * @return
     */
    public char getReferenceAlleleBase() {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getReferenceAlleleByte());
    }
    /**
     * Return the byte value of the alternative allele
     * @return
     */
    public byte getAlternativeAlleleByte() {
        return alternative.getAlleleCoding();
    }

    /**
     * Return the base of the alternative allele
     * @return
     */
    public char getAlternativeAlleleBase () {
        return AlleleEncoder.getAlleleBaseFromCoding(this.getAlternativeAlleleByte());
    }

    /**
     * Return the feature of the reference allele, see {@link AlleleType}
     * @return
     */
    public byte getReferenceAlleleFeature () {
        return this.reference.getAlleleFeature();
    }

    /**
     * Return the feature of the alternative allele, see {@link AlleleType}
     * @return
     */
    public byte getAlternativeAlleleFeature () {
        return this.alternative.getAlleleFeature();
    }

    /**
     *
     * @return the snpID of current snp
     */
    public String getSnpID(){
        return this.snpID;
    }

    /**
     * set the snpID of current snp
     * @param snpID
     */
    public void SetSNP_ID(String snpID){
        this.snpID=snpID;
    }

    /**
     * Return the annotation information of the SNP
     * @return
     */
    public String getSNPInfo () {
        return this.info;
    }

    /**
     * Set the information of the current SNP
     * @param info
     */
    public void setSNPInfo (String info) {
        this.info = info;
    }

    /**
     * Make a copy of the current object, without allele feature replicated (initialized to 0)
     * @return
     */
    public SNP replicateWithoutFeature() {
        return new SNP(this.getChr(), this.getPos(), this.getReferenceAlleleBase(), this.getAlternativeAlleleBase(),
                this.getSnpID(), this.getSNPInfo());
    }

    /**
     * Make a copy of the current object, with allele feature replicated
     * @return
     */
    public SNP replicateWithFeature() {
        SNP s = this.replicateWithoutFeature();
        s.setReferenceAlleleFeature(this.getReferenceAlleleFeature());
        s.setAlternativeAlleleFeature(this.getAlternativeAlleleFeature());
        return s;
    }
}
