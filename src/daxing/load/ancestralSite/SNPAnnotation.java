package daxing.load.ancestralSite;

import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.snp.BiSNP;

public class SNPAnnotation extends BiSNP{

    double maf;
    double[] aaf;
    String daf;
    String[] dafs;
    Region region;
    String variant_type;
    String derived_SIFT;
    String gerp;
    String recombinationRate;

    String effect_snpEff;
    String impact_snpEff;
    String effect_vep;
    String impact_vep;

    public enum Region{
        UTR_5((byte)0), CDS((byte)1), UTR_3((byte)2), Intron((byte)3);

        private byte region;

        Region(byte region){
            this.region=region;
        }
    }

    public enum MethodCallDeleterious{
        SIFT, GERP, SIFT_GERP, STOPGAIN_STARTLOST_SIFT, STOPGAIN_STARTLOST_VEP, STOPGAIN_STARTLOST_SNPEFF, HIGH_VEP,
        HIGH_SNPEFF;

        // SIFT < 0.05
        // GERP > 1

    }

    SNPAnnotation(short chr, int pos, char refBase, char altBase, String transcriptName, char majorBase,
                  String ancestral, double maf, double[] aaf, String daf, String[] dafs, Region region,
                  String variant_type, String derived_SIFT, String gerp, String recombinationRate,
                  String effect_snpEff, String impact_snpEff, String effect_vep, String impact_vep){
        super(chr, pos, refBase, altBase, transcriptName);
        if (this.getReferenceAlleleBase()==majorBase){
            this.setReferenceAlleleType(AlleleType.Major);
            this.setAlternativeAlleleType(AlleleType.Minor);
        }else {
            this.setAlternativeAlleleType(AlleleType.Major);
            this.setReferenceAlleleType(AlleleType.Minor);
        }
        if (!ancestral.equals("NA")){
            char ancestralAllele=ancestral.charAt(0);
            if (this.getReferenceAlleleBase()==ancestralAllele){
                this.setReferenceAlleleType(AlleleType.Ancestral);
                this.setAlternativeAlleleType(AlleleType.Derived);
            }else if (this.getAlternativeAlleleBase()==ancestralAllele){
                this.setAlternativeAlleleType(AlleleType.Ancestral);
                this.setReferenceAlleleType(AlleleType.Derived);
            }
        }
        this.maf=maf;
        this.aaf=aaf;
        this.daf=daf;
        this.dafs=dafs;
        this.region=region;
        this.variant_type=variant_type;
        this.derived_SIFT =derived_SIFT;
        this.recombinationRate=recombinationRate;
        this.gerp=gerp;
        this.recombinationRate=recombinationRate;
        this.effect_vep=effect_vep;
        this.impact_vep=impact_vep;
        this.effect_snpEff=effect_snpEff;
        this.impact_snpEff=impact_snpEff;
    }

    public double getMaf() {
        return maf;
    }

    public String getDerived_SIFT() {
        return derived_SIFT;
    }

    public String getGerp() {
        return gerp;
    }

    public String getVariant_type() {
        return variant_type;
    }

    public double getDaf() {
        return daf.equals("NA") ? -1 : Double.parseDouble(daf);
    }

    public String getEffect_snpEff() {
        return effect_snpEff;
    }

    public String getEffect_vep() {
        return effect_vep;
    }

    public String getImpact_snpEff() {
        return impact_snpEff;
    }

    public String getImpact_vep() {
        return impact_vep;
    }

    public boolean isDeleterious(){
        if(!isNonSyn()) return false;
        if (!hasAncestral()) return false;
        if(this.getDerived_SIFT().equals("NA")) return false;
        double sift=Double.parseDouble(this.getDerived_SIFT());
        if (sift >= 0.05) return false;
        if(this.getGerp().equals("NA")) return false;
        double gerp=Double.parseDouble(this.getGerp());
        return !(gerp <= 1);
//        if (this.getVariant_type().equals("STOP-GAIN") || this.getVariant_type().equals("START-LOST")) return true;
//        return false;
    }

    public boolean isDeleterious(MethodCallDeleterious methodCallDeleterious){
        if (!hasAncestral()) return false;
        double sift, gerp;
        switch (methodCallDeleterious){
            case SIFT:
                if(!isNonSyn()) return false;
                if(this.getDerived_SIFT().equals("NA")) return false;
                sift=Double.parseDouble(this.getDerived_SIFT());
                if (sift >= 0.05) return false;
                return true;
            case GERP:
                if(!isNonSyn()) return false;
                if(this.getGerp().equals("NA")) return false;
                gerp=Double.parseDouble(this.getGerp());
                return !(gerp < 2.14);
            case SIFT_GERP:
                if(!isNonSyn()) return false;
                if(this.getDerived_SIFT().equals("NA")) return false;
                sift=Double.parseDouble(this.getDerived_SIFT());
                if (sift >= 0.05) return false;
                if(this.getGerp().equals("NA")) return false;
                gerp=Double.parseDouble(this.getGerp());
                return !(gerp < 2.14);
            case STOPGAIN_STARTLOST_SIFT:
                if (this.getVariant_type().equals("STOP-GAIN") || this.getVariant_type().equals("START-LOST")) return true;
                return false;
            case STOPGAIN_STARTLOST_VEP:
                if (this.getEffect_vep().equals("stop_gained") || this.getEffect_vep().equals("start_lost")) return true;
                return false;
            case STOPGAIN_STARTLOST_SNPEFF:
                if (this.getEffect_snpEff().equals("stop_gained") || this.getEffect_snpEff().equals("start_lost")) return true;
                return false;
            case HIGH_VEP:
                return this.getImpact_vep().equals("HIGH");
            case HIGH_SNPEFF:
                return this.getImpact_snpEff().equals("HIGH");
        }
        return false;
    }

    public boolean isNonSyn(){
        return this.getVariant_type().equals("NONSYNONYMOUS");
    }

    public boolean isSyn(){
        return this.getVariant_type().equals("SYNONYMOUS");
    }

    public boolean hasAncestral(){
        return this.isReferenceAlleleTypeOf(AlleleType.Ancestral) || this.isAlternativeAlleleTypeOf(AlleleType.Ancestral);
    }

    public char getAncestral(){
        if (this.isRefAlleleAncestral()){
            return this.reference.getAlleleBase();
        }else {
            return this.alternative.getAlleleBase();
        }
    }

    public boolean isRefAlleleAncestral(){
        if (!hasAncestral()){
            System.out.println("error, program quit");
            System.exit(1);
        }
        return this.isReferenceAlleleTypeOf(AlleleType.Ancestral);
    }

    public boolean check(short chr, int pos){
        return this.getChromosome() == chr && this.getPosition() == pos;
    }

    public int getPos(){
        return this.getPosition();
    }
}
