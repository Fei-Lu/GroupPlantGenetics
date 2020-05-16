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
    String alt_SIFT;
    String gerp;
    String recombinationRate;

    public enum Region{
        UTR_5((byte)0), CDS((byte)1), UTR_3((byte)2);

        private byte region;

        Region(byte region){
            this.region=region;
        }
    }

    SNPAnnotation(short chr, int pos, char refBase, char altBase, String transcriptName, char majorBase,
                  String ancestral, double maf, double[] aaf, String daf, String[] dafs, Region region,
                  String variant_type, String alt_SIFT, String gerp, String recombinationRate){
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
        this.alt_SIFT=alt_SIFT;
        this.recombinationRate=recombinationRate;
        this.gerp=gerp;
        this.recombinationRate=recombinationRate;
    }

    public double getMaf() {
        return maf;
    }

    public String getAlt_SIFT() {
        return alt_SIFT;
    }

    public String getGerp() {
        return gerp;
    }

    public String getVariant_type() {
        return variant_type;
    }

    public boolean isDeleterious(){
        if(this.getAlt_SIFT().equals("NA")) return false;
        float sift=Float.parseFloat(this.getAlt_SIFT());
        if (sift > 0.05) return false;
        if(this.getGerp().equals("NA")) return false;
        float gerp=Float.parseFloat(this.getGerp());
        if (gerp < 1) return false;
        return true;
    }

    public boolean isNonSyn(){
        if (!this.getVariant_type().equals("NONSYNONYMOUS")) return false;
        return true;
    }

    public boolean isSyn(){
        if (!this.getVariant_type().equals("SYNONYMOUS")) return false;
        return true;
    }

    public boolean hasAncestral(){
        if(this.isReferenceAlleleTypeOf(AlleleType.Ancestral) || this.isAlternativeAlleleTypeOf(AlleleType.Ancestral)) return true;
        return false;
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
        if (this.isReferenceAlleleTypeOf(AlleleType.Ancestral)) return true;
        return false;
    }

    public boolean check(short chr, int pos){
        if (this.getChromosome()==chr && this.getPosition()==pos) return true;
        return false;
    }

    public int getPos(){
        return this.getPosition();
    }
}
