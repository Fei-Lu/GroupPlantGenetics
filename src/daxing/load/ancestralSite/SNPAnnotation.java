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

    public enum Region{
        UTR_5((byte)0), CDS((byte)1), UTR_3((byte)2);

        private byte region;

        Region(byte region){
            this.region=region;
        }
    }

    SNPAnnotation(short chr, int pos, char refBase, char altBase, String transcriptName, char majorBase,
                  String ancestral, double maf, double[] aaf, String daf, String[] dafs, Region region,
                  String variant_type, String derived_SIFT, String gerp, String recombinationRate){
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

    public boolean isDeleterious(){
        if(!isNonSyn()) return false;
        if (!hasAncestral()) return false;
        if(this.getDerived_SIFT().equals("NA")) return false;
        double sift=Double.parseDouble(this.getDerived_SIFT());
        if (sift > 0.05) return false;
        if(this.getGerp().equals("NA")) return false;
        double gerp=Double.parseDouble(this.getGerp());
        return !(gerp < 1);
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
