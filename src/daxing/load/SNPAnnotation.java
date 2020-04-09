package daxing.load;

import pgl.infra.dna.allele.AlleleType;
import pgl.infra.dna.snp.BiSNP;

public class SNPAnnotation {

    BiSNP snp;
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
        BiSNP biSNP=new BiSNP(chr, pos, refBase, altBase, transcriptName);
        if (biSNP.getReferenceAlleleBase()==majorBase){
            biSNP.setReferenceAlleleType(AlleleType.Major);
            biSNP.setAlternativeAlleleType(AlleleType.Minor);
        }else {
            biSNP.setAlternativeAlleleType(AlleleType.Major);
            biSNP.setReferenceAlleleType(AlleleType.Minor);
        }
        if (!ancestral.equals("NA")){
            char ancestralAllele=ancestral.charAt(0);
            if (biSNP.getReferenceAlleleBase()==ancestralAllele){
                biSNP.setReferenceAlleleType(AlleleType.Ancestral);
                biSNP.setAlternativeAlleleType(AlleleType.Derived);
            }else if (biSNP.getAlternativeAlleleBase()==ancestralAllele){
                biSNP.setAlternativeAlleleType(AlleleType.Ancestral);
                biSNP.setReferenceAlleleType(AlleleType.Derived);
            }
        }
        this.snp=biSNP;
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

    public BiSNP getSnp() {
        return snp;
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
        BiSNP snp=this.getSnp();
        if(snp.isReferenceAlleleTypeOf(AlleleType.Ancestral) || snp.isAlternativeAlleleTypeOf(AlleleType.Ancestral)) return true;
        return false;
    }

    public boolean isRefAlleleAncestral(){
        if (!hasAncestral()){
            System.out.println("error, program quit");
            System.exit(1);
        }
        if (this.snp.isReferenceAlleleTypeOf(AlleleType.Ancestral)) return true;
        return false;
    }

    public boolean check(short chr, int pos){
        BiSNP snp=this.getSnp();
        if (snp.getChromosome()==chr && snp.getPosition()==pos) return true;
        return false;
    }

    public int getPos(){
        return this.snp.getPosition();
    }
}
