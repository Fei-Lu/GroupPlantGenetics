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

    String gerp16way;
    int aaPos;
    int aaSize;
    double list_s2;
    String phylopP_refMask;
    double alleleAgeJ;
    double alleleAgeM;
    double alleleAgeR;
    String derived_PolyPhen2_prediction;
    String derived_PolyPhen2_class;
    double derived_PolyPhen2_prob;


    public enum Region{
        UTR_5((byte)0), CDS((byte)1), UTR_3((byte)2), Intron((byte)3);

        private byte region;

        Region(byte region){
            this.region=region;
        }
    }

    public enum MethodCallDeleterious{
        SIFT, GERP, SIFT_GERP, STOPGAIN_STARTLOST_SIFT, STOPGAIN_STARTLOST_VEP, STOPGAIN_STARTLOST_SNPEFF, HIGH_VEP,
        HIGH_SNPEFF, PHYLOP, LIST_S2, PPH2, GERP_PPH2

        // SIFT < 0.05
        // GERP >= 1.5
        // PhyloP >= 1.5
        // LIST-S2 >= 0.85
        // PPH2 >=0.5
    }


    SNPAnnotation(short chr, int pos, char refBase, char altBase, String transcriptName, char majorBase,
                  String ancestral, double maf, double[] aaf, String daf, String[] dafs, Region region,
                  String variant_type, String derived_SIFT, String gerp, String recombinationRate,
                  String effect_snpEff, String impact_snpEff, String effect_vep, String impact_vep, String gerp16way,
                  String aaPos, String aaSize, String list_s2, String phylopP_refMask, String alleleAgeJ, String alleleAgeM
            , String alleleAgeR, String derived_PolyPhen2_prediction, String derived_PolyPhen2_class,
                  String derived_PolyPhen2_prob){
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
        this.gerp=gerp;
        this.recombinationRate=recombinationRate;
        this.effect_vep=effect_vep;
        this.impact_vep=impact_vep;
        this.effect_snpEff=effect_snpEff;
        this.impact_snpEff=impact_snpEff;

        this.gerp16way = gerp16way;
        this.aaPos= aaPos.equals("NA") ? -1 : Integer.parseInt(aaPos);
        this.aaSize= aaSize.equals("NA") ? -1 : Integer.parseInt(aaSize);
        this.list_s2=list_s2.equals("NA") ? -1 : Double.parseDouble(list_s2);
        this.phylopP_refMask=phylopP_refMask;
        this.alleleAgeJ=alleleAgeJ.equals("NA") ? -1 : Double.parseDouble(alleleAgeJ);
        this.alleleAgeM=alleleAgeM.equals("NA") ? -1 : Double.parseDouble(alleleAgeM);
        this.alleleAgeR=alleleAgeR.equals("NA") ? -1 : Double.parseDouble(alleleAgeR);
        this.derived_PolyPhen2_prediction= derived_PolyPhen2_prediction;
        this.derived_PolyPhen2_class=derived_PolyPhen2_class;
        this.derived_PolyPhen2_prob=derived_PolyPhen2_prob.equals("NA") ? -1 : Double.parseDouble(derived_PolyPhen2_prob);
    }

    public double getMaf() {
        return maf;
    }

    public String getDerived_SIFT() {
        return derived_SIFT;
    }

    public String getGerp16way(){
        return gerp16way;
    }

    public double getList_s2() {
        return list_s2;
    }

    public double getDerived_PolyPhen2_prob() {
        return derived_PolyPhen2_prob;
    }

    public String getPhylopP_refMask() {
        return phylopP_refMask;
    }

    public int getAaPos() {
        return aaPos;
    }

    public int getAaSize() {
        return aaSize;
    }

    public double getAlleleAgeJ() {
        return alleleAgeJ;
    }

    public double getAlleleAgeM() {
        return alleleAgeM;
    }

    public double getAlleleAgeR() {
        return alleleAgeR;
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
        if(this.getGerp16way().equals("NA")) return false;
        double gerp=Double.parseDouble(this.getGerp16way());
        return !(gerp <= 1);
//        if (this.getVariant_type().equals("STOP-GAIN") || this.getVariant_type().equals("START-LOST")) return true;
//        return false;
    }

    public boolean isDeleterious(MethodCallDeleterious methodCallDeleterious){
        if (!hasAncestral()) return false;
        double sift, gerp16wayRefMask, phylop, listS2, pph2;
        double gerpThreshold=1.5;
        double siftThreshold=0.05;
        double phylopThreshold=1.5;
        double listS2Threshold=0.85;
        double pph2Threshold=0.5;
        switch (methodCallDeleterious){
            case SIFT:
                if(!isNonSyn()) return false;
                if(this.getDerived_SIFT().equals("NA")) return false;
                sift=Double.parseDouble(this.getDerived_SIFT());
                if (sift >= siftThreshold) return false;
                return true;
            case GERP:
                if(!isNonSyn()) return false;
                if(this.getGerp16way().equals("NA")) return false;
                gerp16wayRefMask=Double.parseDouble(this.getGerp16way());
                return !(gerp16wayRefMask < gerpThreshold);
            case SIFT_GERP:
                if(!isNonSyn()) return false;
                if(this.getDerived_SIFT().equals("NA")) return false;
                sift=Double.parseDouble(this.getDerived_SIFT());
                if (sift >= siftThreshold) return false;
                if(this.getGerp16way().equals("NA")) return false;
                gerp16wayRefMask=Double.parseDouble(this.getGerp16way());
                return !(gerp16wayRefMask < gerpThreshold);
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
            case PHYLOP:
                if(!isNonSyn()) return false;
                if(this.getPhylopP_refMask().equals("NA")) return false;
                phylop=Double.parseDouble(this.getPhylopP_refMask());
                return !(phylop < phylopThreshold);
            case LIST_S2:
                if(!isNonSyn()) return false;
                if(this.getList_s2() < 0) return false;
                listS2=this.getList_s2();
                return !(listS2 < listS2Threshold);
            case PPH2:
                if(!isNonSyn()) return false;
                if(this.getDerived_PolyPhen2_prob() < 0) return false;
                pph2=this.getDerived_PolyPhen2_prob();
                return !(pph2 < pph2Threshold);
            case GERP_PPH2:
                if(!isNonSyn()) return false;
                if(this.getGerp16way().equals("NA")) return false;
                gerp16wayRefMask=Double.parseDouble(this.getGerp16way());
                if (gerp16wayRefMask < gerpThreshold) return false;
                if(this.getDerived_PolyPhen2_prob() < 0) return false;
                pph2=this.getDerived_PolyPhen2_prob();
                return !(pph2 < pph2Threshold);
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
