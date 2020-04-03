package daxing.load;

import java.util.Arrays;

public class SNPAnnotation {
    int chr;
    int pos;
    String[] allels=new String[2];
    int refAlleleIndex;
    int majorAlleleIndex;
    //-1 表示ancestral为NA或vmapII bialleles不包含ancestra allele, O为alleles[0], 1为alleles[1]
    int ancestralIndex;
    String transcriptName;

    public enum Region{
        UTR_5, CDS, UTR_3;
    }

    public enum Variant_type{
        NONCODING, NONSYNONYMOUS, START_LOST, STOP_GAIN, STOP_LOSS, SYNONYMOUS;
    }

    Region region;
    Variant_type variant_type;

    // -1表示NA
    double alt_SIFT;
    double ref_SIFT;
    double derived_SIFT;
    double gerp;
    double recombinationRate;

    SNPAnnotation(int chr, int pos, String[] allels, String refAllele, String majorAllele, String ancestral,
                  String transcriptName, String region, String variant_type, String alt_SIFT, String ref_SIFT,
                  String derived_SIFT, String gerp, String recombinationRate){
        this.chr=chr;
        this.pos=pos;
        Arrays.sort(allels);
        this.allels=allels;
        refAlleleIndex=Arrays.binarySearch(allels, refAllele);
        majorAlleleIndex=Arrays.binarySearch(allels, majorAllele);
        if (ancestral.equals("NA")){
            ancestralIndex=-1;
        }else {
            int index=Arrays.binarySearch(allels, ancestral);
            ancestralIndex=index < 0 ? -1 : index;
        }
        this.transcriptName=transcriptName;
        this.region=Region.valueOf(region);
        this.variant_type=Variant_type.valueOf(variant_type);
        this.alt_SIFT=alt_SIFT.equals("NA") ? -1 : Double.parseDouble(alt_SIFT);
        this.ref_SIFT=ref_SIFT.equals("NA") ? -1 : Double.parseDouble(ref_SIFT);
        this.derived_SIFT=derived_SIFT.equals("NA") ? -1 : Double.parseDouble(derived_SIFT);
        this.gerp=gerp.equals("NA") ? Double.NaN : Double.parseDouble(gerp);
        this.recombinationRate=recombinationRate.equals("NA") ? -1 : Double.parseDouble(recombinationRate);
    }

    public String getID(){
        StringBuilder sb=new StringBuilder();
        sb.append(chr).append("-").append(pos);
        return sb.toString();
    }

    public String getRefAllele(){
        return allels[refAlleleIndex];
    }

    public String getAltAllele(){
        int altAlleleIndex=refAlleleIndex==0 ? 1 : 0;
        return allels[altAlleleIndex];
    }

    public String getRegion(){
        return region.name();
    }

    public String getVariant_type(){
        return variant_type.name();
    }

    public boolean hasGerp(){
        return Double.isNaN(gerp) ? false : true;
    }

    public double getGerp(){
        return gerp;
    }

}
