package daxing.load.ancestralSite;

import daxing.common.IOTool;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.util.*;

public class ChrSNPAnnoDB {

    private List<SNPAnnotation> snpAnnotationList;

    public ChrSNPAnnoDB(File exonSNPAnnoFile){
        this.snpAnnotationList =this.getGeneSNPAnno(exonSNPAnnoFile);
        Collections.sort(snpAnnotationList);
    }

    private List<SNPAnnotation> getGeneSNPAnno(File exonSNPAnnoFile){
        List<SNPAnnotation> geneAnno=new ArrayList<>();
        try (BufferedReader br = IOTool.getReader(exonSNPAnnoFile)) {
            br.readLine();
            String line;
            while ((line=br.readLine())!=null){
                geneAnno.add(getSNPAnnotation(line));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return geneAnno;
    }

    private SNPAnnotation getSNPAnnotation(String line){
        List<String> temp= PStringUtils.fastSplit(line);
        short chr=Short.parseShort(temp.get(1));
        int pos=Integer.parseInt(temp.get(2));
        char refBase=temp.get(3).charAt(0);
        char altBase=temp.get(4).charAt(0);
        char majorBase=temp.get(5).charAt(0);
        double[] aaf=new double[2];
        aaf[0]=Double.parseDouble(temp.get(8));
        aaf[1]=Double.parseDouble(temp.get(9));
        String geneName=temp.get(10).substring(0,18);
        SNPAnnotation.Region region=SNPAnnotation.Region.valueOf(temp.get(11));
        String variant_type=temp.get(12);
        String alt_SIFT=temp.get(13);
        String gerp=temp.get(18);
        String daf=temp.get(15);
        double maf=Double.parseDouble(temp.get(7));
        String[] dafs=new String[2];
        dafs[0]=temp.get(16);
        dafs[1]=temp.get(17);
        String ancestral=temp.get(14);
        String recombinationRate=null;
        return new SNPAnnotation(chr, pos, refBase, altBase, geneName, majorBase, ancestral, maf
                , aaf, daf, dafs, region, variant_type, alt_SIFT, gerp, recombinationRate);
    }

    /**
     *
     * @param chr
     * @param pos
     * @return negative value, not in DB
     */
    public int binarySearch(int chr, int pos){
        return Collections.binarySearch(this.snpAnnotationList, new ChrPos((short) chr, pos));
    }

    public SNPAnnotation getSNPAnnotation(int index){
        return this.snpAnnotationList.get(index);
    }

    public String getGeneName(int chr, int pos){
        int snpIndex= this.binarySearch(chr, pos);
        if (snpIndex < 0){
            System.out.println("error");
            System.exit(1);
        }
        return this.getSNPAnnotation(snpIndex).getSNPInfo();
    }

    public String[] getGeneName(){
        List<SNPAnnotation> snpAnnotations=this.snpAnnotationList;
        Set<String> geneNames=new HashSet<>();
        String geneName;
        for (int i = 0; i < snpAnnotations.size(); i++) {
            geneName=snpAnnotations.get(i).getSNPInfo();
            geneNames.add(geneName);
        }
        return geneNames.stream().sorted().toArray(String[]::new);
    }

    public SNPAnnotation getSNP(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        return this.snpAnnotationList.get(index);
    }

    public boolean contain(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        if (index < 0) return false;
        return true;
    }

    /**
     *
     * @param chr
     * @param pos
     * @return
     */
    public String getVariantType(int chr, int pos){
        int index=this.binarySearch(chr, pos);
        if (index < 0) return "INTERGENIC";
        return this.getSNP(chr, pos).variant_type;
    }

    public boolean hasAncestral(int chr, int pos){
        return this.getSNP(chr, pos).hasAncestral();
    }

    public boolean isRefAlleleAncestral(int chr, int pos){
        return this.getSNP(chr, pos).isRefAlleleAncestral();
    }

    public boolean isSyn(int chr, int pos){
        return this.getSNP(chr, pos).isSyn();
    }

    public boolean isNonsyn(int chr, int pos){
        return this.getSNP(chr, pos).isNonSyn();
    }

    public boolean isDeleterious(int chr, int pos){
        return this.getSNP(chr, pos).isDeleterious();
    }
}
