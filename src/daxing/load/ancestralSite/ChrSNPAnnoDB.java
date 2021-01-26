package daxing.load.ancestralSite;

import daxing.common.IOTool;
import pgl.infra.pos.ChrPos;
import pgl.infra.utils.PStringUtils;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Predicate;
import java.util.stream.Collectors;

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
        double maf=Double.parseDouble(temp.get(7));
        aaf[0]=Double.parseDouble(temp.get(8));
        aaf[1]=Double.parseDouble(temp.get(9));
        String geneName=temp.get(10).substring(0,18);
        SNPAnnotation.Region region=SNPAnnotation.Region.valueOf(temp.get(11));
        String variant_type=temp.get(12);
        String ancestral=temp.get(15);
        String derived_SIFT=temp.get(16);
        String daf=temp.get(17);
        String[] dafs=new String[2];
        dafs[0]=temp.get(18);
        dafs[1]=temp.get(19);
        String gerp=temp.get(20);
        String recombinationRate=null;
        return new SNPAnnotation(chr, pos, refBase, altBase, geneName, majorBase, ancestral, maf
                , aaf, daf, dafs, region, variant_type, derived_SIFT, gerp, recombinationRate);
    }

    public int getSize(){
        return this.snpAnnotationList.size();
    }

    public List<SNPAnnotation> getSnpAnnotationList() {
        return snpAnnotationList;
    }

    /**
     *
     * @param chr chr
     * @param pos pos
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
        for (SNPAnnotation snpAnnotation : snpAnnotations) {
            geneName = snpAnnotation.getSNPInfo();
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
        return index >= 0;
    }

    /**
     *
     * @param chr chr
     * @param pos pos
     * @return variantType
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

    public int getDelNum(){
        List<SNPAnnotation> snpAnnotations=this.snpAnnotationList;
        Predicate<SNPAnnotation> p = SNPAnnotation::isDeleterious;
        return (int) snpAnnotations.stream().filter(p).count();
    }

    public int getNum(Predicate<SNPAnnotation> p){
        return (int) this.snpAnnotationList.stream().filter(p).count();
    }

    public ChrSNPAnnoDB filter(Predicate<SNPAnnotation> p){
        this.snpAnnotationList= this.snpAnnotationList.stream().filter(p).collect(Collectors.toList());
        return this;
    }

    public void write(String outFile){
        StringBuilder sb=new StringBuilder();
        SNPAnnotation snpAnnotation;
        try (BufferedWriter bw = IOTool.getWriter(outFile)) {
            bw.write("Chr\tPos");
            bw.newLine();
            for (SNPAnnotation annotation : snpAnnotationList) {
                snpAnnotation = annotation;
                sb.setLength(0);
                sb.append(snpAnnotation.getChromosome()).append("\t").append(snpAnnotation.getPos());
                bw.write(sb.toString());
                bw.newLine();
            }
            bw.flush();
        } catch (IOException e) {
            e.printStackTrace();
        }

    }
}
