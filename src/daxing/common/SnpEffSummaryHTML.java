package daxing.common;

import gnu.trove.list.array.TIntArrayList;
import org.apache.commons.lang.StringUtils;
import pgl.infra.utils.IOUtils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.function.Predicate;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * SnpEff version: SnpEff 4.3t (build 2017-11-24 10:18), by Pablo Cingolani
 * @author Daxing Xu
 */
public class SnpEffSummaryHTML {

    private int numberOfVariants;
    private int numberOfEffects;
    private int chrLen;
    private int[] numberOfEffectsByImpact;   //HIGH, LOW, MODERATE, MODIFIER
    private int[] numberOfEffectsByFunctionalClass; //MISSENSE, NONSENSE, SILENT
    private int[] numberOfEffectsByType;  //3_prime_UTR_variant, 5_prime_UTR_premature_start_codon_gain_variant, ...
    private int[] numberOfEffectsByRegion;  //DOWNSTREAM, EXON, INTERGENIC, INTRON, ...

    public SnpEffSummaryHTML(File inputFileOfSummaryHTML){
        this.parseSummaryHTML(inputFileOfSummaryHTML);
    }

    private void parseSummaryHTML(File inputFileOfSummaryHTML){
        try (BufferedReader br = IOUtils.getTextReader(inputFileOfSummaryHTML.getAbsolutePath())) {
            String line;
            boolean first=false;
            TIntArrayList numberOfEffectsByImpactList=new TIntArrayList();
            TIntArrayList numberOfEffectsByFunctionalClassList=new TIntArrayList();
            TIntArrayList numberOfEffectsByTypeList=new TIntArrayList();
            TIntArrayList numberOfEffectsByRegionList=new TIntArrayList();
            List<String> texts;
            while ((line = br.readLine()) != null) {
                if (line.contains("Number of variants (before filter)")){
                    line=br.readLine();
                    texts=StringTool.parseLineInHtmlFormat(line, "body");
                    numberOfVariants=StringTool.parseInt(texts.get(0));
                    first=true;
                }
                if (line.contains("Number of effects") && first){
                    line=br.readLine();
                    texts=StringTool.parseLineInHtmlFormat(line, "body");
                    numberOfEffects=StringTool.parseInt(texts.get(0));
                    first=false;
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();br.readLine();br.readLine();
                    br.readLine();
                    line=br.readLine();
                    texts=StringTool.parseLineInHtmlFormat(line, "body");
                    chrLen=StringTool.parseInt(texts.get(0));
                }
                if (line.contains("HIGH")){
                    numberOfEffectsByImpactList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByImpactList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("MISSENSE")){
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByFunctionalClassList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("3_prime_UTR_variant")) {
                    numberOfEffectsByTypeList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByTypeList.add(this.extractKeyWord_5(br));
                }
                if (line.contains("DOWNSTREAM")){
                    numberOfEffectsByRegionList.add(this.extractKeyWord_1(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                    numberOfEffectsByRegionList.add(this.extractKeyWord_5(br));
                }
            }
            numberOfEffectsByImpact=numberOfEffectsByImpactList.toArray();
            numberOfEffectsByFunctionalClass=numberOfEffectsByFunctionalClassList.toArray();
            numberOfEffectsByType=numberOfEffectsByTypeList.toArray();
            numberOfEffectsByRegion=numberOfEffectsByRegionList.toArray();
        }catch (Exception e){
            e.printStackTrace();
        }
    }

    private int extractKeyWord_1(BufferedReader br){
        int num=Integer.MIN_VALUE;
        try {
            String line;
            List<String> temp=new ArrayList<>();
            br.readLine();
            line=br.readLine();
            List<String> texts=StringTool.parseLineInHtmlFormat(line, "body");
            return StringTool.parseInt(texts.get(0));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return num;
    }

    private int extractKeyWord_5(BufferedReader br){
        int num=Integer.MIN_VALUE;
        try {
            String line;
            List<String> temp=new ArrayList<>();
            br.readLine();
            br.readLine();
            br.readLine();
            br.readLine();
            br.readLine();
            line=br.readLine();
            List<String> texts=StringTool.parseLineInHtmlFormat(line, "body");
            return StringTool.parseInt(texts.get(0));
        } catch (IOException e) {
            e.printStackTrace();
        }
        return num;
    }

    public int getNumberOfVariants() {
        return numberOfVariants;
    }

    public int getNumberOfEffects() {
        return numberOfEffects;
    }

    public int getChrLen() {
        return chrLen;
    }

    public int[] getNumberOfEffectsByImpact() {
        return numberOfEffectsByImpact;
    }

    public int[] getNumberOfEffectsByFunctionalClass() {
        return numberOfEffectsByFunctionalClass;
    }

    public int[] getNumberOfEffectsByType() {
        return numberOfEffectsByType;
    }

    public int[] getNumberOfEffectsByRegion() {
        return numberOfEffectsByRegion;
    }

    public double getVariantsRate(){
        return (double) chrLen/numberOfVariants;
    }

    public double getMissenseSilentRatio(){
        return (double) numberOfEffectsByFunctionalClass[0]/numberOfEffectsByFunctionalClass[2];
    }

    public double[] getImpactPercent(){
        double sum= Arrays.stream(this.numberOfEffectsByImpact).sum();
        return Arrays.stream(this.numberOfEffectsByImpact).mapToDouble(e->e/sum).toArray();
    }

    public double[] getFunctionPercent(){
        double sum= Arrays.stream(numberOfEffectsByFunctionalClass).sum();
        return Arrays.stream(this.numberOfEffectsByFunctionalClass).mapToDouble(e->e/sum).toArray();
    }

    public double[] getEffectsTypePercent(){
        double sum= Arrays.stream(numberOfEffectsByType).sum();
        return Arrays.stream(this.numberOfEffectsByType).mapToDouble(e->e/sum).toArray();
    }

    public double[] getEffectsRegionPercent(){
        double sum= Arrays.stream(numberOfEffectsByRegion).sum();
        return Arrays.stream(this.numberOfEffectsByRegion).mapToDouble(e->e/sum).toArray();
    }

    /**
     * 解析snpEff软件产生的summary.html文件，统计Effects by type and region
     * @param summaryHtmlFileDir 包含HTML文件的目录
     * @param outFileDirOfEffectsByTypeAndRegion
     */
    public static void getA_B_D_lineageTypeAndRegion(String summaryHtmlFileDir,
                                                     String outFileDirOfEffectsByTypeAndRegion){
        File[] input=IOUtils.listRecursiveFiles(new File(summaryHtmlFileDir));
        Predicate<File> p=File::isHidden;
        File[] files=Arrays.stream(input).filter(p.negate()).filter(f->f.getName().endsWith("html")).toArray(File[]::new);
        List<SnpEffSummaryHTML> snpEffSummaryHTMLS=Arrays.stream(files).map(SnpEffSummaryHTML::new).collect(Collectors.toList());
        int[] a_LineageIndex= IntStream.iterate(1, n->n+6).map(n->n-1).limit(7).toArray();
        int[] b_LineageIndex=IntStream.iterate(3, n->n+6).map(n->n-1).limit(7).toArray();
        int[] d_LineageIndex=IntStream.iterate(5, n->n+6).map(n->n-1).limit(7).toArray();
        int[] typeA, typeB, typeD;
        int[] typeAll_A=snpEffSummaryHTMLS.get(a_LineageIndex[0]).getNumberOfEffectsByType();
        int[] typeAll_B=snpEffSummaryHTMLS.get(b_LineageIndex[0]).getNumberOfEffectsByType();
        int[] typeAll_D=snpEffSummaryHTMLS.get(d_LineageIndex[0]).getNumberOfEffectsByType();
        double[] percentTypeAll_A, percentTypeAll_B, percentTypeAll_D;
        for (int i = 1; i < a_LineageIndex.length; i=i+1) {
            typeA=snpEffSummaryHTMLS.get(a_LineageIndex[i]).getNumberOfEffectsByType();
            typeAll_A=ArrayTool.add(typeAll_A, typeA);
        }
        for (int i = 1; i < b_LineageIndex.length; i=i+1) {
            typeB=snpEffSummaryHTMLS.get(b_LineageIndex[i]).getNumberOfEffectsByType();
            typeAll_B=ArrayTool.add(typeAll_B, typeB);
        }
        for (int i = 1; i < d_LineageIndex.length; i=i+1) {
            typeD=snpEffSummaryHTMLS.get(d_LineageIndex[i]).getNumberOfEffectsByType();
            typeAll_D=ArrayTool.add(typeAll_D, typeD);
        }
        int[] regionA, regionB, regionD;
        int[] regionAll_A=snpEffSummaryHTMLS.get(a_LineageIndex[0]).getNumberOfEffectsByRegion();
        int[] regionAll_B=snpEffSummaryHTMLS.get(b_LineageIndex[0]).getNumberOfEffectsByRegion();
        int[] regionAll_D=snpEffSummaryHTMLS.get(d_LineageIndex[0]).getNumberOfEffectsByRegion();
        double[] percentRegionAll_A, percentRegionAll_B, percentRegionAll_D;
        for (int i = 1; i < a_LineageIndex.length; i++) {
            regionA=snpEffSummaryHTMLS.get(a_LineageIndex[i]).getNumberOfEffectsByRegion();
            regionAll_A=ArrayTool.add(regionAll_A, regionA);
        }
        for (int i = 1; i < b_LineageIndex.length; i++) {
            regionB=snpEffSummaryHTMLS.get(b_LineageIndex[i]).getNumberOfEffectsByRegion();
            regionAll_B=ArrayTool.add(regionAll_B, regionB);
        }
        for (int i = 1; i < d_LineageIndex.length; i++) {
            regionD=snpEffSummaryHTMLS.get(d_LineageIndex[0]).getNumberOfEffectsByRegion();
            regionAll_D=ArrayTool.add(regionAll_D, regionD);
        }
        percentTypeAll_A=ArrayTool.getElementPercent(typeAll_A);
        percentTypeAll_B=ArrayTool.getElementPercent(typeAll_B);
        percentTypeAll_D=ArrayTool.getElementPercent(typeAll_D);
        percentRegionAll_A=ArrayTool.getElementPercent(regionAll_A);
        percentRegionAll_B=ArrayTool.getElementPercent(regionAll_B);
        percentRegionAll_D=ArrayTool.getElementPercent(regionAll_D);
        String header="A_count\tA_percent\tB_count\tB_percent\tD_count\tD_percent";
        String columnForType="3_prime_UTR_variant\t5_prime_UTR_premature_start_codon_gain_variant\t5_prime_UTR_variant\tdownstream_gene_variant\tinitiator_codon_variant\tintergenic_region\tintron_variant\tmissense_variant\tsplice_acceptor_variant\tsplice_donor_variant\tsplice_region_variant\tstart_lost\tstop_gained\tstop_lost\tstop_retained_variant\tsynonymous_variant\tupstream_gene_variant";
        String columnForRegion="DOWNSTREAM\tEXON\tINTERGENIC\tINTRON\tSPLICE_SITE_ACCEPTOR\tSPLICE_SITE_DONOR\tSPLICE_SITE_REGION\tUPSTREAM\tUTR_3_PRIME\tUTR_5_PRIME";
        String[] typeColumn= StringUtils.split(columnForType, "\t");
        String[] regionColumn=StringUtils.split(columnForRegion, "\t");
        try(BufferedWriter bwType=IOUtils.getTextWriter(new File(outFileDirOfEffectsByTypeAndRegion, "type.txt").getAbsolutePath());
            BufferedWriter bwRegion=IOUtils.getTextWriter(new File(outFileDirOfEffectsByTypeAndRegion, "region.txt").getAbsolutePath())){
            bwType.write("Type"+"\t"+header);
            bwRegion.write("Region"+"\t"+header);
            bwType.newLine();
            bwRegion.newLine();
            StringBuilder sb;
            for (int i = 0; i < typeAll_A.length; i++) {
                sb=new StringBuilder();
                sb.append(typeColumn[i]).append("\t").append(typeAll_A[i]).append("\t").append(percentTypeAll_A[i]).append("\t")
                        .append(typeAll_B[i]).append("\t").append(percentTypeAll_B[i]).append("\t").append(typeAll_D[i]).append("\t")
                        .append(percentTypeAll_D[i]);
                bwType.write(sb.toString());
                bwType.newLine();
            }
            bwType.flush();
            for (int i = 0; i < regionAll_A.length; i++) {
                sb=new StringBuilder();
                sb.append(regionColumn[i]).append("\t").append(regionAll_A[i]).append("\t").append(percentRegionAll_A[i]).append("\t")
                        .append(regionAll_B[i]).append("\t").append(percentRegionAll_B[i]).append("\t").append(regionAll_D[i]).append("\t")
                        .append(percentRegionAll_D[i]);
                bwRegion.write(sb.toString());
                bwRegion.newLine();
            }
            bwRegion.flush();
        }catch (Exception e){
            e.printStackTrace();
        }
    }


}
